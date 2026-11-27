#include <assert.h>
#include <gkyl_vlasov_priv.h>
#include <gkyl_wv_euler.h>

// Temperature source of an Euler fluid: the moments of the BGK-form source of a kinetic species,
// S = nu_Q (f_M(n, u, T_Q) - f), which leave the density and the momentum unchanged and relax the
// pressure, d_t p = nu_Q (n T_Q - p), i.e. d_t E = (nu_Q T_Q/m rho - nu_Q p)/(gas_gamma - 1).
// The explicit form is in the right-hand side of every stage; the implicit form is one
// backward-Euler step in each stage, at the Gauss-Legendre nodes (the fluid part of the implicit
// source of a PKPM species, with the thermal energy p/(gas_gamma - 1)).

// Context of the nodal evaluation of the coefficients nu_Q T_Q/m and nu_Q.
struct vm_fluid_temp_source_ctx {
  const struct gkyl_vlasov_source_bgk *info; // source input
  double mass; // mass of the species
};

static void
eval_temp_source_coeffs(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct vm_fluid_temp_source_ctx *c = ctx;
  double nu = 0.0, temp = 0.0;
  c->info->rate_profile(t, xn, &nu, c->info->rate_profile_ctx);
  c->info->temp_shape(t, xn, &temp, c->info->temp_shape_ctx);
  fout[0] = nu * temp / c->mass;
  fout[1] = nu;
}

// Explicit source: d_t E = (nu_Q T_Q/m rho - nu_Q p)/(gas_gamma - 1) of the physical state and
// pressure, and its rate in the time step.
static void
fluid_temp_source_rhs_explicit(
  gkyl_vlasov_app *app, struct vm_fluid_species *f, const struct gkyl_array *fluid,
  struct gkyl_array *rhs
)
{
  struct vm_fluid_temp_source *tsrc = &f->tsrc;
  int nb = app->basis.num_basis;
  gkyl_dg_mul_op_range(&app->basis, 0, tsrc->buf, 0, tsrc->nu_temp, 0, fluid, &app->local);
  gkyl_array_accumulate_offset_range(rhs, tsrc->inv_gamma_m1, tsrc->buf, 4 * nb, &app->local);
  gkyl_dg_mul_op_range(&app->basis, 0, tsrc->buf, 0, tsrc->nu, 0, f->p, &app->local);
  gkyl_array_accumulate_offset_range(rhs, -tsrc->inv_gamma_m1, tsrc->buf, 4 * nb, &app->local);
  gkyl_array_accumulate_range(f->cflrate, 1.0, tsrc->nu_cfl, &app->local);
}

// No source, or the implicit one: nothing enters the right-hand side.
static void
fluid_temp_source_rhs_none(
  gkyl_vlasov_app *app, struct vm_fluid_species *f, const struct gkyl_array *fluid,
  struct gkyl_array *rhs
)
{
}

// Implicit source: the rate and the temperature over the mass at the Gauss-Legendre nodes, and
// the updater of the backward-Euler step.
static void
fluid_temp_source_init_implicit(
  struct gkyl_vlasov_app *app, struct vm_fluid_species *f, struct vm_fluid_temp_source *tsrc
)
{
  const struct gkyl_vlasov_source_bgk *info = &f->info.source_bgk;
  int nb = app->basis.num_basis;
  tsrc->rate =
    vm_pkpm_species_coll_profile_new(app, info->rate_profile, info->rate_profile_ctx, 1.0, true);
  tsrc->vtsq = vm_pkpm_species_coll_profile_new(
    app, info->temp_shape, info->temp_shape_ctx, 1.0 / f->mass, true
  );
  tsrc->shape = mkarr(false, nb, app->local_ext.volume);
  tsrc->temp_ratio = mkarr(false, nb, app->local_ext.volume);
  tsrc->up_temp = gkyl_pkpm_temp_source_inew(&(struct gkyl_pkpm_temp_source_inp){
    .conf_basis = &app->basis,
    .energy_fac = tsrc->inv_gamma_m1,
    .use_gpu = false,
  });
}

void
vm_fluid_species_temp_source_init(
  struct gkyl_vlasov_app *app, struct vm_fluid_species *f, struct vm_fluid_temp_source *tsrc
)
{
  const struct gkyl_vlasov_source_bgk *info = &f->info.source_bgk;
  *tsrc = (struct vm_fluid_temp_source){.is_on = false, .rhs_func = fluid_temp_source_rhs_none};
  if (info->source_bgk_id == GKYL_NO_SOURCE_BGK) {
    return;
  }
  if (info->source_bgk_id != GKYL_SOURCE_BGK_HEATING || f->eqn_type != GKYL_EQN_EULER ||
      !info->rate_profile || !info->temp_shape) {
    gkyl_exit(
      "vm_fluid_species: a fluid temperature source (GKYL_SOURCE_BGK_HEATING, with rate_profile and temp_shape) needs an Euler fluid."
    );
  }
  tsrc->is_on = true;
  tsrc->is_implicit = info->is_implicit;
  tsrc->inv_gamma_m1 = 1.0 / (gkyl_wv_euler_gas_gamma(f->equation) - 1.0);
  if (tsrc->is_implicit) {
    // polynomial order 1 on the CPU: checked by the app, with the collisions of fluid species
    fluid_temp_source_init_implicit(app, f, tsrc);
    return;
  }
  tsrc->rhs_func = fluid_temp_source_rhs_explicit;

  int nb = app->basis.num_basis;
  struct gkyl_array *coeffs_ho = mkarr(false, 2 * nb, app->local_ext.volume);
  struct vm_fluid_temp_source_ctx ctx = {.info = info, .mass = f->mass};
  // At the nodes (continuous and non-negative where nu_Q is), at the physical coordinates.
  struct gkyl_eval_on_nodes *ev = gkyl_eval_on_nodes_inew(&(struct gkyl_eval_on_nodes_inp){
    .grid = &app->grid,
    .basis = &app->basis,
    .num_ret_vals = 2,
    .eval = eval_temp_source_coeffs,
    .ctx = &ctx,
    .c2p_func = vm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });
  gkyl_eval_on_nodes_advance(ev, 0.0, &app->local_ext, coeffs_ho);
  gkyl_eval_on_nodes_release(ev);

  struct gkyl_array *nu_temp_ho = mkarr(false, nb, app->local_ext.volume);
  struct gkyl_array *nu_ho = mkarr(false, nb, app->local_ext.volume);
  gkyl_array_set_offset(nu_temp_ho, 1.0, coeffs_ho, 0);
  gkyl_array_set_offset(nu_ho, 1.0, coeffs_ho, nb);
  tsrc->nu_temp = mkarr(app->use_gpu, nb, app->local_ext.volume);
  tsrc->nu = mkarr(app->use_gpu, nb, app->local_ext.volume);
  gkyl_array_copy(tsrc->nu_temp, nu_temp_ho);
  gkyl_array_copy(tsrc->nu, nu_ho);
  tsrc->nu_cfl = mkarr(app->use_gpu, 1, app->local_ext.volume);
  gkyl_dg_calc_average_range(&app->basis, 0, tsrc->nu_cfl, 0, tsrc->nu, app->local);
  tsrc->buf = mkarr(app->use_gpu, nb, app->local_ext.volume);
  gkyl_array_release(coeffs_ho);
  gkyl_array_release(nu_temp_ho);
  gkyl_array_release(nu_ho);
}

void
vm_fluid_species_temp_source_rhs(
  gkyl_vlasov_app *app, struct vm_fluid_species *f, const struct gkyl_array *fluid,
  struct gkyl_array *rhs
)
{
  f->tsrc.rhs_func(app, f, fluid, rhs);
}

void
vm_fluid_species_temp_source_implicit_advance(
  gkyl_vlasov_app *app, struct vm_fluid_species *f, double dt, struct gkyl_array *fluid
)
{
  struct vm_fluid_temp_source *tsrc = &f->tsrc;
  // the state may carry the position-map Jacobian: the step is homogeneous in it
  gkyl_pkpm_temp_source_implicit_advance(
    tsrc->up_temp, &app->local, dt, tsrc->rate, tsrc->vtsq, fluid, tsrc->shape, tsrc->temp_ratio
  );
  // The step keeps the thermal energy positive at the nodes, not at the vertices, so the vertex
  // limiter and the boundary conditions follow it.
  vm_fluid_species_apply_positivity(app, f, fluid);
  vm_fluid_species_apply_bc(app, f, fluid);
}

void
vm_fluid_species_temp_source_release(
  const struct gkyl_vlasov_app *app, const struct vm_fluid_temp_source *tsrc
)
{
  if (!tsrc->is_on) {
    return;
  }
  if (tsrc->is_implicit) {
    gkyl_array_release(tsrc->rate);
    gkyl_array_release(tsrc->vtsq);
    gkyl_array_release(tsrc->shape);
    gkyl_array_release(tsrc->temp_ratio);
    gkyl_pkpm_temp_source_release(tsrc->up_temp);
    return;
  }
  gkyl_array_release(tsrc->nu_temp);
  gkyl_array_release(tsrc->nu);
  gkyl_array_release(tsrc->nu_cfl);
  gkyl_array_release(tsrc->buf);
}
