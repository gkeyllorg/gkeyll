#include <assert.h>
#include <gkyl_pkpm_priv.h>

// Evaluate a p=1 configuration-space expansion at the computational coordinate xn.
static double
eval_conf_field(
  const struct pkpm_bgk_eval_ctx *c, const struct gkyl_array *arr, int comp, const double *xn
)
{
  int idx[GKYL_MAX_CDIM];
  double xc[GKYL_MAX_CDIM], eta[GKYL_MAX_CDIM];
  gkyl_rect_grid_coord_idx(c->conf_grid, xn, idx);
  for (int d = 0; d < c->cdim; ++d) {
    idx[d] = GKYL_MAX2(GKYL_MIN2(idx[d], c->local_ext->upper[d]), c->local_ext->lower[d]);
  }
  gkyl_rect_grid_cell_center(c->conf_grid, idx, xc);
  for (int d = 0; d < c->cdim; ++d) {
    eta[d] = 2.0 * (xn[d] - xc[d]) / c->conf_grid->dx[d];
  }
  const double *coeff = gkyl_array_cfetch(arr, gkyl_range_idx(c->local_ext, idx));
  return c->conf_basis->eval_expand(eta, &coeff[comp * c->conf_basis->num_basis]);
}

// The BGK target in the normalized variables: a Maxwellian of variance r = T_0/(m xi^2)
// in w (r = 1 when relaxing to the local temperature) and its isotropic partner G = r F.
static void
eval_bgk_target(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct pkpm_bgk_eval_ctx *c = ctx;
  double w = xn[c->cdim];
  double r = 1.0;
  if (c->fixed_temp_relax) {
    double xisq = eval_conf_field(c, c->xi, 1, xn),
           vtsq0 = eval_conf_field(c, c->fixed_vtsq, 0, xn);
    r = (xisq > 0.0 && vtsq0 > 0.0) ? vtsq0 / xisq : 1.0;
  }
  double F = exp(-0.5 * w * w / r) / sqrt(2.0 * M_PI * r);
  fout[0] = F;
  fout[1] = r * F;
}

// nu (f_target - f) at the nodes: all blocks relax to the corrected unit Maxwellian.
static void
pkpm_species_bgk_kinetic_rhs_nodal(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  gkyl_pkpm_bgk_collisions_advance(
    bgk->up_bgk, &app->local, &s->local_vel, &s->local, bgk->rate, app->pos_map->jacob_pos_gauss,
    bgk->f_lte, fin, rhs, s->cflrate_f
  );
}

// The implicit operator adds nothing to the right-hand side of the explicit update.
static void
pkpm_species_bgk_kinetic_rhs_none(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
}

// nu (f_target - f) in the modal form, for the fixed-temperature target and GPU runs; the target
// is gyrotropic ([F, G] only), so the first harmonic [M0_perp, N_perp] relaxes to zero (before
// the closures: with a fixed temperature the target's energy differs from the state's and the
// closures keep the constraints exact while the fluid energy carries the change).
static void
pkpm_species_bgk_kinetic_rhs_modal(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  if (app->use_gpu) {
    gkyl_array_copy(bgk->xi_host, s->pkpm_xi);
    bgk->eval_ctx.xi = bgk->xi_host;
  } else {
    bgk->eval_ctx.xi = s->pkpm_xi;
  }
  gkyl_proj_on_basis_advance(bgk->proj_target, app->tcurr, &s->local, bgk->f_target_host);
  if (app->use_gpu) {
    gkyl_array_copy(bgk->f_target, bgk->f_target_host);
  }
  gkyl_vlasov_position_map_rescale_jacobpos(
    app->pos_map, &app->basis, &s->local, bgk->f_target, bgk->f_target
  );
  pkpm_species_accumulate_conf_phase(app, s, rhs, 1.0, bgk->nu, bgk->f_target);
  pkpm_species_accumulate_conf_phase(app, s, rhs, -1.0, bgk->nu, fin);
}

// The collision frequency at the nodes and the single velocity-space target of the operator:
// the unit Maxwellian in w, with its discrete moments corrected to [1, 0, 1].
static void
pkpm_species_bgk_init_nodal(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk,
  const struct gkyl_array *nu_host
)
{
  bgk->rate = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  gkyl_pkpm_coll_nu *nu_calc = gkyl_pkpm_coll_nu_new(&app->confBasis, false, 0.0, false);
  gkyl_pkpm_coll_nu_to_nodal(nu_calc, &app->local_ext, nu_host, bgk->rate);
  gkyl_pkpm_coll_nu_release(nu_calc);

  bgk->f_lte = mkarr(false, 3, s->local_ext_vel.volume);
  gkyl_pkpm_lte *lte = gkyl_pkpm_lte_inew(&(struct gkyl_pkpm_lte_inp){
    .vel_grid = &s->grid_vel,
    .vel_range = &s->local_vel,
    .num_quad = 4,
    .eps = 1.0e-14,
    .max_iter = 100,
  });
  double moms[3] = {1.0, 0.0, 1.0};
  struct gkyl_pkpm_lte_status status = gkyl_pkpm_lte_correct(lte, moms, bgk->f_lte);
  gkyl_pkpm_lte_release(lte);
  if (!status.converged) {
    gkyl_pkpm_app_cout(
      app, stderr, "WARNING: species %s: the BGK target moments reached a relative error of %g\n",
      s->info.name, status.error
    );
  }

  bgk->up_bgk = gkyl_pkpm_bgk_collisions_new(&app->confBasis, &app->basis, false);
  bgk->is_implicit = s->info.collisions.is_implicit;
  if (bgk->is_implicit) {
    bgk->shape = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    bgk->kinetic_rhs = pkpm_species_bgk_kinetic_rhs_none;
  } else {
    bgk->kinetic_rhs = pkpm_species_bgk_kinetic_rhs_nodal;
  }
}

static void
pkpm_species_bgk_init_modal(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk,
  const struct gkyl_array *nu_host
)
{
  int nc = app->confBasis.num_basis, np = app->basis.num_basis;
  bgk->nu = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_copy(bgk->nu, nu_host);

  bgk->fixed_vtsq = mkarr(app->use_gpu, nc, app->local_ext.volume);
  gkyl_array_clear(bgk->fixed_vtsq, 0.0);
  bgk->fixed_vtsq_host = app->use_gpu ? mkarr(false, nc, app->local_ext.volume) :
                                        gkyl_array_acquire(bgk->fixed_vtsq);
  bgk->xi_host = app->use_gpu ? mkarr(false, 3 * nc, app->local_ext.volume) : 0;
  bgk->f_target = mkarr(app->use_gpu, 2 * np, s->local_ext.volume);
  bgk->f_target_host = app->use_gpu ? mkarr(false, 2 * np, s->local_ext.volume) :
                                      gkyl_array_acquire(bgk->f_target);
  gkyl_array_clear(bgk->f_target, 0.0);
  bgk->E_target = mkarr(app->use_gpu, nc, app->local_ext.volume);
  bgk->scratch = mkarr(app->use_gpu, nc, app->local_ext.volume);
  bgk->uu = mkarr(app->use_gpu, nc, app->local_ext.volume);
  bgk->eval_ctx = (struct pkpm_bgk_eval_ctx){
    .cdim = app->cdim,
    .conf_grid = &app->grid,
    .conf_basis = &app->confBasis,
    .local_ext = &app->local_ext,
    .fixed_temp_relax = bgk->fixed_temp_relax,
    .fixed_vtsq = bgk->fixed_vtsq_host,
  };
  bgk->proj_target = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &s->grid,
    .basis = &app->basis,
    .num_quad = 4,
    .num_ret_vals = 2,
    .eval = eval_bgk_target,
    .ctx = &bgk->eval_ctx,
  });
  bgk->kinetic_rhs = pkpm_species_bgk_kinetic_rhs_modal;
}

void
pkpm_species_bgk_init(gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk)
{
  // Collision frequency (projected once, like the LBO's).
  struct gkyl_array *nu_host = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &app->grid,
    .basis = &app->confBasis,
    .num_quad = app->poly_order + 1,
    .num_ret_vals = 1,
    .eval = s->info.collisions.self_nu,
    .ctx = s->info.collisions.ctx,
    .c2p_func = pkpm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local_ext, nu_host);
  gkyl_proj_on_basis_release(proj);
  double nu_frac = s->info.collisions.nuFrac > 0.0 ? s->info.collisions.nuFrac : 1.0;
  gkyl_array_scale(nu_host, nu_frac);

  bgk->fixed_temp_relax = s->info.collisions.fixed_temp_relax;
  // The fixed-temperature target and GPU runs take the explicit modal path.
  bgk->use_nodal = !bgk->fixed_temp_relax && !app->use_gpu;
  bgk->is_implicit = false;
  if (bgk->use_nodal) {
    pkpm_species_bgk_init_nodal(app, s, bgk, nu_host);
  } else {
    pkpm_species_bgk_init_modal(app, s, bgk, nu_host);
  }
  gkyl_array_release(nu_host);
}

// The fixed temperature is the initial one: xi^2 = T_0/m of the initial fluid state (call
// after the initial condition, with the PKPM variables of that state computed).
void
pkpm_species_bgk_set_fixed_temp(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk
)
{
  int nc = app->confBasis.num_basis;
  gkyl_array_set_offset(bgk->fixed_vtsq, 1.0, s->pkpm_xi, nc);
  gkyl_array_copy(bgk->fixed_vtsq_host, bgk->fixed_vtsq);
}

void
pkpm_species_bgk_kinetic_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  struct timespec tm = gkyl_wall_clock();
  bgk->kinetic_rhs(app, s, bgk, fin, rhs);
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
}

// Backward Euler for self-collisions, S_q = 1/(1 + nu_q dt): n, u and the energy are unchanged,
// and the step is a convex combination of the state and the target, so it keeps the constraints of
// the state and takes F, G to the target and the first harmonic to zero as nu dt grows.
void
pkpm_species_bgk_implicit_advance(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk, double dt
)
{
  struct timespec tm = gkyl_wall_clock();
  int nc = app->confBasis.num_basis;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local, iter.idx);
    const double *rate_d = gkyl_array_cfetch(bgk->rate, loc);
    double *shape_d = gkyl_array_fetch(bgk->shape, loc);
    for (int q = 0; q < nc; ++q) {
      shape_d[q] = 1.0 / (1.0 + rate_d[q] * dt);
    }
  }
  gkyl_pkpm_bgk_collisions_implicit_advance(
    bgk->up_bgk, &app->local, &s->local_vel, &s->local, bgk->shape, app->pos_map->jacob_pos_gauss,
    bgk->f_lte, s->f
  );
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
}

// With a fixed temperature the fluid energy relaxes to rho |u|^2/2 + 3 rho xi_0^2/2 at the
// rate nu (density and momentum are conserved by the BGK operator).
void
pkpm_species_bgk_fluid_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk,
  struct gkyl_array *rhs
)
{
  if (!bgk->fixed_temp_relax) {
    return;
  }
  int nc = app->confBasis.num_basis;
  gkyl_array_clear(bgk->uu, 0.0);
  for (int i = 0; i < 3; ++i) {
    gkyl_dg_mul_op_range(&app->confBasis, 0, bgk->scratch, i, s->pkpm_u, i, s->pkpm_u, &app->local);
    gkyl_array_accumulate(bgk->uu, 1.0, bgk->scratch);
  }
  gkyl_dg_mul_op_range(&app->confBasis, 0, bgk->E_target, 0, s->pkpm_rho, 0, bgk->uu, &app->local);
  gkyl_dg_mul_op_range(
    &app->confBasis, 0, bgk->scratch, 0, s->pkpm_rho, 0, bgk->fixed_vtsq, &app->local
  );
  gkyl_array_scale(bgk->E_target, 0.5);
  gkyl_array_accumulate(bgk->E_target, 1.5, bgk->scratch);
  // physical E of the stage: block 4 of fluid_no_J
  gkyl_array_accumulate_offset(bgk->E_target, -1.0, s->fluid_no_J, 4 * nc);
  // the evolved fluid carries the position-map Jacobian
  gkyl_vlasov_position_map_rescale_jacobpos_conf(
    app->pos_map, &app->local, bgk->E_target, bgk->E_target
  );
  gkyl_dg_mul_op_range(&app->confBasis, 0, bgk->scratch, 0, bgk->nu, 0, bgk->E_target, &app->local);
  gkyl_array_accumulate_offset_range(rhs, 1.0, bgk->scratch, 4 * nc, &app->local);
}

void
pkpm_species_bgk_release(const gkyl_pkpm_app *app, const struct pkpm_bgk_collisions *bgk)
{
  if (bgk->use_nodal) {
    if (bgk->is_implicit) {
      gkyl_array_release(bgk->shape);
    }
    gkyl_array_release(bgk->rate);
    gkyl_array_release(bgk->f_lte);
    gkyl_pkpm_bgk_collisions_release(bgk->up_bgk);
    return;
  }
  gkyl_proj_on_basis_release(bgk->proj_target);
  gkyl_array_release(bgk->nu);
  gkyl_array_release(bgk->fixed_vtsq);
  gkyl_array_release(bgk->fixed_vtsq_host);
  gkyl_array_release(bgk->f_target);
  gkyl_array_release(bgk->f_target_host);
  gkyl_array_release(bgk->E_target);
  gkyl_array_release(bgk->scratch);
  gkyl_array_release(bgk->uu);
  if (app->use_gpu) {
    gkyl_array_release(bgk->xi_host);
  }
}
