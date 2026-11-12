#include <assert.h>
#include <gkyl_pkpm_priv.h>

struct gkyl_array *
pkpm_species_coll_profile_new(
  gkyl_pkpm_app *app, void (*eval)(double t, const double *xn, double *fout, void *ctx), void *ctx,
  double scale, bool is_nodal
)
{
  struct gkyl_array *modal = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &app->grid,
    .basis = &app->confBasis,
    .num_quad = app->poly_order + 1,
    .num_ret_vals = 1,
    .eval = eval,
    .ctx = ctx,
    .c2p_func = pkpm_pos_c2p,
    .c2p_func_ctx = app->pos_map,
  });
  gkyl_proj_on_basis_advance(proj, 0.0, &app->local_ext, modal);
  gkyl_proj_on_basis_release(proj);
  gkyl_array_scale(modal, scale);
  if (!is_nodal) {
    return modal;
  }

  struct gkyl_array *nodal = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  gkyl_dg_conf_nodal_conf_to_nodal(app->conf_nodal, &app->local_ext, modal, nodal);
  gkyl_array_release(modal);
  return nodal;
}

// Node values of zero stand in for a self-collision frequency that is computed every stage.
struct gkyl_array *
pkpm_species_coll_self_nu_new(gkyl_pkpm_app *app, const struct pkpm_species *s, double nu_frac)
{
  const struct gkyl_pkpm_collisions *coll = &s->info.collisions;
  if (coll->self_nu) {
    return pkpm_species_coll_profile_new(app, coll->self_nu, coll->ctx, nu_frac, true);
  }
  struct gkyl_array *nu = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  gkyl_array_clear(nu, 0.0);
  return nu;
}

struct gkyl_array *
pkpm_species_bgk_target_new(gkyl_pkpm_app *app, struct pkpm_species *s)
{
  struct gkyl_array *f_lte = mkarr(false, 3, s->local_ext_vel.volume);
  gkyl_pkpm_lte *lte = gkyl_pkpm_lte_inew(&(struct gkyl_pkpm_lte_inp){
    .vel_grid = &s->grid_vel,
    .vel_range = &s->local_vel,
    .num_quad = 4,
    .eps = 1.0e-14,
    .max_iter = 100,
  });
  double moms[3] = {1.0, 0.0, 1.0};
  struct gkyl_pkpm_lte_status status = gkyl_pkpm_lte_correct(lte, moms, f_lte);
  gkyl_pkpm_lte_release(lte);
  if (!status.converged) {
    gkyl_pkpm_app_cout(
      app, stderr, "WARNING: species %s: the BGK target moments reached a relative error of %g\n",
      s->info.name, status.error
    );
  }
  return f_lte;
}

// The unit Maxwellian in w and its isotropic partner G = F (the modal form of the target).
static void
eval_bgk_target(double t, const double *xn, double *restrict fout, void *ctx)
{
  const int *cdim = ctx;
  double w = xn[*cdim];
  double F = exp(-0.5 * w * w) / sqrt(2.0 * M_PI);
  fout[0] = F;
  fout[1] = F;
}

// nu (f_target - f) accumulated onto the right-hand side, with the collision frequency at the
// nodes: all blocks relax to the corrected unit Maxwellian.
static void
pkpm_species_bgk_kinetic_rhs_explicit(
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

// nu (f_target - f) in the modal form, for GPU runs. The target is gyrotropic ([F, G] only), so
// the first harmonic [M0_perp, N_perp] relaxes to zero. This form adds nothing to the CFL rate.
static void
pkpm_species_bgk_kinetic_rhs_modal(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk,
  const struct gkyl_array *fin, struct gkyl_array *rhs
)
{
  pkpm_species_accumulate_conf_phase(app, s, rhs, 1.0, bgk->nu, bgk->f_target);
  pkpm_species_accumulate_conf_phase(app, s, rhs, -1.0, bgk->nu, fin);
}

// Backward Euler for self-collisions, S_q = 1/(1 + nu_q dt): n, u and the energy are unchanged.
static void
pkpm_species_bgk_shape_self(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk, double dt
)
{
  int num_nodes = app->confBasis.num_basis;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local, iter.idx);
    const double *rate_d = gkyl_array_cfetch(bgk->rate_self, loc);
    double *shape_d = gkyl_array_fetch(bgk->shape, loc);
    for (int q = 0; q < num_nodes; ++q) {
      shape_d[q] = 1.0 / (1.0 + rate_d[q] * dt);
    }
  }
}

// Backward Euler with partners: S_q = min(1, (T_q / T_q') / (1 + (nu_q + sum_r nu_sr_q) dt)) with
// the temperatures before and after the exchange of the stage, so that every pressure component
// takes the backward-Euler step at the total frequency.
static void
pkpm_species_bgk_shape_cross(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk, double dt
)
{
  pkpm_species_coll_cross_shape(app, s, dt, bgk->shape);
}

// The collision frequency at the nodes and the single velocity-space target of the operator:
// the unit Maxwellian in w, with its discrete moments corrected to [1, 0, 1].
static void
pkpm_species_bgk_init_nodal(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk, double nu_frac
)
{
  bgk->rate_self = pkpm_species_coll_self_nu_new(app, s, nu_frac);
  // Without partners the self frequency is the rate of the operator. With partners the rate is
  // rebuilt from it every stage.
  if (s->info.collisions.num_cross_collisions) {
    bgk->rate = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    gkyl_array_copy(bgk->rate, bgk->rate_self);
  } else {
    bgk->rate = gkyl_array_acquire(bgk->rate_self);
  }
  bgk->f_lte = pkpm_species_bgk_target_new(app, s);
  bgk->up_bgk = gkyl_pkpm_bgk_collisions_new(&app->confBasis, &app->basis, false);
  bgk->is_implicit = s->info.collisions.is_implicit;
  if (bgk->is_implicit) {
    bgk->shape = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
    bgk->kinetic_rhs = pkpm_species_bgk_kinetic_rhs_none;
    bgk->calc_shape = s->info.collisions.num_cross_collisions ? pkpm_species_bgk_shape_cross :
                                                                pkpm_species_bgk_shape_self;
  } else {
    bgk->kinetic_rhs = pkpm_species_bgk_kinetic_rhs_explicit;
  }
}

// The modal collision frequency and the projected target J [F_t, G_t], which does not change
// in time.
static void
pkpm_species_bgk_init_modal(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk, double nu_frac
)
{
  // The computed frequency is not available in this form.
  assert(s->info.collisions.self_nu);
  struct gkyl_array *nu_host = pkpm_species_coll_profile_new(
    app, s->info.collisions.self_nu, s->info.collisions.ctx, nu_frac, false
  );
  bgk->nu = mkarr(app->use_gpu, nu_host->ncomp, nu_host->size);
  gkyl_array_copy(bgk->nu, nu_host);
  gkyl_array_release(nu_host);

  int num_blocks = 2, cdim = app->cdim;
  struct gkyl_array *f_target_host =
    mkarr(false, num_blocks * app->basis.num_basis, s->local_ext.volume);
  gkyl_proj_on_basis *proj = gkyl_proj_on_basis_inew(&(struct gkyl_proj_on_basis_inp){
    .grid = &s->grid,
    .basis = &app->basis,
    .num_quad = 4,
    .num_ret_vals = num_blocks,
    .eval = eval_bgk_target,
    .ctx = &cdim,
  });
  gkyl_array_clear(f_target_host, 0.0);
  gkyl_proj_on_basis_advance(proj, 0.0, &s->local, f_target_host);
  gkyl_proj_on_basis_release(proj);
  bgk->f_target = mkarr(app->use_gpu, f_target_host->ncomp, f_target_host->size);
  gkyl_array_copy(bgk->f_target, f_target_host);
  gkyl_array_release(f_target_host);
  gkyl_vlasov_position_map_rescale_jacobpos(
    app->pos_map, &app->basis, &s->local, bgk->f_target, bgk->f_target
  );
  bgk->kinetic_rhs = pkpm_species_bgk_kinetic_rhs_modal;
}

void
pkpm_species_bgk_init(gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk)
{
  double nu_frac = s->info.collisions.nuFrac > 0.0 ? s->info.collisions.nuFrac : 1.0;

  // GPU runs take the explicit modal form, which has no implicit step.
  assert(!(app->use_gpu && s->info.collisions.is_implicit));
  bgk->use_nodal = !app->use_gpu;
  bgk->is_implicit = false;
  if (bgk->use_nodal) {
    pkpm_species_bgk_init_nodal(app, s, bgk, nu_frac);
  } else {
    pkpm_species_bgk_init_modal(app, s, bgk, nu_frac);
  }
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

// Backward Euler of the operator: the step is a convex combination of the state and the target,
// so it keeps the constraints of the state and takes F, G to the target and the first harmonic
// to zero as nu dt grows.
void
pkpm_species_bgk_implicit_advance(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_bgk_collisions *bgk, double dt,
  struct gkyl_array *f
)
{
  struct timespec tm = gkyl_wall_clock();
  bgk->calc_shape(app, s, bgk, dt);
  gkyl_pkpm_bgk_collisions_implicit_advance(
    bgk->up_bgk, &app->local, &s->local_vel, &s->local, bgk->shape, app->pos_map->jacob_pos_gauss,
    bgk->f_lte, f
  );
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
}

void
pkpm_species_bgk_release(const gkyl_pkpm_app *app, const struct pkpm_bgk_collisions *bgk)
{
  if (bgk->use_nodal) {
    if (bgk->is_implicit) {
      gkyl_array_release(bgk->shape);
    }
    gkyl_array_release(bgk->rate_self);
    gkyl_array_release(bgk->rate);
    gkyl_array_release(bgk->f_lte);
    gkyl_pkpm_bgk_collisions_release(bgk->up_bgk);
    return;
  }
  gkyl_array_release(bgk->nu);
  gkyl_array_release(bgk->f_target);
}
