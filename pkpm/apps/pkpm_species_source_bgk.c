#include <assert.h>
#include <gkyl_pkpm_priv.h>

// The temperature source S = nu_Q (f_M(n, u, T_Q) - f) in the normalized variables, at the
// configuration-space nodes: the thermal energy of the fluid relaxes to (3/2) n T_Q at the rate
// nu_Q, and the distributions relax to the unit Maxwellian at the rate nu_Q T_Q/T, the BGK
// operator at that rate. Density, momentum and the constraints are untouched, so the closures
// receive nothing from it.

// Explicit update, fluid half: the energy relaxation of the stage into the array the fluid
// right-hand side adds, with its CFL rate, and the rate of the distributions.
static void
pkpm_species_source_bgk_calc_fluid_explicit(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src,
  const struct gkyl_array *fluidin
)
{
  gkyl_array_clear(s->fluid->src_bgk_rhs, 0.0);
  gkyl_pkpm_temp_source_advance(
    src->up_temp, &app->local, src->rate, src->vtsq, fluidin, s->fluid->src_bgk_rhs,
    s->fluid->cflrate, src->shape
  );
}

// Explicit update, kinetic half: the relaxation of the distributions onto the right-hand side of
// the kinetic equations, with its CFL rate.
static void
pkpm_species_source_bgk_kinetic_rhs_explicit(
  gkyl_pkpm_app *app, struct pkpm_kinetic_species *kin, struct pkpm_source_bgk *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs_f
)
{
  gkyl_pkpm_bgk_collisions_advance(
    src->up_bgk, &app->local, &kin->local_vel, &kin->local, src->shape,
    app->pos_map->jacob_pos_gauss, src->f_lte, fin, rhs_f, kin->cflrate
  );
}

// No source, or the implicit one: nothing enters the right-hand sides of the explicit update.
static void
pkpm_species_source_bgk_calc_fluid_none(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src,
  const struct gkyl_array *fluidin
)
{
}

static void
pkpm_species_source_bgk_kinetic_rhs_none(
  gkyl_pkpm_app *app, struct pkpm_kinetic_species *kin, struct pkpm_source_bgk *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs_f
)
{
}

void
pkpm_species_source_bgk_init(gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src)
{
  const struct gkyl_vlasov_source_bgk *info = &s->kinetic->info.source_bgk;
  src->source_bgk_id = info->source_bgk_id;
  src->is_implicit = false;
  src->calc_fluid = pkpm_species_source_bgk_calc_fluid_none;
  src->kinetic_rhs = pkpm_species_source_bgk_kinetic_rhs_none;
  if (src->source_bgk_id == GKYL_NO_SOURCE_BGK) {
    return;
  }
  // The prescribed-temperature source is the only model, and its updaters are not on the GPU.
  assert(src->source_bgk_id == GKYL_SOURCE_BGK_HEATING);
  assert(info->rate_profile && info->temp_shape);
  assert(!app->use_gpu);

  src->rate =
    pkpm_species_coll_profile_new(app, info->rate_profile, info->rate_profile_ctx, 1.0, true);
  src->vtsq =
    pkpm_species_coll_profile_new(app, info->temp_shape, info->temp_shape_ctx, 1.0 / s->mass, true);
  src->shape = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  src->temp_ratio = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  src->f_lte = pkpm_species_bgk_target_new(app, s);
  src->up_temp = gkyl_pkpm_temp_source_inew(&(struct gkyl_pkpm_temp_source_inp){
    .conf_basis = &app->confBasis,
    .use_gpu = false,
  });
  src->up_bgk = gkyl_pkpm_bgk_collisions_inew(&(struct gkyl_pkpm_bgk_collisions_inp){
    .conf_basis = &app->confBasis,
    .phase_basis = &app->basis,
    .use_gpu = false,
  });

  src->is_implicit = info->is_implicit;
  if (!src->is_implicit) {
    src->calc_fluid = pkpm_species_source_bgk_calc_fluid_explicit;
    src->kinetic_rhs = pkpm_species_source_bgk_kinetic_rhs_explicit;
    // the energy relaxation of the stage, for the fluid right-hand side
    s->fluid->src_bgk_rhs = mkarr(false, 5 * app->confBasis.num_basis, app->local_ext.volume);
    s->fluid->has_src_bgk_rhs = true;
  }
}

void
pkpm_species_source_bgk_calc_fluid(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src,
  const struct gkyl_array *fluidin
)
{
  struct timespec tm = gkyl_wall_clock();
  src->calc_fluid(app, s, src, fluidin);
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
}

void
pkpm_species_source_bgk_kinetic_rhs(
  gkyl_pkpm_app *app, struct pkpm_kinetic_species *kin, struct pkpm_source_bgk *src,
  const struct gkyl_array *fin, struct gkyl_array *rhs_f
)
{
  struct timespec tm = gkyl_wall_clock();
  src->kinetic_rhs(app, kin, src, fin, rhs_f);
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
}

// Backward Euler: the thermal energy first, eps <- (eps + nu_Q dt eps_Q)/(1 + nu_Q dt), then
// the distributions with S = 1/(1 + nu_Q dt eps_Q/eps), eps being the energy before the step.
// S times the new energy over the old one is 1/(1 + nu_Q dt), so every pressure component, the
// product of the energy and a normalized moment, takes the backward-Euler step of the source.
// Both updates are convex combinations at the nodes.
void
pkpm_species_source_bgk_implicit_advance(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src, double dt,
  struct gkyl_array *f, struct gkyl_array *fluid
)
{
  struct timespec tm = gkyl_wall_clock();
  gkyl_pkpm_temp_source_implicit_advance(
    src->up_temp, &app->local, dt, src->rate, src->vtsq, fluid, src->shape, src->temp_ratio
  );
  gkyl_pkpm_bgk_collisions_implicit_advance(
    src->up_bgk, &app->local, &s->kinetic->local_vel, &s->kinetic->local, src->shape,
    app->pos_map->jacob_pos_gauss, src->f_lte, f
  );
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
}

void
pkpm_species_source_bgk_release(const gkyl_pkpm_app *app, const struct pkpm_source_bgk *src)
{
  if (src->source_bgk_id == GKYL_NO_SOURCE_BGK) {
    return;
  }
  gkyl_array_release(src->rate);
  gkyl_array_release(src->vtsq);
  gkyl_array_release(src->shape);
  gkyl_array_release(src->temp_ratio);
  gkyl_array_release(src->f_lte);
  gkyl_pkpm_temp_source_release(src->up_temp);
  gkyl_pkpm_bgk_collisions_release(src->up_bgk);
}
