#include <assert.h>
#include <gkyl_pkpm_priv.h>

// The temperature source S = nu_Q (f_M(n, u, T_Q) - f) in the normalized variables, at the
// configuration-space nodes: the thermal energy of the fluid relaxes to (3/2) n T_Q at the rate
// nu_Q, and the distributions relax to the unit Maxwellian at the rate nu_Q T_Q/T, the BGK
// operator at that rate. Density, momentum and the constraints are untouched, so the closures
// receive nothing from it.

// Explicit update: the energy relaxation onto the fluid right-hand side and the relaxation of
// the distributions onto the kinetic one, each with its CFL rate.
static void
pkpm_species_source_bgk_rhs_explicit(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, struct gkyl_array *rhs_f,
  struct gkyl_array *rhs_fluid
)
{
  gkyl_pkpm_temp_source_advance(
    src->up_temp, &app->local, src->rate, src->vtsq, fluidin, rhs_fluid, s->cflrate_fluid,
    src->shape
  );
  gkyl_pkpm_bgk_collisions_advance(
    src->up_bgk, &app->local, &s->local_vel, &s->local, src->shape, app->pos_map->jacob_pos_gauss,
    src->f_lte, fin, rhs_f, s->cflrate_f
  );
}

// No source, or the implicit one: nothing enters the right-hand sides of the explicit update.
static void
pkpm_species_source_bgk_rhs_none(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, struct gkyl_array *rhs_f,
  struct gkyl_array *rhs_fluid
)
{
}

void
pkpm_species_source_bgk_init(gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src)
{
  const struct gkyl_pkpm_source_bgk *info = &s->info.source_bgk;
  src->source_bgk_id = info->source_bgk_id;
  src->is_implicit = false;
  src->rhs = pkpm_species_source_bgk_rhs_none;
  if (src->source_bgk_id == GKYL_NO_SOURCE_BGK) {
    return;
  }
  // The prescribed-temperature source is the only model, and its updaters are not on the GPU.
  assert(src->source_bgk_id == GKYL_SOURCE_BGK_HEATING);
  assert(info->rate_profile && info->temp_shape);
  assert(!app->use_gpu);

  src->rate =
    pkpm_species_coll_profile_new(app, info->rate_profile, info->rate_profile_ctx, 1.0, true);
  src->vtsq = pkpm_species_coll_profile_new(
    app, info->temp_shape, info->temp_shape_ctx, 1.0 / s->info.mass, true
  );
  src->shape = mkarr(false, app->confBasis.num_basis, app->local_ext.volume);
  src->f_lte = pkpm_species_bgk_target_new(app, s);
  src->up_temp = gkyl_pkpm_temp_source_new(&app->confBasis, false);
  src->up_bgk = gkyl_pkpm_bgk_collisions_new(&app->confBasis, &app->basis, false);

  src->is_implicit = info->is_implicit;
  if (!src->is_implicit) {
    src->rhs = pkpm_species_source_bgk_rhs_explicit;
  }
}

void
pkpm_species_source_bgk_rhs(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src,
  const struct gkyl_array *fin, const struct gkyl_array *fluidin, struct gkyl_array *rhs_f,
  struct gkyl_array *rhs_fluid
)
{
  struct timespec tm = gkyl_wall_clock();
  src->rhs(app, s, src, fin, fluidin, rhs_f, rhs_fluid);
  app->stat.species_coll_tm += gkyl_time_diff_now_sec(tm);
}

// Backward Euler: the thermal energy first, eps <- (eps + nu_Q dt eps_Q)/(1 + nu_Q dt), then
// the distributions with S = 1/(1 + nu_Q dt eps_Q/eps), eps being the energy before the step.
// S times the new energy over the old one is 1/(1 + nu_Q dt), so every pressure component, the
// product of the energy and a normalized moment, takes the backward-Euler step of the source.
// Both updates are convex combinations at the nodes.
void
pkpm_species_source_bgk_implicit_advance(
  gkyl_pkpm_app *app, struct pkpm_species *s, struct pkpm_source_bgk *src, double dt
)
{
  struct timespec tm = gkyl_wall_clock();
  gkyl_pkpm_temp_source_implicit_advance(
    src->up_temp, &app->local, dt, src->rate, src->vtsq, s->fluid, src->shape
  );
  gkyl_pkpm_bgk_collisions_implicit_advance(
    src->up_bgk, &app->local, &s->local_vel, &s->local, src->shape, app->pos_map->jacob_pos_gauss,
    src->f_lte, s->f
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
  gkyl_array_release(src->f_lte);
  gkyl_pkpm_temp_source_release(src->up_temp);
  gkyl_pkpm_bgk_collisions_release(src->up_bgk);
}
