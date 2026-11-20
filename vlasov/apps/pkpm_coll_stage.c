#include <gkyl_pkpm_priv.h>

// Implicit collisions and sources inside the SSP-RK3 stages: the cross-species momentum and
// energy exchange of every collision operator, the implicit BGK relaxation of the distributions
// and the implicit temperature source. After the forward-Euler step and the combination of
// stage i, whose state is in f1 and fluid1 of each species, they take one backward-Euler step
// over the time weight of the stage: dt, dt/4 and 2 dt/3. The step has three phases over all
// species: the moments of the stage state (the exchange, which needs every species, and the
// energy step and frequencies of each species), the implicit update of the distributions, then
// the limiter and the boundary conditions. A species with an implicit temperature source and
// implicit BGK collisions takes the two in one step.
//
// An implicit relaxation of the distributions that carries a temperature change (implicit BGK
// with partners, implicit temperature source) rescales the normalized distributions in the
// implicit step, after the combination. For those species the combination of stages 2 and 3 is
// weighted by the temperature factor lam the implicit steps of this time step have applied, so
// that it is linear in the physical second moments for any temperature change
// (gkyl_pkpm_stage_weights). Operators whose forward-Euler step carries the temperature change
// in its rate (LBO, explicit BGK) keep the SSP-RK3 weights.
static const double pkpm_coll_stage_weight[3] = {1.0, 1.0 / 4.0, 2.0 / 3.0};
// SSP-RK3 weight of the forward-Euler output in the combination of each stage.
static const double pkpm_coll_stage_coef_fe[3] = {1.0, 1.0 / 4.0, 2.0 / 3.0};

// The weighted combination of a stage (2 or 3) from the plain one in f1, with the temperature
// factor of the state behind the forward-Euler output fnew; at stage 1 there is no combination
// and the factor starts from 1.
static void
pkpm_coll_stage_combine(gkyl_pkpm_app *app, struct pkpm_species *s, int stage)
{
  if (stage == 1) {
    gkyl_array_clear_range(s->kinetic->stage_lam, 1.0, &app->local);
    return;
  }
  gkyl_pkpm_stage_weights_advance(
    s->kinetic->up_stage_weights, &app->local, &s->kinetic->local,
    pkpm_coll_stage_coef_fe[stage - 1], s->kinetic->f, s->kinetic->fnew, s->kinetic->stage_lam,
    s->kinetic->f1
  );
}

static bool
coll_stage_has_implicit_bgk(const struct pkpm_species *s)
{
  return s->kinetic->collision_id == GKYL_BGK_COLLISIONS && s->kinetic->bgk.is_implicit;
}

void
pkpm_coll_stage_advance(gkyl_pkpm_app *app, int stage, double dt)
{
  double dt_stage = pkpm_coll_stage_weight[stage - 1] * dt;
  struct gkyl_array *fluid[GKYL_MAX_SPECIES];
  for (int i = 0; i < app->num_species; ++i) {
    fluid[i] = app->species[i].fluid->fluid1;
  }
  for (int i = 0; i < app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    if (s->kinetic->has_stage_weights) {
      pkpm_coll_stage_combine(app, s, stage);
    }
  }
  // Moments of the stage state. The exchange of the species that collide with each other comes
  // first: it also forms the node values and frequencies the implicit steps read. Then the
  // energy step of the implicit temperature source and the shape factor of the implicit step of
  // each species.
  pkpm_coll_cross_advance(app, dt_stage, fluid);
  for (int i = 0; i < app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    if (s->kinetic->src_bgk.is_implicit) {
      pkpm_species_source_bgk_implicit_moms(app, s, &s->kinetic->src_bgk, dt_stage, fluid[i]);
    } else if (coll_stage_has_implicit_bgk(s)) {
      pkpm_species_bgk_implicit_moms(app, s, &s->kinetic->bgk, dt_stage);
    }
  }

  // The implicit update of the distributions: the source with the implicit collisions of the
  // species, or the collisions alone.
  for (int i = 0; i < app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    if (s->kinetic->src_bgk.is_implicit) {
      pkpm_species_source_bgk_implicit_advance(app, s, &s->kinetic->src_bgk, s->kinetic->f1);
    } else if (coll_stage_has_implicit_bgk(s)) {
      pkpm_species_bgk_implicit_advance(app, s, &s->kinetic->bgk, s->kinetic->f1);
    }
  }

  // The temperature factors of the stage for the weighted combination, the vertex limiter where
  // the temperature source changed the energy at the nodes, and the boundary conditions.
  for (int i = 0; i < app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    bool has_implicit_bgk = coll_stage_has_implicit_bgk(s);
    if (has_implicit_bgk && s->kinetic->has_stage_weights &&
        s->kinetic->cross.num_cross_collisions > 0) {
      gkyl_pkpm_stage_weights_scale(
        s->kinetic->up_stage_weights, &app->local, s->kinetic->cross.temp_ratio,
        s->kinetic->stage_lam
      );
    }
    if (s->kinetic->src_bgk.is_implicit) {
      gkyl_pkpm_stage_weights_scale(
        s->kinetic->up_stage_weights, &app->local, s->kinetic->src_bgk.temp_ratio,
        s->kinetic->stage_lam
      );
      pkpm_species_apply_positivity(app, s->fluid, fluid[i]);
      pkpm_fluid_species_apply_bc(app, s->fluid, fluid[i]);
    }
    if (has_implicit_bgk || s->kinetic->src_bgk.is_implicit) {
      pkpm_species_apply_bc(app, s->kinetic, s->kinetic->f1);
    }
  }
}
