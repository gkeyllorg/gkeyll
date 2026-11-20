#include <gkyl_vlasov_priv.h>

// The implicit operators of an SSP-RK3 stage. They act on the stage state (the
// f1 buffer) once the explicit part of the stage, its boundary conditions and
// the stage combination are done.

// Implicit BGK collisions and implicit temperature sources, one backward-Euler
// step of both per species, in three phases over all species: moments of the
// stage state, the implicit update (which forms cross moments, so it needs all
// species' moments first), then BCs and the copy back into the stage state.
static void
stage_implicit_species(gkyl_vlasov_app *app, double dt)
{
  int num_species = app->num_species;
  struct gkyl_array *fstage[num_species];
  struct gkyl_array *fout[num_species];
  const struct gkyl_array *fluidin[num_species];
  struct gkyl_array *fluidout[num_species];
  vlasov_species_gather_rk_state(app, VM_RK_F1, VM_RK_FNEW, fstage, fout, fluidin, fluidout);

  // Kinetic species only: the implicit operators of PKPM species have their own step.
  for (int i = 0; i < num_species; ++i) {
    if (app->species[i].type == GKYL_SPECIES_VLASOV) {
      vlasov_species_calc_implicit_moms(app, &app->species[i], fstage[i], dt);
    }
  }

  for (int i = 0; i < num_species; ++i) {
    if (app->species[i].type == GKYL_SPECIES_VLASOV) {
      vlasov_species_rhs_implicit(app, &app->species[i], fstage[i], fout[i], dt);
    }
  }

  for (int i = 0; i < num_species; ++i) {
    if (app->species[i].type == GKYL_SPECIES_VLASOV) {
      vlasov_species_finish_implicit_update(app, &app->species[i], fout[i], fstage[i], app->tcurr);
    }
  }
}

// ---- PKPM species ----

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
static const double vm_pkpm_coll_stage_weight[3] = {1.0, 1.0 / 4.0, 2.0 / 3.0};
// SSP-RK3 weight of the forward-Euler output in the combination of each stage.
static const double vm_pkpm_coll_stage_coef_fe[3] = {1.0, 1.0 / 4.0, 2.0 / 3.0};

// The weighted combination of a stage (2 or 3) from the plain one in f1, with the temperature
// factor of the state behind the forward-Euler output fnew; at stage 1 there is no combination
// and the factor starts from 1.
static void
vm_pkpm_coll_stage_combine(gkyl_vlasov_app *app, struct vlasov_species *s, int stage)
{
  if (stage == 1) {
    gkyl_array_clear_range(s->kinetic->pkpm.stage_lam, 1.0, &app->local);
    return;
  }
  gkyl_pkpm_stage_weights_advance(
    s->kinetic->pkpm.up_stage_weights, &app->local, &s->kinetic->local,
    vm_pkpm_coll_stage_coef_fe[stage - 1], s->kinetic->f, s->kinetic->fnew,
    s->kinetic->pkpm.stage_lam, s->kinetic->f1
  );
}

static bool
coll_stage_has_implicit_bgk(const struct vlasov_species *s)
{
  return s->kinetic->pkpm.collision_id == GKYL_BGK_COLLISIONS && s->kinetic->pkpm.bgk.is_implicit;
}

void
vm_pkpm_coll_stage_advance(gkyl_vlasov_app *app, int stage, double dt)
{
  double dt_stage = vm_pkpm_coll_stage_weight[stage - 1] * dt;
  struct gkyl_array *fluid[GKYL_MAX_SPECIES];
  for (int i = 0; i < app->num_species; ++i) {
    fluid[i] = app->species[i].fluid->fluid1;
  }
  for (int i = 0; i < app->num_species; ++i) {
    struct vlasov_species *s = &app->species[i];
    if (s->kinetic->pkpm.has_stage_weights) {
      vm_pkpm_coll_stage_combine(app, s, stage);
    }
  }
  // Moments of the stage state. The exchange of the species that collide with each other comes
  // first: it also forms the node values and frequencies the implicit steps read. Then the
  // energy step of the implicit temperature source and the shape factor of the implicit step of
  // each species.
  vm_pkpm_coll_cross_advance(app, dt_stage, fluid);
  for (int i = 0; i < app->num_species; ++i) {
    struct vlasov_species *s = &app->species[i];
    if (s->kinetic->pkpm.src_bgk.is_implicit) {
      vm_pkpm_species_source_bgk_implicit_moms(
        app, s, &s->kinetic->pkpm.src_bgk, dt_stage, fluid[i]
      );
    } else if (coll_stage_has_implicit_bgk(s)) {
      vm_pkpm_species_bgk_implicit_moms(app, s, &s->kinetic->pkpm.bgk, dt_stage);
    }
  }

  // The implicit update of the distributions: the source with the implicit collisions of the
  // species, or the collisions alone.
  for (int i = 0; i < app->num_species; ++i) {
    struct vlasov_species *s = &app->species[i];
    if (s->kinetic->pkpm.src_bgk.is_implicit) {
      vm_pkpm_species_source_bgk_implicit_advance(app, s, &s->kinetic->pkpm.src_bgk, s->kinetic->f1);
    } else if (coll_stage_has_implicit_bgk(s)) {
      vm_pkpm_species_bgk_implicit_advance(app, s, &s->kinetic->pkpm.bgk, s->kinetic->f1);
    }
  }

  // The temperature factors of the stage for the weighted combination, the vertex limiter where
  // the temperature source changed the energy at the nodes, and the boundary conditions.
  for (int i = 0; i < app->num_species; ++i) {
    struct vlasov_species *s = &app->species[i];
    bool has_implicit_bgk = coll_stage_has_implicit_bgk(s);
    if (has_implicit_bgk && s->kinetic->pkpm.has_stage_weights &&
        s->kinetic->pkpm.cross.num_cross_collisions > 0) {
      gkyl_pkpm_stage_weights_scale(
        s->kinetic->pkpm.up_stage_weights, &app->local, s->kinetic->pkpm.cross.temp_ratio,
        s->kinetic->pkpm.stage_lam
      );
    }
    if (s->kinetic->pkpm.src_bgk.is_implicit) {
      gkyl_pkpm_stage_weights_scale(
        s->kinetic->pkpm.up_stage_weights, &app->local, s->kinetic->pkpm.src_bgk.temp_ratio,
        s->kinetic->pkpm.stage_lam
      );
      vm_pkpm_species_apply_positivity(app, s->fluid, fluid[i]);
      vm_pkpm_fluid_species_apply_bc(app, s->fluid, fluid[i]);
    }
    if (has_implicit_bgk || s->kinetic->pkpm.src_bgk.is_implicit) {
      vm_pkpm_species_apply_bc(app, s->kinetic, s->kinetic->f1);
    }
  }
}

// ---- the implicit operators of a stage ----

// Share of the time-step each stage adds to its state: stage 1 is a full Euler
// step, and the combinations of stages 2 and 3 carry 1/4 and 2/3 of theirs.
static const double stage_dt_frac[3] = {1.0, 1.0 / 4.0, 2.0 / 3.0};

void
vlasov_stage_implicit(gkyl_vlasov_app *app, int stage, double dt)
{
  // The implicit operators of the species take one backward-Euler step over the
  // stage's share of the time-step, on the state the stage produced.
  if (app->has_implicit_species_ops) {
    stage_implicit_species(app, stage_dt_frac[stage - 1] * dt);
  }
  // PKPM species: the exchange between them, implicit BGK collisions and implicit temperature
  // sources, with the weighted stage combination where the implicit step changes a temperature.
  if (app->has_pkpm) {
    vm_pkpm_coll_stage_advance(app, stage, dt);
  }
}
