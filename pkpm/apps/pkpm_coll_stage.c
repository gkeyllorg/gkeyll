#include <gkyl_pkpm_priv.h>

// Implicit collisions and sources inside the SSP-RK3 stages (PKPM_COLLISIONS_PLAN.md 8.4, scheme
// A): the cross-species momentum and energy exchange of every collision operator, the implicit
// BGK relaxation of the distributions and the implicit temperature source. After the
// forward-Euler step and the combination of stage i, whose state is in f1 and fluid1 of each
// species, each takes one backward-Euler step over the time weight of the stage: dt, dt/4 and
// 2 dt/3. Every stage is then a convex combination of forward-Euler steps and the targets for any
// dt, and it has the balance of transport and collisions as its fixed point, so the relaxed
// quantities the transport sees do not depend on dt (no conduction or friction of order dt).
static const double pkpm_coll_stage_weight[3] = {1.0, 1.0 / 4.0, 2.0 / 3.0};

void
pkpm_coll_stage_advance(gkyl_pkpm_app *app, int stage, double dt)
{
  double dt_stage = pkpm_coll_stage_weight[stage - 1] * dt;
  struct gkyl_array *fluid[GKYL_MAX_SPECIES];
  for (int i = 0; i < app->num_species; ++i) {
    fluid[i] = app->species[i].fluid1;
  }
  // Exchange of the species that collide with each other, which also forms the node values and
  // frequencies the implicit BGK steps read.
  pkpm_coll_cross_advance(app, dt_stage, fluid);

  // The implicit operators act at the nodes on the same target, so they commute. The
  // temperature source changes the energy at the nodes, so the vertex limiter follows it.
  for (int i = 0; i < app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    bool is_implicit_bgk = s->collision_id == GKYL_BGK_COLLISIONS && s->bgk.is_implicit;
    if (is_implicit_bgk) {
      pkpm_species_bgk_implicit_advance(app, s, &s->bgk, dt_stage, s->f1);
    }
    if (s->src_bgk.is_implicit) {
      pkpm_species_source_bgk_implicit_advance(app, s, &s->src_bgk, dt_stage, s->f1, s->fluid1);
      pkpm_species_apply_positivity(app, s, s->fluid1);
      pkpm_fluid_species_apply_bc(app, s, s->fluid1);
    }
    if (is_implicit_bgk || s->src_bgk.is_implicit) {
      pkpm_species_apply_bc(app, s, s->f1);
    }
  }
}
