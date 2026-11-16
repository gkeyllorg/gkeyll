#include <gkyl_vlasov_priv.h>

// The implicit operators of an SSP-RK3 stage. They act on the stage state (the
// f1 buffer) once the explicit part of the stage, its boundary conditions and
// the stage combination are done.

// Implicit BGK collisions in three phases over all species: moments of the
// stage state, the implicit update (which forms cross moments, so it needs all
// species' moments first), then BCs and the copy back into the stage state.
static void
stage_implicit_coll(gkyl_vlasov_app *app, double dt)
{
  int num_species = app->num_species;
  struct gkyl_array *fstage[num_species];
  struct gkyl_array *fout[num_species];
  const struct gkyl_array *fluidin[num_species];
  struct gkyl_array *fluidout[num_species];
  vlasov_species_gather_rk_state(app, VM_RK_F1, VM_RK_FNEW, fstage, fout, fluidin, fluidout);

  for (int i = 0; i < num_species; ++i) {
    vlasov_species_calc_implicit_moms(app, &app->species[i], fstage[i]);
  }

  for (int i = 0; i < num_species; ++i) {
    vlasov_species_rhs_implicit(app, &app->species[i], fstage[i], fout[i], dt);
  }

  for (int i = 0; i < num_species; ++i) {
    vlasov_species_finish_implicit_update(app, &app->species[i], fout[i], fstage[i], app->tcurr);
  }
}

// Share of the time-step each stage adds to its state: stage 1 is a full Euler
// step, and the combinations of stages 2 and 3 carry 1/4 and 2/3 of theirs.
static const double stage_dt_frac[3] = {1.0, 1.0 / 4.0, 2.0 / 3.0};

void
vlasov_stage_implicit(gkyl_vlasov_app *app, int stage, double dt)
{
  // The implicit BGK operator takes one backward-Euler step over the stage's
  // share of the time-step, on the state the stage produced.
  if (app->has_implicit_coll_scheme) {
    stage_implicit_coll(app, stage_dt_frac[stage - 1] * dt);
  }
}
