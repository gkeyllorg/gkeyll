#include <gkyl_pkpm_priv.h>

// Take a forward Euler step of the PKPM system of equations with the suggested time-step dt.
// Note: this may not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by stability.
// The actual time-step and dt_suggested are returned in the status object.
void
pkpm_forward_euler(
  gkyl_pkpm_app *app, int stage, double tcurr, double dt, struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *fout[],
  struct gkyl_array *fluidout[], struct gkyl_array *emout, struct gkyl_update_status *st
)
{
  app->stat.nfeuler += 1;

  double dtmin = DBL_MAX;

  // Compute applied acceleration, external EM fields, and applied currents if present.
  // Note these are only computed here *if* momentum-EM coupling is explicit
  // Otherwise, applies acceleration, external EM fields, and applied currents are computed in implicit update.
  // Also, these quantities use proj_on_basis so we copy to GPU every call if app->use_gpu = true.
  if (app->field->ext_em_evolve && app->use_explicit_sources) {
    pkpm_field_calc_ext_em(app, app->field, tcurr);
  }
  if (app->field->app_current_evolve && app->use_explicit_sources) {
    pkpm_field_calc_app_current(app, app->field, tcurr);
  }
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].fluid->app_accel_evolve && app->use_explicit_sources) {
      pkpm_species_calc_app_accel(app, &app->species[i], tcurr);
    }
  }

  // Compute the magnetic field unit vector b = B/|B|
  pkpm_field_calc_bvar(app, app->field, emin);

  // RHS of Maxwell's equations first: it does not depend on the species, and its magnetic part
  // (stored physical in field->dt_em_no_J) is the d_t B of the first-harmonic force terms.
  double dtmin_field = pkpm_field_rhs(app, app->field, stage, emin, emout);

  // Everything of the stage that reads a species and the field alone: the moments, the
  // primitive variables, the collision moments, the face fluxes and the wall losses.
  for (int i = 0; i < app->num_species; ++i) {
    pkpm_species_calc_self_moms(app, &app->species[i], fin[i], fluidin[i]);
  }

  // What reads other species (the rate of the kinetic operator with its partners, the rate of
  // the boundary-flux sources) after all species' self moments, hence a separate loop.
  for (int i = 0; i < app->num_species; ++i) {
    pkpm_species_calc_coupled_vars(app, &app->species[i], fin[i], fluidin[i]);
  }

  // compute RHS of pkpm equations
  for (int i = 0; i < app->num_species; ++i) {
    // the species see the physical fields of the stage (computed by pkpm_field_calc_bvar)
    double dt1 = pkpm_species_rhs(
      app, &app->species[i], fin[i], fluidin[i], app->field->em_no_J, fout[i], fluidout[i]
    );
    dtmin = fmin(dtmin, dt1);
  }

  dtmin = fmin(dtmin, dtmin_field);

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt - dtmin) / dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff) {
    dtmin = dt;
  }

  // compute minimum time-step across all processors
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  dtmin = dtmin_global;

  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete update of the species (distributions and fluid state)
  for (int i = 0; i < app->num_species; ++i) {
    pkpm_species_step_f(&app->species[i], dta, fin[i], fluidin[i], fout[i], fluidout[i]);
  }
  // then what follows the step: the integrating-factor stage of the first harmonic and the
  // positivity limiter, and the boundary conditions
  for (int i = 0; i < app->num_species; ++i) {
    pkpm_species_finish_stage(app, &app->species[i], stage, dta, fout[i], fluidout[i]);
    pkpm_species_apply_bc(app, app->species[i].kinetic, fout[i]);
    pkpm_fluid_species_apply_bc(app, app->species[i].fluid, fluidout[i]);
  }

  // Explicit accumulation current contribution from momentum to electric field terms
  if (app->use_explicit_sources) {
    pkpm_field_explicit_accumulate_current(app, app->field, fluidin, emout);
  }
  // complete update of field (even when field is static, it is
  // safest to do this accumulate as it ensure emout = emin)
  gkyl_array_accumulate(gkyl_array_scale(emout, dta), 1.0, emin);
  pkpm_field_apply_bc(app, app->field, emout);
}