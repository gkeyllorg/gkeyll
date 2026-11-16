#include <gkyl_vlasov_priv.h>

// Take a forward Euler step of the Vlasov-Maxwell system of equations
// with the suggested time-step dt. Also supports just Maxwell's equations
// and fluid equations (Euler's) with potential Vlasov-fluid coupling.
// Note: this may not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by stability.
// The actual time-step and dt_suggested are returned in the status object.
// Species loops run over the overall species count; the RK-state arrays are
// NULL where a species lacks that aspect. 'stage' is the SSP-RK3 stage (1, 2
// or 3) this step belongs to.
void
vlasov_forward_euler(
  gkyl_vlasov_app *app, int stage, double tcurr, double dt, struct gkyl_array *fin[],
  const struct gkyl_array *fluidin[], const struct gkyl_array *emin, struct gkyl_array *fout[],
  struct gkyl_array *fluidout[], struct gkyl_array *emout, struct gkyl_update_status *st
)
{
  app->stat.nfeuler += 1;

  double dtmin = DBL_MAX;
  int num_species = app->num_species;

  // The field reads the stage distributions without changing them.
  const struct gkyl_array *fin_field[num_species];
  for (int i = 0; i < num_species; ++i) {
    fin_field[i] = fin[i];
  }

  // Compute external EM field or applied currents if present and time-dependent.
  // Note: external EM field and  applied currents use proj_on_basis
  // so does copy to GPU every call if app->use_gpu = true.
  if (app->field->app_current_evolve && !app->has_fluid_em_coupling) {
    vlasov_field_calc_app_current(app, tcurr);
  }
  if (app->field->ext_em_evolve) {
    vlasov_field_calc_ext_em(app, tcurr);
  }
  if (app->field->ext_pot_evolve) {
    vlasov_field_calc_ext_pot(app, tcurr);
  }
  // Compute applied acceleration if present and time-dependent.
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_calc_app_accel(app, &app->species[i], tcurr);
  }

  // Update the field at the start of the step: the Maxwell RHS, or the Poisson
  // solve for the potential at the current time (read by the species RHS below).
  double dt1_field = vlasov_field_update(app, stage, tcurr, fin_field, emin, emout);
  dtmin = fmin(dtmin, dt1_field); // null field returns DBL_MAX (no constraint)

  // Compute what each species needs of itself and the field before its RHS
  // (self-collision moments and boundary corrections).
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_calc_self_moms(app, &app->species[i], fin[i], fluidin[i]);
  }

  // Compute the coupled variables (LBO cross-collision moments, fluid primitive
  // variables) after all species' self moments, hence a separate loop.
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_calc_coupled_vars(app, &app->species[i], fin[i], fluidin[i]);
  }

  // Compute RHS of the Vlasov/fluid equations.
  for (int i = 0; i < num_species; ++i) {
    double dt1 =
      vlasov_species_rhs(app, &app->species[i], fin[i], fluidin[i], emin, fout[i], fluidout[i]);
    dtmin = fmin(dtmin, dt1);
  }

  // Compute source term. Done here as the RHS update for all species should be
  // complete in case we need a bflux calculation for the source species.
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_calc_source_moms(app, &app->species[i], fin[i]);
  }
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_source_rhs(app, &app->species[i], tcurr, fout[i], fluidout[i]);
  }

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

  // Complete the update of the species (distribution and/or fluid).
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_step_f(&app->species[i], dta, fin[i], fluidin[i], fout[i], fluidout[i]);
  }
  // Then what a species does to its stepped state before the boundary
  // conditions are applied to it.
  for (int i = 0; i < num_species; ++i) {
    vlasov_species_finish_stage(app, &app->species[i], stage, dta, fout[i], fluidout[i]);
  }

  // Complete the field update: for Vlasov-Maxwell, accumulate the species
  // current onto the RHS and finalize emout = emin + dta*RHS.
  vlasov_field_complete_update(app, dta, fin_field, fluidin, emin, emout);
}
