#include <gkyl_pkpm_priv.h>

// Take a forward Euler step of the PKPM system of equations with the suggested time-step dt. 
// Note: this may not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by stability. 
// The actual time-step and dt_suggested are returned in the status object.
void
pkpm_forward_euler(gkyl_pkpm_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[], const struct gkyl_array *emin,
  struct gkyl_array *fout[], struct gkyl_array *fluidout[], struct gkyl_array *emout, 
  struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;

  double dtmin = DBL_MAX;

  // Compute applied acceleration, external EM fields, and applied currents if present.
  // Note these are only computed here *if* momentum-EM coupling is explicit
  // Otherwise, applies acceleration, external EM fields, and applied currents are computed in implicit update.
  // Also, these quantities use proj_on_basis so we copy to GPU every call if app->use_gpu = true.
  if (app->field->ext_em_evolve && app->use_explicit_source) {
    pkpm_field_calc_ext_em(app, app->field, tcurr);
  }
  if (app->field->app_current_evolve && app->use_explicit_source) {
    pkpm_field_calc_app_current(app, app->field, tcurr); 
  }
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].app_accel_evolve && app->use_explicit_source) {
      pkpm_species_calc_app_accel(app, &app->species[i], tcurr);
    }
  }

  // Compute magnetic field unit vector and tensor, and div(b)
  pkpm_field_calc_bvar(app, app->field, emin);  

  // Two separate loops over number of species to compute needed auxiliary quantities.
  for (int i=0; i<app->num_species; ++i) {
    // Compute parallel-kinetic-perpendicular moment (pkpm) variables.
    // These are the coupling moments [rho, p_par, p_perp], the self-consistent
    // pressure force (div(p_par b_hat)), and the primitive variables
    pkpm_species_calc_pkpm_vars(app, &app->species[i], fin[i], fluidin[i]);
    // compute necessary moments and boundary corrections for collisions
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      pkpm_species_lbo_moms(app, &app->species[i], &app->species[i].lbo, fin[i]);
    }
  }
  for (int i=0; i<app->num_species; ++i) {
    // compute necessary moments for cross-species collisions
    // needs to be done after self-collisions moments, so separate loop over species
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS
      && app->species[i].lbo.num_cross_collisions) {
      pkpm_species_lbo_cross_moms(app, &app->species[i], &app->species[i].lbo, fin[i]);
    }
    // Finish computing parallel-kinetic-perpendicular moment (pkpm) variables.
    // These are the update variables including the acceleration variables in the kinetic
    // equation and the source distribution functions for Laguerre couplings.
    // Needs to be done after all collisional moment computations for collisional sources
    // in Laguerre couplings.
    pkpm_species_calc_pkpm_update_vars(app, &app->species[i], fin[i]); 
  }

  // compute RHS of pkpm equations
  for (int i=0; i<app->num_species; ++i) {
    double dt1 = pkpm_species_rhs(app, &app->species[i], fin[i], fluidin[i], emin, fout[i], fluidout[i]);
    dtmin = fmin(dtmin, dt1);
  }

  // compute RHS of Maxwell equations
  double dt1 = pkpm_field_rhs(app, app->field, emin, emout);
  dtmin = fmin(dtmin, dt1);

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // compute minimum time-step across all processors
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  dtmin = dtmin_global;
  
  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete update of distribution function
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dta), 1.0, fin[i]);
    gkyl_array_accumulate(gkyl_array_scale(fluidout[i], dta), 1.0, fluidin[i]);
    pkpm_species_apply_bc(app, &app->species[i], fout[i]);
    pkpm_fluid_species_apply_bc(app, &app->species[i], fluidout[i]);
  }

  // Explicit accumulation current contribution from momentum to electric field terms
  if (app->use_explicit_source) {
    pkpm_field_explicit_accumulate_current(app, app->field, fluidin, emout);
  }
  // complete update of field (even when field is static, it is
  // safest to do this accumulate as it ensure emout = emin)
  gkyl_array_accumulate(gkyl_array_scale(emout, dta), 1.0, emin);
  pkpm_field_apply_bc(app, app->field, emout);
}