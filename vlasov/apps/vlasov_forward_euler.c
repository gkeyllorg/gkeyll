#include <gkyl_vlasov_priv.h>

// Take a forward Euler step of the Vlasov-Maxwell system of equations 
// with the suggested time-step dt. Also supports just Maxwell's equations
// and fluid equations (Euler's) with potential Vlasov-fluid coupling. 
// Note: this may not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by stability. 
// The actual time-step and dt_suggested are returned in the status object.
void
vlasov_forward_euler(gkyl_vlasov_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *fluidin[], const struct gkyl_array *emin,
  struct gkyl_array *fout[], struct gkyl_array *fluidout[], struct gkyl_array *emout, 
  struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;

  double dtmin = DBL_MAX;

  // Compute external EM field or applied currents if present and time-dependent.
  // Note: external EM field and  applied currents use proj_on_basis 
  // so does copy to GPU every call if app->use_gpu = true.
  if (app->has_field) {
    if (app->field->app_current_evolve && !app->has_fluid_em_coupling) {
      vm_field_calc_app_current(app, app->field, tcurr); 
    }
    if (app->field->ext_em_evolve) {
      vm_field_calc_ext_em(app, app->field, tcurr);
    }
    if (app->field->ext_pot_evolve) {
      vp_field_calc_ext_pot(app, app->field, tcurr);
    }
  }
  // Compute applied acceleration if if present and time-dependent.
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].app_accel_evolve) {
      vm_species_calc_app_accel(app, &app->species[i], tcurr);
    }
  }

  // compute necessary moments and boundary corrections for collisions
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      vm_species_lbo_moms(app, &app->species[i], &app->species[i].lbo, fin[i]);
    }
    else if (app->species[i].collision_id == GKYL_BGK_COLLISIONS && !app->has_implicit_coll_scheme) {
      vm_species_bgk_moms(app, &app->species[i], 
        &app->species[i].bgk, fin[i]);
    }
    else if (app->species[i].collision_id == GKYL_FPO_COLLISIONS) {
      vm_species_fpo_drag_diff_coeffs(app, &app->species[i], 
        &app->species[i].fpo, fin[i]);
    }
  }

  // compute necessary moments for cross-species collisions
  // needs to be done after self-collisions moments, so separate loop over species
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS
      && app->species[i].lbo.num_cross_collisions) {
      vm_species_lbo_cross_moms(app, &app->species[i], &app->species[i].lbo, fin[i]);
    }
  }

  // Compute primitive moments for fluid species evolution
  for (int i=0; i<app->num_fluid_species; ++i) 
    vm_fluid_species_prim_vars(app, &app->fluid_species[i], fluidin[i]);

  // compute RHS of Vlasov equations
  for (int i=0; i<app->num_species; ++i) {
    double dt1 = vm_species_rhs(app, &app->species[i], fin[i], emin, fout[i]);
    dtmin = fmin(dtmin, dt1);
  }
  for (int i=0; i<app->num_fluid_species; ++i) {
    double dt1 = vm_fluid_species_rhs(app, &app->fluid_species[i], fluidin[i], emin, fluidout[i]);
    dtmin = fmin(dtmin, dt1);
  }
  // compute source term
  // done here as the RHS update for all species should be complete before
  // bflux calculation of the source species
  for (int i=0; i<app->num_species; ++i) {
    if (app->species[i].source_id) {
      vm_species_source_rhs(app, &app->species[i], &app->species[i].src, fin, fout);
    }
  }
  for (int i=0; i<app->num_fluid_species; ++i) {
    if (app->fluid_species[i].source_id) {
      vm_fluid_species_source_rhs(app, &app->fluid_species[i], &app->fluid_species[i].src, fluidin, fluidout);
    }
  }
  // compute RHS of Maxwell equations
  if (app->has_field) {
    if (app->field->field_id == GKYL_FIELD_E_B) {
      double dt1 = vm_field_rhs(app, app->field, emin, emout);
      dtmin = fmin(dtmin, dt1);
    }
  }

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
  }

  // complete update of fluid species
  for (int i=0; i<app->num_fluid_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fluidout[i], dta), 1.0, fluidin[i]);
  }

  if (app->has_field) {
    if (app->field->field_id == GKYL_FIELD_E_B) {
      struct timespec wst = gkyl_wall_clock();

      // (can't accumulate current when field is static)
      if (!app->field->info.is_static) {
        // accumulate current contribution from kinetic species to electric field terms
        vm_field_accumulate_current(app, fin, fluidin, emout);
        app->stat.current_tm += gkyl_time_diff_now_sec(wst);
      }

      // complete update of field (even when field is static, it is
      // safest to do this accumulate as it ensure emout = emin)
      gkyl_array_accumulate(gkyl_array_scale(emout, dta), 1.0, emin);
    }
  }
}
