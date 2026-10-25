#include <gkyl_app.h>
#include <gkyl_flash_priv.h>

// Take a single forward Euler timestep of the hybrid Euler-PKPM equation system
// with suggested timestep dt0. The actual timestep taken may be smaller than the
// suggested timestep, and both the suggested and actual timesteps are returned
// in the status structure.
void
flash_forward_euler(gkyl_flash_app *app, double tcurr, double dt0,
  const struct gkyl_array *fin[], const struct gkyl_array *fv_fluidin[],
  struct gkyl_array *fout[], struct gkyl_array *fv_fluidout[],
  struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;
  struct timespec tm = gkyl_wall_clock();

  double dtmin = DBL_MAX;

  flash_field_calc_bvar(app, app->field, app->field->em);

  // Compute required variables for PKPM update, interpolate bulk velocity onto kinetic grid,
  // and calculate flux contributions for Euler equations.
  for (int i=0; i<app->num_species; ++i) {
    // Compute PKPM moments, [rho, p_par, p_perp].
    flash_species_calc_pkpm_moms(app, &app->species[i], fin[i]);

    // Interpolate fluid momentum density vector onto kinetic grid.
    flash_species_interp_rhou_to_dg(app, &app->species[i], &app->fluid_species[i], fv_fluidin[i]);

    // Compute coupling moments, self-consistent pressure
    // force (div(p_par b_hat)), and primitive variables
    flash_species_calc_pkpm_vars(app, &app->species[i], fin[i]);

    // Compute remaining PKPM variables: acceleration variables in kinetic equation and 
    // source distribution functions for Laguerre couplings.
    flash_species_calc_pkpm_update_vars(app, &app->species[i], fin[i]);

    // Compute flux contributions to Euler system from kinetic closure.
    flash_species_calc_coupling_fluxes(app, &app->species[i], fin[i]);
  }

  // Compute RHS of PKPM kinetic system.
  for (int i=0; i<app->num_species; ++i) {
    double dt1 = flash_species_rhs(app, &app->species[i], fin[i], fout[i]);
    dtmin = fmin(dtmin, dt1);
  }

  // Advance Euler system with wave propagation algorithm,
  // storing fluid timestep to determine if CFL was violated.
  double dt_flu = DBL_MAX;
  for (int i=0; i<app->num_species; ++i) {
    double dt1 = flash_fluid_species_update(app, tcurr, dt0, &app->species[i], &app->fluid_species[i], fv_fluidin[i], fv_fluidout[i]);
    dt_flu = fmin(dt_flu, dt1);
  }
  dtmin = fmin(dtmin, dt_flu);

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  // Not applied if the fluid update failed.
  double dt_rel_diff = (dt0-dtmin)/dt0;
  if (dt_flu >= dt0 && dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt0;

  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt0 < dtmin ? dt0 : dtmin;
  st->dt_suggested = dtmin;

  // If wave propagation update failed earlier due to violating CFL, re-run with the proper timestep.
  if (dta < dt0) {
    for (int i=0; i<app->num_species; ++i)
      flash_fluid_species_update(app, tcurr, dta, &app->species[i], &app->fluid_species[i], fv_fluidin[i], fv_fluidout[i]);
  }

  // Complete forward Euler update of distribution functions,
  // and apply BCs to distribution functions and Euler variables.
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dta), 1.0, fin[i]);

    flash_species_apply_bc(app, &app->species[i], fout[i]);
    flash_fluid_species_apply_bc(app, &app->fluid_species[i], fv_fluidout[i]);
  }
  app->stat.fwd_euler_tm += gkyl_time_diff_now_sec(tm);
}

// Take a single forward Euler timestep of the hybrid Euler-PKPM equation system
// with suggested timestep dt0. The actual timestep taken may be smaller than the
// suggested timestep, and both the suggested and actual timesteps are returned
// in the status structure.
void
flash_vlasov_forward_euler(gkyl_flash_app *app, double tcurr, double dt0,
  const struct gkyl_array *fin[], struct gkyl_array *fout[],
  struct gkyl_update_status *st)
{
  app->stat.nfeuler += 1;
  struct timespec tm = gkyl_wall_clock();

  double dtmin = DBL_MAX;

  flash_field_calc_bvar(app, app->field, app->field->em);

  // Compute required variables for PKPM update, interpolate bulk velocity onto kinetic grid,
  // and calculate flux contributions for Euler equations.
  for (int i=0; i<app->num_species; ++i) {
    // Compute PKPM moments, [rho, p_par, p_perp].
    flash_species_calc_pkpm_moms(app, &app->species[i], fin[i]);

    // Compute coupling moments, self-consistent pressure
    // force (div(p_par b_hat)), and primitive variables
    flash_species_calc_pkpm_vars(app, &app->species[i], fin[i]);

    // Compute remaining PKPM variables: acceleration variables in kinetic equation and 
    // source distribution functions for Laguerre couplings.
    flash_species_calc_pkpm_update_vars(app, &app->species[i], fin[i]);

    // // Compute flux contributions to Euler system from kinetic closure.
    // flash_species_calc_coupling_fluxes(app, &app->species[i], fin[i]);
  }

  // Compute RHS of PKPM kinetic system.
  for (int i=0; i<app->num_species; ++i) {
    double dt1 = flash_species_rhs(app, &app->species[i], fin[i], fout[i]);
    dtmin = fmin(dtmin, dt1);
  }

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt0-dtmin)/dt0;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt0;

  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt0 < dtmin ? dt0 : dtmin;
  st->dt_suggested = dtmin;

  // Complete forward Euler update of distribution functions,
  // and apply BCs to distribution functions and Euler variables
  for (int i=0; i<app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dta), 1.0, fin[i]);

    flash_species_apply_bc(app, &app->species[i], fout[i]);
  }
  app->stat.fwd_euler_tm += gkyl_time_diff_now_sec(tm);
}
