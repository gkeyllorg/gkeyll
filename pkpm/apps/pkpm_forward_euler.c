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
  if (app->field->ext_em_evolve && app->use_explicit_source) {
    pkpm_field_calc_ext_em(app, app->field, tcurr);
  }
  if (app->field->app_current_evolve && app->use_explicit_source) {
    pkpm_field_calc_app_current(app, app->field, tcurr);
  }
  for (int i = 0; i < app->num_species; ++i) {
    if (app->species[i].app_accel_evolve && app->use_explicit_source) {
      pkpm_species_calc_app_accel(app, &app->species[i], tcurr);
    }
  }

  // Compute the magnetic field unit vector b = B/|B|
  pkpm_field_calc_bvar(app, app->field, emin);

  // Every kernel of the stage sees first-harmonic blocks perpendicular to the stage's b: project
  // M0_perp and N_perp of the stage state at the Gauss-Legendre nodes of b (plan 17.11).
  for (int i = 0; i < app->num_species; ++i) {
    pkpm_species_project_transverse(app, &app->species[i], fin[i]);
  }

  // RHS of Maxwell's equations first: it does not depend on the species, and its magnetic part
  // (stored physical in field->dt_em_no_J) is the d_t B of the first-harmonic force terms.
  double dtmin_field = pkpm_field_rhs(app, app->field, stage, emin, emout);

  // Two separate loops over number of species to compute needed auxiliary quantities.
  for (int i = 0; i < app->num_species; ++i) {
    // Compute parallel-kinetic-perpendicular moment (pkpm) variables:
    // the moments [rho, p_par, p_perp, M1], the flow velocity and T_perp/m.
    pkpm_species_calc_pkpm_vars(app, &app->species[i], fin[i], fluidin[i]);
    // The distributions of the stage at the configuration-space nodes, for the collision
    // operators and the temperature source.
    pkpm_species_nodal_dist(app, &app->species[i], fin[i]);
    // compute necessary moments and boundary corrections for collisions
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS) {
      pkpm_species_lbo_moms(app, &app->species[i], &app->species[i].lbo, fin[i]);
    }
    // node values of the fluid state for the computed frequencies and cross-species collisions
    pkpm_species_coll_cross_prim(app, &app->species[i], fluidin[i]);
  }
  for (int i = 0; i < app->num_species; ++i) {
    // rate of the kinetic operator with the partner species: needs the node values (and the
    // self-collision frequencies) of every species, so separate loop over species
    pkpm_species_coll_cross_moms(app, &app->species[i]);
  }

  // Configuration-space face fluxes of the stage for every species (kinetic, fluid and their
  // rates) and the particle loss rates they imply through absorbing walls; then the rates of
  // the boundary-flux sources (which read another species' losses).
  for (int i = 0; i < app->num_species; ++i) {
    struct timespec tm = gkyl_wall_clock();
    gkyl_array_clear(app->species[i].cflrate_f, 0.0);
    gkyl_array_clear(app->species[i].cflrate_fluid, 0.0);
    pkpm_species_face_flux(app, &app->species[i], fin[i]);
    pkpm_species_calc_bflux(app, &app->species[i]);
    app->species[i].fluid_tm += gkyl_time_diff_now_sec(tm);
  }
  for (int i = 0; i < app->num_species; ++i) {
    pkpm_species_source_calc_rate(app, &app->species[i], &app->species[i].src);
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

  // complete update of distribution function; at the nodes the gyrophase mixing damps, the
  // first-harmonic blocks are then replaced by the integrating-factor stage (plan 17.44)
  for (int i = 0; i < app->num_species; ++i) {
    gkyl_array_accumulate(gkyl_array_scale(fout[i], dta), 1.0, fin[i]);
    pkpm_species_harm_exp_advance(app, &app->species[i], stage, dta, fout[i]);
    gkyl_array_accumulate(gkyl_array_scale(fluidout[i], dta), 1.0, fluidin[i]);
    pkpm_species_apply_positivity(app, &app->species[i], fluidout[i]);
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