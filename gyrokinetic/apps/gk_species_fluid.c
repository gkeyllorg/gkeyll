#include <gkyl_gk_species_priv.h>

void
gk_species_fluid_file_import_init(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gkyl_gyrokinetic_ic_import inp
)
{
  if (inp.type != GKYL_IC_IMPORT_F || inp.enforce_positivity || inp.conf_scale ||
      inp.jacobtot_inv_file_name[0] || inp.jacobvel_file_name[0]) {
    gkyl_exit("Fluid initial-state import supports unmodified configuration-space moments only.");
  }
  enum gkyl_array_rio_status status = gk_species_read(app, gks, inp.file_name);
  if (status != GKYL_ARRAY_RIO_SUCCESS) {
    fprintf(
      stderr, "Unable to import fluid species '%s' from '%s': %s\n", gks->info.name, inp.file_name,
      gkyl_array_rio_status_msg(status)
    );
    gkyl_exit("Fluid initial-state import failed.");
  }
}

// No explicit fluid evolution equations are implemented yet: df/dt = 0.
static double
gk_species_fluid_rhs_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, const struct gkyl_array *fin,
  const struct gkyl_array *fbar_in, struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  gkyl_array_clear(rhs, 0.0);
  return DBL_MAX;
}

// No implicit fluid terms are implemented yet: f_new = f_in.
static double
gk_species_fluid_rhs_implicit_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt
)
{
  gkyl_array_copy(rhs, fin);
  return DBL_MAX;
}

// Release the independently allocated stage arrays and dynamic diagnostics.
static void
gk_species_fluid_release_dynamic(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  gkyl_array_release(gks->f1);
  gkyl_array_release(gks->fnew);
  gk_species_release_integrated_mom(app, gks);
}

// Initialize time stepping and diagnostics for a non-static fluid species.
static void
gk_species_fluid_init_dynamic(gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  gks->f1 = mkarr(app->use_gpu, gks->f->ncomp, gks->f->size);
  gks->fnew = mkarr(app->use_gpu, gks->f->ncomp, gks->f->size);

  gk_species_init_integrated_mom(app, gks);

  gks->rhs_func = gk_species_fluid_rhs_dynamic;
  gks->rhs_implicit_func = gk_species_fluid_rhs_implicit_dynamic;
  gks->bc_func = gk_species_apply_bc_static; // No spatial fluid operators yet.
  gks->release_func = gk_species_fluid_release_dynamic;
  gks->step_f_func = gk_species_step_f_dynamic;
  gks->combine_func = gk_species_combine_dynamic;
  gks->copy_func = gk_species_copy_range_dynamic;
  gks->write_func = gk_species_write_dynamic;
  gks->write_cfl_func = gk_species_write_cfl_disabled;
  gks->write_mom_func = gk_species_write_mom_dynamic;
  gks->calc_integrated_mom_func = gk_species_calc_integrated_mom_dynamic;
  gks->write_integrated_mom_func = gk_species_write_integrated_mom_dynamic;
  gks->calc_L2norm_func = gk_species_calc_L2norm_static;
  gks->write_L2norm_func = gk_species_write_L2norm_static;
  gks->calc_int_mom_dt_func = gk_species_calc_int_mom_dt_disabled;
}

static void
gk_species_fluid_release(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  gks->release_func(app, gks);

  if (gks->info.init_from_file.type == GKYL_IC_IMPORT_NONE) {
    gk_species_projection_release(app, &gks->proj_init);
  }

  gk_species_moment_release(app, &gks->m0);
  gk_species_moment_release(app, &gks->lte.moms);

  for (int m = 0; m < gks->info.num_diag_moments; ++m) {
    gk_species_moment_release(app, &gks->moms[m]);
  }
  gkyl_free(gks->moms);

  gk_species_scaling_release(app, &gks->sca);
  gk_species_lbo_release(app, &gks->lbo);
  gk_species_bgk_release(app, &gks->bgk);

  gkyl_array_release(gks->m0_gyroavg);
  gkyl_array_release(gks->gyro_phi);
  gkyl_array_release(gks->f);
  if (app->use_gpu) {
    gkyl_array_release(gks->f_host);
  }

  gkyl_comm_release(gks->comm);

  gkyl_msgpack_map_elem_release(gks->io_meta_basic_len, gks->io_meta_basic);
  gkyl_msgpack_map_elem_release(gks->io_meta_conf_len, gks->io_meta_conf);
  gkyl_msgpack_map_elem_release(gks->io_meta_phase_len, gks->io_meta_phase);
}

void
gk_species_fluid_init(struct gkyl_gk *inp, gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Fluid species currently have algebraic density closures, not transport equations.
  if (gks->info.collisionless.type || gks->info.flr.type || gks->info.source.source_id ||
      gks->info.source_bgk.source_bgk_id || gks->info.radiation.radiation_id ||
      gks->info.react.num_react || gks->info.react_neut.num_react || gks->info.positivity.type ||
      gks->info.positivity.quasineutrality_rescale || gks->info.damping.type ||
      gks->info.anomalous_diffusion.anomalous_diff_id ||
      gks->info.time_rate_multiplier.num_multipliers || gks->info.scale_with_polarization ||
      gks->info.num_time_rate_diagnostics || gks->info.omegaH_dt_diagnostic ||
      gks->info.write_omega_cfl || gks->info.mapc2p.mapping ||
      gks->info.boundary_flux_diagnostics.num_diag_moments ||
      gks->info.boundary_flux_diagnostics.num_integrated_diag_moments) {
    gkyl_exit("Fluid gyrokinetic species: requested kinetic-only operator or diagnostic.");
  }
  if (gks->info.collisions.collision_id && !gks->info.collisions.not_in_dfdt) {
    gkyl_exit("Fluid gyrokinetic collisions require not_in_dfdt=true.");
  }
  if (gks->info.scaling.type == GKYL_GK_SPECIES_SCALING_RECYCLING_IZ_BALANCE) {
    gkyl_exit("Recycling/ionization balance is not supported for charged fluid species.");
  }
  for (int d = 0; d < GKYL_MAX_VDIM; ++d) {
    if (gks->info.cells[d] || gks->info.lower[d] || gks->info.upper[d]) {
      gkyl_exit("Fluid gyrokinetic species must omit the velocity grid entirely.");
    }
  }

  // Fluid moments always represent parallel and perpendicular velocity coordinates.
  gks->info.vdim = 2;

  gks->num_moments = 4; // rho, Mpar, E, pperp/B (all multiplied by J).

  for (int b = 0; b < 2 * app->cdim; ++b) {
    if (gks->info.bcs[b].type != GKYL_BC_GK_SKIP &&
        gks->info.bcs[b].type != GKYL_BC_GK_SPECIES_PERIODIC) {
      gkyl_exit("Kinetic boundary conditions are not supported for fluid gyrokinetic species.");
    }
  }
  if (gks->info.mass <= 0.0 || gks->info.charge == 0.0) {
    gkyl_exit("Invalid fluid gyrokinetic species mass or charge.");
  }

  gks->release_species_func = gk_species_fluid_release;

  gks->comm = gkyl_comm_acquire(app->comm); // Species comm, same as app/conf comm.

  // Set grid, basis, and ranges (same as app/conf ones).
  gks->grid = app->grid;
  gks->basis = app->basis;
  gks->basis_on_dev = app->basis_on_dev;
  gks->local = app->local;
  gks->local_ext = app->local_ext;
  gks->global = app->global;
  gks->global_ext = app->global_ext;

  int ghost[GKYL_MAX_CDIM] = {1, 1, 1};
  for (int d = 0; d < app->cdim; ++d) {
    gkyl_skin_ghost_ranges(
      &gks->local_lower_skin[d], &gks->local_lower_ghost[d], d, GKYL_LOWER_EDGE, &gks->local_ext,
      ghost
    );
    gkyl_skin_ghost_ranges(
      &gks->local_upper_skin[d], &gks->local_upper_ghost[d], d, GKYL_UPPER_EDGE, &gks->local_ext,
      ghost
    );
    gkyl_skin_ghost_ranges(
      &gks->global_lower_skin[d], &gks->global_lower_ghost[d], d, GKYL_LOWER_EDGE, &gks->global_ext,
      ghost
    );
    gkyl_skin_ghost_ranges(
      &gks->global_upper_skin[d], &gks->global_upper_ghost[d], d, GKYL_UPPER_EDGE, &gks->global_ext,
      ghost
    );
  }

  // Allocte array holding fluid moments.
  gks->f = mkarr(app->use_gpu, gks->num_moments * gks->basis.num_basis, gks->local_ext.volume);
  gks->f_host = app->use_gpu ? mkarr(false, gks->f->ncomp, gks->f->size) : gks->f;

  gk_species_init_io(app, gks);

  // Allocate moment operators to compute M0 and Maxwellian moments.
  gk_species_moment_init(app, gks, &gks->m0, GKYL_F_MOMENT_M0, false);
  gk_species_moment_init(app, gks, &gks->lte.moms, GKYL_F_MOMENT_MAXWELLIAN, false);

  gks->m0_gyroavg = gkyl_array_acquire(gks->m0.marr);
  gks->gyro_phi = gkyl_array_acquire(app->field->phi_smooth);
  gks->gyroaverage = gk_species_gyroaverage_disabled;

  gks->moms = gkyl_malloc(sizeof(struct gk_species_moment) * gks->info.num_diag_moments);
  for (int m = 0; m < gks->info.num_diag_moments; ++m) {
    gk_species_moment_init(app, gks, &gks->moms[m], gks->info.diag_moments[m], false);
  }

  if (gks->info.init_from_file.type == GKYL_IC_IMPORT_NONE) {
    // Set ICs by projecting functions in input file.
    gk_species_projection_init(app, gks, gks->info.projection, &gks->proj_init);
  } else {
    // Set ICs from file.
    gk_species_file_import_init(app, gks, gks->info.init_from_file);
  }

  // Disable collisionless transport.
  gk_species_collisionless_init(app, gks, &gks->collisionless);

  // Disable anomalous diffusion.
  gk_species_anomalous_diff_init(app, gks, &gks->anom_diff);

  // Disable damping.
  gk_species_damping_init(app, gks, &gks->damping);

  // Disable time-rate multipliers.
  gk_species_fdot_multiplier_init(
    app, gks, &gks->fdot_mult, &gks->info.time_rate_multiplier, "fdot_multiplier"
  );

  // Disable the BGK source.
  gk_species_source_bgk_init(app, gks, &gks->bgk_src);

  // Disable positivity correction.
  gk_species_positivity_init(app, gks, &gks->positivity);

  // Disable reactions with charged species.
  gk_species_react_init(app, gks, gks->info.react, &gks->react, true);

  // Disable reactions with neutral species.
  gk_species_react_init(app, gks, gks->info.react_neut, &gks->react_neut, false);

  // Disable boundary-flux calculations.
  gk_species_bflux_init(
    app, gks, &gks->bflux, GK_SPECIES_BFLUX_NONE, (struct gkyl_phase_diagnostics_inp){}
  );

  // Set up density scaling.
  gk_species_scaling_init(app, gks, &gks->sca);

  // Set up LBO collision moments and partners.
  gk_species_lbo_init(app, gks, &gks->lbo);

  // Set up BGK collision moments and partners.
  gk_species_bgk_init(app, gks, &gks->bgk);

  if (gks->info.is_static) {
    gk_species_init_static(inp, app, gks);
  } else {
    gk_species_fluid_init_dynamic(app, gks);
  }
}
