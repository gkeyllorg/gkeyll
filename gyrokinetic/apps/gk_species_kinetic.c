#include <gkyl_gk_species_priv.h>

double
gk_species_omegaH_frequency(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, const struct gkyl_array *fin
)
{
  // Compute omega_H without a CFL factor; the time-step constraint uses
  // dt_omegaH = omegaH_CFL/omega_H.
  // Each species computes its own omega_H as:
  //   omega_H = q_e*sqrt(jacobgeo*n_{s0}/m_s) * omegaH_gf
  // where
  //   - n_{s0} is either a reference, average or max density.
  //   - omegaH_gf = (cmag/(jacobgeo*B^_\parallel))*kpar_max /
  //                 min(sqrt(k_x^2*eps_xx+k_x*k_y*eps_xy+k_y^2*eps_yy+)).
  // and k_x,k_y,k_par are wavenumbers in computational space, and eps_ij is
  // the polarization weight in our field equation.

  if (!(app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN || app->field->has_adiabatic_species)) {
    // Obtain the maximum density (using cell centers).
    gk_species_moment_calc(app, &gks->m0, gks->local, app->local, fin);
    gkyl_array_reduce_range(gks->m0_max, gks->m0.marr, GKYL_MAX, &app->local);

    double m0_max[1];
    if (app->use_gpu) {
      gkyl_cu_memcpy(m0_max, gks->m0_max, sizeof(double), GKYL_CU_MEMCPY_D2H);
    } else {
      m0_max[0] = gks->m0_max[0];
    }
    m0_max[0] *= 1.0 / pow(sqrt(2.0), app->cdim);

    double time_dilation_scale_const =
      gk_fdot_multiplier_get_time_dilation_scale_const(app, &gks->fdot_mult);

    double omegaH = fabs(gks->info.charge) * sqrt(GKYL_MAX2(0.0, m0_max[0]) / gks->info.mass) *
                    app->omegaH_gf * time_dilation_scale_const * gks->collisionless.scale_fac;

    return omegaH;
  } else {
    return 0.0;
  }
}

static double
gk_species_kinetic_omegaH_dt(
  gkyl_gyrokinetic_app *app, struct gk_species *gks, const struct gkyl_array *fin
)
{
  double omegaH = gk_species_omegaH_frequency(app, gks, fin);
  return omegaH > 1e-20 ? app->cfl_omegaH / omegaH : DBL_MAX;
}

static double
gk_species_kinetic_rhs_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  // Gyroaverage the potential if needed.
  species->gyroaverage(app, species, app->field->phi_smooth, species->gyro_phi);

  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Collisionless terms.
  gk_species_collisionless_rhs(app, species, &species->collisionless, fin, rhs);

  // Damping term.
  gk_species_damping_advance(
    app, species, &species->damping, app->field->phi_smooth, fin, species->lte.f_lte, rhs,
    species->cflrate
  );

  // LBO Collisions.
  gk_species_lbo_rhs(app, species, &species->lbo, fin, rhs);

  // BGK collisions.
  gk_species_bgk_rhs(app, species, &species->bgk, fin, rhs);

  // Anomalous diffusion.
  gk_species_anomalous_diff_rhs(app, species, &species->anom_diff, fin, rhs);

  // Line radiation.
  gk_species_radiation_rhs(app, species, &species->rad, fin, rhs);

  // Reactions with charged species.
  gk_species_react_rhs(app, species, &species->react, fin, rhs);

  // Reactions with neutral species.
  gk_species_react_rhs(app, species, &species->react_neut, fin, rhs);

  // BGK source.
  gk_species_source_bgk_rhs(app, species, &species->bgk_src, fin, rhs);

  // Compute and store (in the ghost cell of rhs and in an array in bflux) the boundary fluxes.
  gk_species_bflux_rhs(app, &species->bflux, fin, rhs);

  // Compute moments of the boundary fluxes.
  gk_species_bflux_calc_moms(app, &species->bflux, rhs, bflux_moms);

  // Multiply CFL rate by the df/dt multiplier.
  gk_species_fdot_multiplier_advance_times_cfl(
    app, species, &species->fdot_mult, app->field->phi_smooth, fin, species->cflrate
  );

  // Reduce the CFL frequency and compute stable dt needed by this species.
  app->stat.n_species_omega_cfl += 1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    omega_cfl_ho[0] = species->omega_cfl[0];
  }

  double dt_out = app->cfl / omega_cfl_ho[0];

  // Apply omega_H constraint on dt.
  species->dt_omegaH = gk_species_kinetic_omegaH_dt(app, species, fin);
  dt_out = fmin(dt_out, species->dt_omegaH);

  species->dt_cfl_global_ho = dt_out;

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return dt_out;
}

static double
gk_species_kinetic_rhs_implicit_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_species *species, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt
)
{
  double omega_cfl = 1 / DBL_MAX;
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Compute implicit update and update rhs to new time step
  gk_species_bgk_rhs_implicit(app, species, &species->bgk, fin, dt, rhs);

  gkyl_array_accumulate(gkyl_array_scale(rhs, dt), 1.0, fin);

  app->stat.n_species_omega_cfl += 1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    omega_cfl_ho[0] = species->omega_cfl[0];
  }
  omega_cfl = omega_cfl_ho[0];

  app->stat.species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl / omega_cfl;
}

static void
gk_species_kinetic_apply_bc_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_species *species, double tm, struct gkyl_array *f
)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = species->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(
    species->comm, &species->local, &species->local_ext, num_periodic_dir, species->periodic_dirs, f
  );

  for (int d = 0; d < cdim; ++d) {
    if (species->bc_is_np[d]) {
      switch (species->lower_bc[d].type) {
        case GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING:
          gk_species_phi_wall_advance(app, &species->phi_wall_lo, tm);
          gkyl_bc_sheath_gyrokinetic_advance(
            species->bc_sheath_lo, app->field->phi_smooth, species->phi_wall_lo.phi, NULL, NULL,
            NULL, NULL, f, &app->local
          );
          break;
        case GKYL_BC_GK_SPECIES_SHEATH_SURROGATE:
          gk_species_phi_wall_advance(app, &species->phi_wall_lo, tm);
          gk_species_moment_calc(app, &species->sheath_moms, species->local, app->local, f);
          gk_species_moment_diag_jacobgeo_div(
            app, &species->sheath_moms, species->sheath_moms.marr, species->sheath_moms.marr
          );
          gkyl_array_set_offset_range(
            species->dens_sheath, 1.0, species->sheath_moms.marr, 0, &app->local
          );
          gkyl_array_set_offset_range(
            species->temp_sheath, species->info.mass, species->sheath_moms.marr,
            2 * app->basis.num_basis, &app->local
          );
          gkyl_bc_sheath_gyrokinetic_advance(
            species->bc_sheath_lo, app->field->phi_smooth, species->phi_wall_lo.phi,
            species->dens_sheath, species->temp_sheath, app->gk_geom->geo_int.bmag,
            app->gk_geom->geo_surf[d].bimpactangle, f, &app->local
          );
          break;
        case GKYL_BC_GK_SPECIES_TWISTSHIFT:
          gkyl_bc_twistshift_advance(species->bc_ts_lo, f, f);
          break;
        case GKYL_BC_GK_SPECIES_COPY:
        case GKYL_BC_GK_SPECIES_REFLECT:
        case GKYL_BC_GK_SPECIES_ABSORB:
          gkyl_bc_basic_gyrokinetic_advance(species->bc_lo[d], species->bc_buffer, f);
          break;
        case GKYL_BC_GK_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_gyrokinetic_advance(species->bc_lo[d], species->bc_buffer_lo_fixed, f);
          break;
        case GKYL_BC_GK_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }

      switch (species->upper_bc[d].type) {
        case GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING:
          gk_species_phi_wall_advance(app, &species->phi_wall_up, tm);
          gkyl_bc_sheath_gyrokinetic_advance(
            species->bc_sheath_up, app->field->phi_smooth, species->phi_wall_up.phi, NULL, NULL,
            NULL, NULL, f, &app->local
          );
          break;
        case GKYL_BC_GK_SPECIES_SHEATH_SURROGATE:
          gk_species_phi_wall_advance(app, &species->phi_wall_up, tm);
          gk_species_moment_calc(app, &species->sheath_moms, species->local, app->local, f);
          gk_species_moment_diag_jacobgeo_div(
            app, &species->sheath_moms, species->sheath_moms.marr, species->sheath_moms.marr
          );
          gkyl_array_set_offset_range(
            species->dens_sheath, 1.0, species->sheath_moms.marr, 0, &app->local
          );
          gkyl_array_set_offset_range(
            species->temp_sheath, species->info.mass, species->sheath_moms.marr,
            2 * app->basis.num_basis, &app->local
          );
          gkyl_bc_sheath_gyrokinetic_advance(
            species->bc_sheath_up, app->field->phi_smooth, species->phi_wall_up.phi,
            species->dens_sheath, species->temp_sheath, app->gk_geom->geo_int.bmag,
            app->gk_geom->geo_surf[d].bimpactangle, f, &app->local
          );
          break;
        case GKYL_BC_GK_SPECIES_TWISTSHIFT:
          gkyl_bc_twistshift_advance(species->bc_ts_up, f, f);
          break;
        case GKYL_BC_GK_SPECIES_COPY:
        case GKYL_BC_GK_SPECIES_REFLECT:
        case GKYL_BC_GK_SPECIES_ABSORB:
          gkyl_bc_basic_gyrokinetic_advance(species->bc_up[d], species->bc_buffer, f);
          break;
        case GKYL_BC_GK_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_gyrokinetic_advance(species->bc_up[d], species->bc_buffer_up_fixed, f);
          break;
        case GKYL_BC_GK_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }
    }
  }

  if (app->gk_geom->has_LCFS && app->cdim == 3) {
    // Apply twistshift.
    gkyl_bc_twistshift_advance(species->bc_ts_lo, f, f);
    gkyl_bc_twistshift_advance(species->bc_ts_up, f, f);
  }

  gkyl_comm_array_sync(species->comm, &species->local, &species->local_ext, f);

  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);
}

static void
gk_species_kinetic_release_dynamic(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  // Release various arrays and objects for a dynamic species.
  gkyl_array_release(gks->f1);
  gkyl_array_release(gks->fnew);
  gkyl_array_release(gks->bc_buffer);
  gkyl_array_release(gks->bc_buffer_lo_fixed);
  gkyl_array_release(gks->bc_buffer_up_fixed);

  if (gks->info.write_omega_cfl) {
    gkyl_array_release(gks->cflrate_ho);
  }

  int par_dir = app->cdim - 1; // Parallel direction index.

  // Copy BCs are allocated by default. Need to free.
  bool release_surr_aux_var = false;
  for (int d = 0; d < app->cdim; ++d) {
    if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING) {
      gkyl_bc_sheath_gyrokinetic_release(gks->bc_sheath_lo);
      gk_species_phi_wall_release(app, &gks->phi_wall_lo);
    } else if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE) {
      gkyl_bc_sheath_gyrokinetic_release(gks->bc_sheath_lo);
      gk_species_phi_wall_release(app, &gks->phi_wall_lo);
      release_surr_aux_var = true;
    } else if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_TWISTSHIFT) {
      gkyl_bc_twistshift_release(gks->bc_ts_lo);
    } else if ((gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
               (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
               (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
               (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
      gkyl_bc_basic_gyrokinetic_release(gks->bc_lo[d]);
    }

    if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING) {
      gkyl_bc_sheath_gyrokinetic_release(gks->bc_sheath_up);
      gk_species_phi_wall_release(app, &gks->phi_wall_up);
    } else if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE) {
      gkyl_bc_sheath_gyrokinetic_release(gks->bc_sheath_up);
      gk_species_phi_wall_release(app, &gks->phi_wall_up);
      release_surr_aux_var = true;
    } else if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_TWISTSHIFT) {
      gkyl_bc_twistshift_release(gks->bc_ts_up);
    } else if ((gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
               (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
               (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
               (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
      gkyl_bc_basic_gyrokinetic_release(gks->bc_up[d]);
    }
  }
  if (release_surr_aux_var) {
    gk_species_moment_release(app, &gks->sheath_moms);
    gkyl_array_release(gks->dens_sheath);
    gkyl_array_release(gks->temp_sheath);
  }

  if (app->gk_geom->has_LCFS && app->cdim == 3) {
    // Free twishift memory.
    gkyl_bc_twistshift_release(gks->bc_ts_lo);
    gkyl_bc_twistshift_release(gks->bc_ts_up);
  }

  if (app->use_gpu) {
    gkyl_cu_free(gks->omega_cfl);
    gkyl_cu_free(gks->m0_max);
  } else {
    gkyl_free(gks->omega_cfl);
    gkyl_free(gks->m0_max);
  }

  gk_species_release_integrated_mom(app, gks);

  if (gks->info.frequency_diagnostics) {
    gk_species_frequencies_release(app, gks);
  }

  // Release L2 norm memory.
  gkyl_array_integrate_release(gks->integ_wfsq_op);
  gkyl_dynvec_release(gks->L2norm);
  if (app->use_gpu) {
    gkyl_cu_free(gks->L2norm_local);
    gkyl_cu_free(gks->L2norm_global);
  } else {
    gkyl_free(gks->L2norm_local);
    gkyl_free(gks->L2norm_global);
  }
}

static void
gk_species_kinetic_init_dynamic(
  struct gkyl_gk *gk_app_inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks
)
{
  int cdim = app->cdim, vdim = gks->info.vdim;
  int pdim = cdim + vdim;

  int num_ghost[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    num_ghost[d] = 1;
  }
  for (int d = 0; d < vdim; ++d) {
    num_ghost[cdim + d] = 0; // No ghost-cells in velocity space.
  }
  gks->dt_omegaH = DBL_MIN;

  // Allocate distribution function arrays.
  gks->f1 = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);
  gks->fnew = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);

  if (app->use_gpu) {
    gks->omega_cfl = gkyl_cu_malloc(sizeof(double));
    gks->m0_max = gkyl_cu_malloc(app->basis.num_basis * sizeof(double));
  } else {
    gks->omega_cfl = gkyl_malloc(sizeof(double));
    gks->m0_max = gkyl_malloc(app->basis.num_basis * sizeof(double));
  }

  gk_species_init_integrated_mom(app, gks);

  if (gks->info.frequency_diagnostics) {
    gk_species_frequencies_init(app, gks);
  }

  // Objects for L2 norm diagnostic.
  gks->integ_wfsq_op = gkyl_array_integrate_new(
    &gks->grid, &gks->basis, 1, GKYL_ARRAY_INTEGRATE_OP_SQ_WEIGHTED, app->use_gpu
  );
  if (app->use_gpu) {
    gks->L2norm_local = gkyl_cu_malloc(sizeof(double));
    gks->L2norm_global = gkyl_cu_malloc(sizeof(double));
  } else {
    gks->L2norm_local = gkyl_malloc(sizeof(double));
    gks->L2norm_global = gkyl_malloc(sizeof(double));
  }
  gks->L2norm = gkyl_dynvec_new(GKYL_DOUBLE, 1); // L2 norm.
  gks->is_first_L2norm_write_call = true;

  int par_dir = app->cdim - 1; // Parallel direction index.
  if (gk_app_inp->geometry.has_LCFS ||
      (gks->lower_bc[par_dir].type == GKYL_BC_GK_SPECIES_TWISTSHIFT ||
       gks->upper_bc[par_dir].type == GKYL_BC_GK_SPECIES_TWISTSHIFT)) {
    // Make the parallel direction periodic so that we sync before applying TS BC.
    gks->periodic_dirs[gks->num_periodic_dir] = par_dir;
    gks->num_periodic_dir += 1;
  }

  // Allocate buffer needed for BCs.
  bool need_bc_buffer = false;
  bool need_bc_buffer_lo_fixed = false;
  bool need_bc_buffer_up_fixed = false;
  for (int d = 0; d < cdim; ++d) {
    if (gks->bc_is_np[d]) {
      if ((gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
          (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
          (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
          (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
          (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
          (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB)) {
        need_bc_buffer = true;
      }
      if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) {
        need_bc_buffer_lo_fixed = true;
      }
      if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) {
        need_bc_buffer_up_fixed = true;
      }
    }
  }

  long buff_sz = 1;
  for (int dir = 0; dir < cdim; ++dir) {
    long vol = GKYL_MAX2(gks->local_lower_skin[dir].volume, gks->local_upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  gks->bc_buffer = need_bc_buffer ? mkarr(app->use_gpu, gks->basis.num_basis, buff_sz) :
                                    mkarr(app->use_gpu, 1, 1);
  // Buffer arrays for fixed function BCs on distribution function.
  gks->bc_buffer_lo_fixed = need_bc_buffer_lo_fixed ?
                              mkarr(app->use_gpu, gks->basis.num_basis, buff_sz) :
                              mkarr(app->use_gpu, 1, 1);
  gks->bc_buffer_up_fixed = need_bc_buffer_up_fixed ?
                              mkarr(app->use_gpu, gks->basis.num_basis, buff_sz) :
                              mkarr(app->use_gpu, 1, 1);

  if (gks->lower_bc[par_dir].type == GKYL_BC_GK_SPECIES_TWISTSHIFT &&
      gks->upper_bc[par_dir].type == GKYL_BC_GK_SPECIES_TWISTSHIFT) {
    // Local range extended in the BC dir, on the coarse grid.
    int lower_bcdir_ext[pdim], upper_bcdir_ext[pdim];
    for (int i = 0; i < pdim; i++) {
      lower_bcdir_ext[i] = gks->local.lower[i];
      upper_bcdir_ext[i] = gks->local.upper[i];
    }
    lower_bcdir_ext[par_dir] = gks->local_ext.lower[par_dir];
    upper_bcdir_ext[par_dir] = gks->local_ext.upper[par_dir];
    gkyl_sub_range_init(
      &gks->bc_ts_local_par_ext, &gks->local_ext, lower_bcdir_ext, upper_bcdir_ext
    );
  }

  for (int d = 0; d < cdim; ++d) {
    // Lower BC.
    if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING ||
        gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE) {
      struct gkyl_range *sol_skin = gk_app_inp->geometry.has_LCFS ? &gks->local_lower_skin_par_sol :
                                                                    &gks->local_lower_skin[d];
      struct gkyl_range *sol_ghost = gk_app_inp->geometry.has_LCFS ?
                                       &gks->local_lower_ghost_par_sol :
                                       &gks->local_lower_ghost[d];
      gks->bc_sheath_lo = gkyl_bc_sheath_gyrokinetic_new(
        d, GKYL_LOWER_EDGE, gks->basis_on_dev, sol_skin, sol_ghost, gks->vel_map, cdim,
        2.0 * (gks->info.charge / gks->info.mass), gks->lower_bc[d].type, gks->lower_bc[d].aux_str,
        &gks->grid, &gks->local, gks->basis.poly_order, app->use_gpu
      );
      gk_species_phi_wall_init(app, &gks->lower_bc[d], &gks->phi_wall_lo);
    } else if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_TWISTSHIFT) {
      assert(cdim == 3);
      struct gkyl_bc_twistshift_inp tsinp = {
        .bc_dir = d,
        .shift_dir = 1, // y shift.
        .shear_dir = 0, // shift varies with x.
        .edge = GKYL_LOWER_EDGE,
        .cdim = cdim,
        .bcdir_ext_update_r = &gks->bc_ts_local_par_ext,
        .num_ghost = num_ghost,
        .basis = &gks->basis,
        .grid = &gks->grid,
        .use_gpu = app->use_gpu,
        .upsample_factor = app->core_parallel_bcs.ts_upsample_factor,
        .filter_half_width = app->core_parallel_bcs.ts_filter_half_width,
        .filter_cutoff_wavelength = app->core_parallel_bcs.ts_filter_cutoff_wavelength,
      };
      if (app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK) {
        tsinp.shift_dg = app->delta_ts_x_lo;
      } else {
        tsinp.shift_func = app->gk_geom->lower_shift_func;
        tsinp.shift_func_ctx = app->gk_geom->lower_shift_ctx;
      }

      gks->bc_ts_lo = gkyl_bc_twistshift_inew(&tsinp);

    } else if ((gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
               (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
               (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
               (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
      if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) {
        assert(
          d == cdim - 1
        ); // MF 2025/09/25: Reflecting BCs only allowed in the parallel direction.
      }

      gks->bc_lo[d] = gkyl_bc_basic_gyrokinetic_new(
        d, GKYL_LOWER_EDGE, gks->lower_bc[d].type, gks->basis_on_dev, &gks->local_lower_skin[d],
        &gks->local_lower_ghost[d], gks->f->ncomp, app->cdim, app->use_gpu
      );

      if (gks->lower_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) {
        // Project distribution desired in the ghost cell.
        struct gk_proj gk_proj_bc_lo;
        gk_species_projection_init(app, gks, gks->lower_bc[d].projection, &gk_proj_bc_lo);
        gk_species_projection_calc(app, gks, &gk_proj_bc_lo, gks->f1, 0.0); // Temporarily use f1.
        // Fill ghost cell with the boundary value.
        struct gkyl_bc_basic_gyrokinetic *gfss_bc_op = gkyl_bc_basic_gyrokinetic_new(
          d, GKYL_LOWER_EDGE, GKYL_BC_GK_SPECIES_BOUNDARY_VALUE, gks->basis_on_dev,
          &gks->local_lower_skin[d], &gks->local_lower_ghost[d], gks->f->ncomp, app->cdim,
          app->use_gpu
        );
        gkyl_bc_basic_gyrokinetic_advance(gfss_bc_op, gks->bc_buffer_lo_fixed, gks->f1);
        gkyl_bc_basic_gyrokinetic_release(gfss_bc_op);
        // Copy ghost range into buffer.
        gkyl_bc_basic_gyrokinetic_buffer_fixed_func(gks->bc_lo[d], gks->bc_buffer_lo_fixed, gks->f1);
        gkyl_array_clear(gks->f1, 0.0);
        gk_species_projection_release(app, &gk_proj_bc_lo);
      }
    }

    // Upper BC.
    if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH_CONDUCTING ||
        gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE) {
      struct gkyl_range *sol_skin = gk_app_inp->geometry.has_LCFS ? &gks->local_upper_skin_par_sol :
                                                                    &gks->local_upper_skin[d];
      struct gkyl_range *sol_ghost = gk_app_inp->geometry.has_LCFS ?
                                       &gks->local_upper_ghost_par_sol :
                                       &gks->local_upper_ghost[d];
      gks->bc_sheath_up = gkyl_bc_sheath_gyrokinetic_new(
        d, GKYL_UPPER_EDGE, gks->basis_on_dev, sol_skin, sol_ghost, gks->vel_map, cdim,
        2.0 * (gks->info.charge / gks->info.mass), gks->upper_bc[d].type, gks->upper_bc[d].aux_str,
        &gks->grid, &gks->local, gks->basis.poly_order, app->use_gpu
      );
      gk_species_phi_wall_init(app, &gks->upper_bc[d], &gks->phi_wall_up);
    } else if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_TWISTSHIFT) {
      assert(cdim == 3);
      struct gkyl_bc_twistshift_inp tsinp = {
        .bc_dir = d,
        .shift_dir = 1, // y shift.
        .shear_dir = 0, // shift varies with x.
        .edge = GKYL_UPPER_EDGE,
        .cdim = cdim,
        .bcdir_ext_update_r = &gks->bc_ts_local_par_ext,
        .num_ghost = num_ghost,
        .basis = &gks->basis,
        .grid = &gks->grid,
        .use_gpu = app->use_gpu,
        .upsample_factor = app->core_parallel_bcs.ts_upsample_factor,
        .filter_half_width = app->core_parallel_bcs.ts_filter_half_width,
        .filter_cutoff_wavelength = app->core_parallel_bcs.ts_filter_cutoff_wavelength,
      };
      if (app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK) {
        tsinp.shift_dg = app->delta_ts_x_up;
      } else {
        tsinp.shift_func = app->gk_geom->upper_shift_func;
        tsinp.shift_func_ctx = app->gk_geom->upper_shift_ctx;
      }

      gks->bc_ts_up = gkyl_bc_twistshift_inew(&tsinp);
    } else if ((gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
               (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
               (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
               (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
      if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) {
        assert(
          d == cdim - 1
        ); // MF 2025/09/25: Reflecting BCs only allowed in the parallel direction.
      }

      gks->bc_up[d] = gkyl_bc_basic_gyrokinetic_new(
        d, GKYL_UPPER_EDGE, gks->upper_bc[d].type, gks->basis_on_dev, &gks->local_upper_skin[d],
        &gks->local_upper_ghost[d], gks->f->ncomp, app->cdim, app->use_gpu
      );

      if (gks->upper_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) {
        // Project distribution desired in the ghost cell.
        struct gk_proj gk_proj_bc_up;
        gk_species_projection_init(app, gks, gks->upper_bc[d].projection, &gk_proj_bc_up);
        gk_species_projection_calc(app, gks, &gk_proj_bc_up, gks->f1, 0.0); // Temporarily use f1.
        // Fill ghost cell with the boundary value.
        struct gkyl_bc_basic_gyrokinetic *gfss_bc_op = gkyl_bc_basic_gyrokinetic_new(
          d, GKYL_UPPER_EDGE, GKYL_BC_GK_SPECIES_BOUNDARY_VALUE, gks->basis_on_dev,
          &gks->local_upper_skin[d], &gks->local_upper_ghost[d], gks->f->ncomp, app->cdim,
          app->use_gpu
        );
        gkyl_bc_basic_gyrokinetic_advance(gfss_bc_op, gks->bc_buffer_up_fixed, gks->f1);
        gkyl_bc_basic_gyrokinetic_release(gfss_bc_op);
        // Copy ghost range into buffer.
        gkyl_bc_basic_gyrokinetic_buffer_fixed_func(gks->bc_up[d], gks->bc_buffer_up_fixed, gks->f1);
        gkyl_array_clear(gks->f1, 0.0);
        gk_species_projection_release(app, &gk_proj_bc_up);
      }
    }
  }

  if (gk_app_inp->geometry.has_LCFS && app->cdim == 3) {
    // Allocate a twistshift operator in the core.
    struct gkyl_bc_twistshift_inp tsinp_lo = {
      .bc_dir = par_dir,
      .shift_dir = 1, // y shift.
      .shear_dir = 0, // shift varies with x.
      .edge = GKYL_LOWER_EDGE,
      .cdim = cdim,
      .bcdir_ext_update_r = &gks->local_par_ext_core,
      .num_ghost = num_ghost,
      .basis = &gks->basis,
      .grid = &gks->grid,
      .use_gpu = app->use_gpu,
      .upsample_factor = app->core_parallel_bcs.ts_upsample_factor,
      .filter_half_width = app->core_parallel_bcs.ts_filter_half_width,
      .filter_cutoff_wavelength = app->core_parallel_bcs.ts_filter_cutoff_wavelength,
    };
    if (app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK) {
      tsinp_lo.shift_dg = app->delta_ts_x_lo;
    } else {
      tsinp_lo.shift_func = app->gk_geom->lower_shift_func;
      tsinp_lo.shift_func_ctx = app->gk_geom->lower_shift_ctx;
    }
    gks->bc_ts_lo = gkyl_bc_twistshift_inew(&tsinp_lo);

    struct gkyl_bc_twistshift_inp tsinp_up = {
      .bc_dir = par_dir,
      .shift_dir = 1, // y shift.
      .shear_dir = 0, // shift varies with x.
      .edge = GKYL_UPPER_EDGE,
      .cdim = cdim,
      .bcdir_ext_update_r = &gks->local_par_ext_core,
      .num_ghost = num_ghost,
      .basis = &gks->basis,
      .grid = &gks->grid,
      .use_gpu = app->use_gpu,
      .upsample_factor = app->core_parallel_bcs.ts_upsample_factor,
      .filter_half_width = app->core_parallel_bcs.ts_filter_half_width,
      .filter_cutoff_wavelength = app->core_parallel_bcs.ts_filter_cutoff_wavelength,
    };
    if (app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK) {
      tsinp_up.shift_dg = app->delta_ts_x_up;
    } else {
      tsinp_up.shift_func = app->gk_geom->upper_shift_func;
      tsinp_up.shift_func_ctx = app->gk_geom->upper_shift_ctx;
    }
    gks->bc_ts_up = gkyl_bc_twistshift_inew(&tsinp_up);
  }

  bool alloc_surr_aux_var = false;
  for (int d = 0; d < 2 * app->cdim; ++d) {
    alloc_surr_aux_var = alloc_surr_aux_var ||
                         gks->info.bcs[d].type == GKYL_BC_GK_SPECIES_SHEATH_SURROGATE;
  }

  if (alloc_surr_aux_var) {
    gk_species_moment_init(app, gks, &gks->sheath_moms, GKYL_F_MOMENT_MAXWELLIAN, false);
    gks->dens_sheath = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gks->temp_sheath = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  }

  // Set function pointers.
  gks->rhs_func = gk_species_kinetic_rhs_dynamic;
  gks->rhs_implicit_func = gk_species_kinetic_rhs_implicit_dynamic;
  gks->bc_func = gk_species_kinetic_apply_bc_dynamic;
  gks->release_func = gk_species_kinetic_release_dynamic;
  gks->step_f_func = gk_species_step_f_dynamic;
  gks->combine_func = gk_species_combine_dynamic;
  gks->copy_func = gk_species_copy_range_dynamic;
  gks->write_func = gk_species_write_dynamic;
  if (gks->info.write_omega_cfl) {
    gks->cflrate_ho = mkarr(false, gks->cflrate->ncomp, gks->cflrate->size);
    gks->write_cfl_func = gk_species_write_cfl_enabled;
  } else {
    gks->write_cfl_func = gk_species_write_cfl_disabled;
  }
  gks->write_mom_func = gk_species_write_mom_dynamic;
  gks->calc_integrated_mom_func = gk_species_calc_integrated_mom_dynamic;
  gks->write_integrated_mom_func = gk_species_write_integrated_mom_dynamic;
  gks->calc_L2norm_func = gk_species_calc_L2norm_dynamic;
  gks->write_L2norm_func = gk_species_write_L2norm_dynamic;
  if (gks->info.time_rate_diagnostics) {
    gks->calc_int_mom_dt_func = gk_species_calc_int_mom_dt_enabled;
  } else {
    gks->calc_int_mom_dt_func = gk_species_calc_int_mom_dt_disabled;
  }
}

void
gk_species_kinetic_read_distf_from_file(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gkyl_gyrokinetic_ic_import inp,
  struct gkyl_array *fout
)
{
  // Read a distribution function from a file onto fout. Intended options include importing:
  //   1) ICs with same grid.
  //   2) ICs one dimensionality lower (e.g. 2x2v for 3x2v sim).
  //   3) ICs with same grid extents but different resolution (NYI).

  struct gkyl_rect_grid grid = gks->grid;
  int pdim = grid.ndim;
  int cdim = app->cdim;
  int vdim = pdim - cdim;
  int poly_order = app->poly_order;

  struct gkyl_rect_grid grid_do; // Donor grid.
  struct gkyl_array_header_info hdr;
  int pdim_do, vdim_do, cdim_do;
  bool same_res = true;

  // Read the header of the input file, extract needed info an create a grid
  // and other things needed.
  FILE *fp;
  with_file(fp, inp.file_name, "r")
  {
    int status = gkyl_grid_sub_array_header_read_fp(&grid_do, &hdr, fp);

    pdim_do = grid_do.ndim;
    vdim_do = vdim; // Assume velocity space dimensionality is the same.
    cdim_do = pdim_do - vdim_do;

    // Perform some basic checks.
    if (pdim_do == pdim) {
      // Check if the grid resolution is the same.
      for (int d = 0; d < pdim; d++) {
        same_res = same_res && (grid_do.cells[d] == grid.cells[d]);
      }
    } else {
      // Assume the loaded file has one lower conf-space dimension.
      // Primarily meant for loading:
      //   - 1x2v for a 2x2v sim.
      //   - 2x2v for a 3x2v sim.
      assert(pdim_do == pdim - 1);
      double tol_eq = 1e-10;
      for (int d = 0; d < cdim_do - 1; d++) {
        assert(gkyl_compare_double(grid_do.lower[d], grid.lower[d], tol_eq));
        assert(gkyl_compare_double(grid_do.upper[d], grid.upper[d], tol_eq));
        assert(gkyl_compare_double(grid_do.dx[d], grid.dx[d], tol_eq));
        assert(grid_do.cells[d] == grid.cells[d]);
      }
      assert(gkyl_compare_double(grid_do.lower[cdim_do - 1], grid.lower[cdim - 1], tol_eq));
      assert(gkyl_compare_double(grid_do.upper[cdim_do - 1], grid.upper[cdim - 1], tol_eq));
      assert(gkyl_compare_double(grid_do.dx[cdim_do - 1], grid.dx[cdim - 1], tol_eq));
      assert(grid_do.cells[cdim_do - 1] == grid.cells[cdim - 1]);
      for (int d = 0; d < vdim; d++) {
        assert(gkyl_compare_double(grid_do.lower[cdim_do + d], grid.lower[cdim + d], tol_eq));
        assert(gkyl_compare_double(grid_do.upper[cdim_do + d], grid.upper[cdim + d], tol_eq));
        assert(gkyl_compare_double(grid_do.dx[cdim_do + d], grid.dx[cdim + d], tol_eq));
        assert(grid_do.cells[cdim_do + d] == grid.cells[cdim + d]);
      }
    }

    // Read basis info from header, check its consistency.
    struct gkyl_msgpack_map_elem elem_list[] = {
      {.key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0},
      {.key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = 0}
    };
    int elem_list_len = sizeof(elem_list) / sizeof(elem_list[0]);
    gkyl_msgpack_to_map_elem_list(
      &(struct gkyl_msgpack_data){.meta = hdr.meta, .meta_sz = hdr.meta_size}, elem_list_len,
      elem_list
    );
    assert(
      strcmp(
        gks->basis.id, gkyl_msgpack_map_elem_get_string(elem_list_len, elem_list, "basis_type")
      ) == 0
    );
    assert(poly_order == gkyl_msgpack_map_elem_get_uint(elem_list_len, elem_list, "poly_order"));
    gkyl_msgpack_map_elem_release_string(elem_list_len, elem_list, "basis_type");

    gkyl_grid_sub_array_header_release(&hdr);
  }

  // Donor basis.
  struct gkyl_basis basis_do;
  switch (gks->basis.b_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (poly_order > 1) {
        gkyl_cart_modal_serendip(&basis_do, pdim_do, poly_order);
      } else if (poly_order == 1) {
        // p=2 in vparallel
        gkyl_cart_modal_gkhybrid(&basis_do, cdim_do, vdim_do);
      }
      break;
    default:
      assert(false);
      break;
  }

  // Donor global range.
  int ghost_do[pdim_do];
  for (int d = 0; d < cdim_do; d++) {
    ghost_do[d] = 1;
  }
  for (int d = 0; d < vdim_do; d++) {
    ghost_do[cdim_do + d] = 0;
  }
  struct gkyl_range global_ext_do, global_do;
  gkyl_create_grid_ranges(&grid_do, ghost_do, &global_ext_do, &global_do);

  // Create a donor communicator.
  int cuts_tar[GKYL_MAX_DIM] = {-1}, cuts_do[GKYL_MAX_CDIM] = {-1};
  gkyl_rect_decomp_get_cuts(app->decomp, cuts_tar);
  if (cdim_do == cdim - 1) {
    for (int d = 0; d < cdim_do - 1; d++) {
      cuts_do[d] = cuts_tar[d];
    }
    cuts_do[cdim_do - 1] = cuts_tar[cdim - 1];
  } else {
    for (int d = 0; d < cdim; d++) {
      cuts_do[d] = cuts_tar[d];
    }
  }
  // Set velocity space cuts to 1 as we do not use MPI in vel-space.
  for (int d = 0; d < vdim; d++) {
    cuts_do[cdim_do + d] = 1;
  }

  struct gkyl_rect_decomp *decomp_do = gkyl_rect_decomp_new_from_cuts(pdim_do, cuts_do, &global_do);
  struct gkyl_comm *comm_do = gkyl_comm_split_comm(gks->comm, 0, decomp_do);

  // Donor local range.
  int my_rank = 0;
  gkyl_comm_get_rank(comm_do, &my_rank);

  struct gkyl_range local_ext_do, local_do;
  gkyl_create_ranges(&decomp_do->ranges[my_rank], ghost_do, &local_ext_do, &local_do);

  // Donor array.
  struct gkyl_array *fdo = mkarr(app->use_gpu, basis_do.num_basis, local_ext_do.volume);
  struct gkyl_array *fdo_host =
    app->use_gpu ? mkarr(false, basis_do.num_basis, local_ext_do.volume) : gkyl_array_acquire(fdo);

  // Read donor distribution function and Jacobian inverse.
  struct gkyl_app_restart_status rstat;
  rstat.io_status = gkyl_comm_array_read(comm_do, &grid_do, &local_do, fdo_host, inp.file_name);

  bool scale_by_jacobtot = false;
  with_file(fp, inp.jacobtot_inv_file_name, "r")
  {
    // Configuration space donor grid and basis
    struct gkyl_rect_grid conf_grid_do;
    gkyl_rect_grid_init(&conf_grid_do, cdim_do, grid_do.lower, grid_do.upper, grid_do.cells);
    struct gkyl_basis conf_basis_do;
    gkyl_cart_modal_serendip(&conf_basis_do, cdim_do, poly_order);
    // Configuration space donor ranges and decomposition.
    struct gkyl_range conf_local_ext_do, conf_local_do, conf_global_ext_do, conf_global_do;
    gkyl_create_grid_ranges(&conf_grid_do, ghost_do, &conf_global_ext_do, &conf_global_do);
    struct gkyl_rect_decomp *conf_decomp_do =
      gkyl_rect_decomp_new_from_cuts(cdim_do, cuts_do, &conf_global_do);
    gkyl_create_ranges(
      &conf_decomp_do->ranges[my_rank], ghost_do, &conf_local_ext_do, &conf_local_do
    );
    // Read and multiply by the reciprocal of jacobtot.
    struct gkyl_array *jacobtot_inv_do_host =
      mkarr(false, conf_basis_do.num_basis, conf_local_ext_do.volume);
    rstat.io_status = gkyl_comm_array_read(
      comm_do, &conf_grid_do, &conf_local_do, jacobtot_inv_do_host, inp.jacobtot_inv_file_name
    );
    gkyl_dg_mul_conf_phase_op_range(
      &conf_basis_do, &basis_do, fdo_host, jacobtot_inv_do_host, fdo_host, &conf_local_ext_do,
      &local_ext_do
    );
    gkyl_array_release(jacobtot_inv_do_host);
    gkyl_rect_decomp_release(conf_decomp_do);
    scale_by_jacobtot = true;
  }

  bool scale_by_jacobvel = false;
  with_file(fp, inp.jacobvel_file_name, "r")
  {
    // Read and multiply by the reciprocal of jacobvel.
    struct gkyl_array *jacobvel_do_host = mkarr(false, 1, local_ext_do.volume);
    rstat.io_status =
      gkyl_comm_array_read(comm_do, &grid_do, &local_do, jacobvel_do_host, inp.jacobvel_file_name);
    gkyl_array_divide_by_cell(fdo_host, jacobvel_do_host);
    gkyl_array_release(jacobvel_do_host);
    scale_by_jacobvel = true;
  }

  if (app->use_gpu) {
    gkyl_array_copy(fdo, fdo_host);
  }

  if (pdim_do == pdim - 1) {
    struct gkyl_translate_dim *transdim =
      gkyl_translate_dim_new(cdim_do, basis_do, cdim, gks->basis, -1, GKYL_NO_EDGE, app->use_gpu);
    gkyl_translate_dim_advance(transdim, &local_do, &gks->local, fdo, 1, fout);
    gkyl_translate_dim_release(transdim);
  } else {
    if (same_res) {
      gkyl_array_copy(fout, fdo);
    } else {
      // Interpolate the donor distribution to the target grid.
      struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(
        app->cdim, &gks->basis, &grid_do, &grid, &local_do, &gks->local, ghost_do, app->use_gpu
      );
      gkyl_dg_interpolate_advance(interp, fdo, fout);
      gkyl_dg_interpolate_release(interp);
    }
  }

  if (inp.type == GKYL_IC_IMPORT_AF) {
    // Scale fout by a conf-space factor.
    gkyl_proj_on_basis *proj_conf_scale = gkyl_proj_on_basis_new(
      &app->grid, &app->basis, poly_order + 1, 1, inp.conf_scale, inp.conf_scale_ctx
    );
    struct gkyl_array *xfac = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *xfac_ho = app->use_gpu ?
                                   mkarr(false, app->basis.num_basis, app->local_ext.volume) :
                                   gkyl_array_acquire(xfac);
    gkyl_proj_on_basis_advance(proj_conf_scale, 0.0, &app->local, xfac_ho);
    gkyl_array_copy(xfac, xfac_ho);
    gkyl_dg_mul_conf_phase_op_range(
      &app->basis, &gks->basis, fout, xfac, fout, &app->local, &gks->local
    );
    gkyl_proj_on_basis_release(proj_conf_scale);
    gkyl_array_release(xfac_ho);
    gkyl_array_release(xfac);
  }

  // Multiply fout by the Jacobian.
  if (scale_by_jacobtot) {
    gkyl_dg_mul_conf_phase_op_range(
      &app->basis, &gks->basis, fout, app->gk_geom->geo_int.jacobtot, fout, &app->local, &gks->local
    );
  }

  // Multiply fout by the velocity space Jacobian.
  if (scale_by_jacobvel) {
    gkyl_array_scale_by_cell(fout, gks->vel_map->jacobvel);
  }

  gkyl_rect_decomp_release(decomp_do);
  gkyl_comm_release(comm_do);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_host);
}

void
gk_species_kinetic_file_import_init(
  struct gkyl_gyrokinetic_app *app, struct gk_species *gks, struct gkyl_gyrokinetic_ic_import inp
)
{
  gk_species_kinetic_read_distf_from_file(app, gks, inp, gks->f);

  if (inp.enforce_positivity) {
    // Positivity enforcing by shifting f (ps=positivity shift).
    struct gkyl_positivity_shift_gyrokinetic *pos_shift_op = gkyl_positivity_shift_gyrokinetic_new(
      app->basis, gks->basis, gks->grid, gks->info.mass, app->gk_geom, gks->vel_map,
      &app->local_ext, app->use_gpu
    );

    gkyl_positivity_shift_gyrokinetic_advance(
      pos_shift_op, &app->local, &gks->local, gks->f, gks->m0.marr, gks->m0.marr
    );

    gkyl_positivity_shift_gyrokinetic_release(pos_shift_op);
  }
}

static bool
gk_species_kinetic_do_I_recycle(struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Check whether one of the neutral species has recycling BCs that depend on
  // this gyrokinetic species.
  bool recycling_bcs = false;
  int neuts = app->num_neut_species;
  for (int i = 0; i < neuts; ++i) {
    for (int k = 0; k < 2 * app->cdim; k++) {
      const struct gkyl_gyrokinetic_bc *bc = &app->neut_species[i].info.bcs[k];
      if (bc->type == GKYL_BC_GK_SPECIES_RECYCLE) {
        for (int j = 0; j < bc->emission.num_species; j++) {
          recycling_bcs = recycling_bcs || 0 == strcmp(gks->info.name, bc->emission.in_species[j]);
        }
      }
    }
  }
  return recycling_bcs;
}

static bool
gk_species_kinetic_do_I_recycle_react_scale(struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Check whether one of the neutral species has a recycle_react_scale
  // operation thats depend on this gyrokinetic species.
  bool has_sca = false;
  int neuts = app->num_neut_species;
  for (int i = 0; i < neuts; ++i) {
    struct gk_neut_species *gkns = &app->neut_species[i];
    struct gkyl_gyrokinetic_scaling_inp *sca_inp = &gkns->info.scaling;
    if ((sca_inp->num_boundaries > 0) &&
        (0 == strcmp(gks->info.name, sca_inp->impacting_ion_name))) {
      has_sca = true;
      break;
    }
  }
  return has_sca;
}

static bool
gk_species_kinetic_do_I_adapt_src(struct gkyl_gyrokinetic_app *app, struct gk_species *gks)
{
  // Check whether one of the species adapts its source depending on
  // this gyrokinetic species.
  bool adapt_src = false;
  for (int i = 0; i < app->num_species; ++i) {
    struct gkyl_gyrokinetic_source *source = &app->species[i].info.source;
    if (source->num_adapt_sources > 0) {
      for (int j = 0; j < source->num_adapt_sources; j++) {
        adapt_src = adapt_src || 0 == strcmp(gks->info.name, source->adapt[j].adapt_to_species);
      }
    }
  }
  return adapt_src;
}

void
gk_species_kinetic_init(
  struct gkyl_gk *gk_app_inp, struct gkyl_gyrokinetic_app *app, struct gk_species *gks
)
{
  gks->release_species_func = gk_species_kinetic_release;
  int cdim = app->cdim, vdim = gks->info.vdim;
  int pdim = cdim + vdim;

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  assert(vdim > 0); // Ensure user provided vdim in input file.
  assert(gks->info.charge != 0.0); // Charged GK operators require a nonzero charge.

  for (int d = 0; d < cdim; ++d) {
    cells[d] = gk_app_inp->cells[d];
    lower[d] = gk_app_inp->lower[d];
    upper[d] = gk_app_inp->upper[d];
    ghost[d] = 1;
  }
  for (int d = 0; d < vdim; ++d) {
    cells[cdim + d] = gks->info.cells[d];
    lower[cdim + d] = gks->info.lower[d];
    upper[cdim + d] = gks->info.upper[d];
    ghost[cdim + d] = 0; // No ghost-cells in velocity space.

    cells_vel[d] = gks->info.cells[d];
    lower_vel[d] = gks->info.lower[d];
    upper_vel[d] = gks->info.upper[d];
    ghost_vel[d] = 0; // No ghost-cells in velocity space.
  }

  // Define phase basis
  if (app->use_gpu) {
    // allocate device basis if we are using GPUs
    gks->basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  } else {
    gks->basis_on_dev = &gks->basis;
  }

  enum gkyl_basis_type b_type = app->basis.b_type;

  // basis functions
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (app->poly_order > 1) {
        gkyl_cart_modal_serendip(&gks->basis, pdim, app->poly_order);
      } else if (app->poly_order == 1) {
        gkyl_cart_modal_gkhybrid(&gks->basis, cdim, vdim); // p=2 in vparallel
      }

      if (app->use_gpu) {
        if (app->poly_order > 1) {
          gkyl_cart_modal_serendip_cu_dev(gks->basis_on_dev, pdim, app->poly_order);
        } else if (app->poly_order == 1) {
          gkyl_cart_modal_gkhybrid_cu_dev(gks->basis_on_dev, cdim, vdim); // p=2 in vparallel
        }
      }
      break;
    default:
      assert(false);
      break;
  }

  // Full phase space grid.
  gkyl_rect_grid_init(&gks->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&gks->grid, ghost, &gks->global_ext, &gks->global);

  // Velocity space grid.
  gkyl_rect_grid_init(&gks->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&gks->grid_vel, ghost_vel, &gks->local_ext_vel, &gks->local_vel);

  // Phase-space communicator.
  gks->comm = gkyl_comm_extend_comm(app->comm, &gks->local_vel);

  // Create local and local_ext from app local range.
  struct gkyl_range local;
  // Local = conf-local X local_vel.
  gkyl_range_ten_prod(&local, &app->local, &gks->local_vel);
  gkyl_create_ranges(&local, ghost, &gks->local_ext, &gks->local);

  // Velocity space mapping.
  gks->vel_map = gkyl_velocity_map_new(
    gks->info.mapc2p, gks->grid, gks->grid_vel, gks->local, gks->local_ext, gks->local_vel,
    gks->local_ext_vel, app->use_gpu
  );

  // Write out the velocity space mapping and its Jacobian.
  gkyl_velocity_map_write(gks->vel_map, gks->comm, app->name, gks->info.name);

  // Keep a copy of num_periodic_dir and periodic_dirs in species so we can
  // add the parallel direction in case TS BCs are needed.
  gks->num_periodic_dir = app->num_periodic_dir;
  for (int d = 0; d < gks->num_periodic_dir; ++d) {
    gks->periodic_dirs[d] = app->periodic_dirs[d];
  }

  for (int d = 0; d < app->cdim; ++d) {
    gks->bc_is_np[d] = true;
  }
  for (int d = 0; d < gks->num_periodic_dir; ++d) {
    gks->bc_is_np[gks->periodic_dirs[d]] = false;
  }

  // Store the BCs from the input file.
  for (int d = 0; d < app->cdim; ++d) {
    if (gks->bc_is_np[d]) {
      struct gkyl_gyrokinetic_bc *bc_lo =
        gk_fetch_bc_with_dir_edge(gks->info.bcs, 2 * app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        gks->lower_bc[d] = *bc_lo;
      } else {
        gks->lower_bc[d].type = GKYL_BC_GK_SKIP;
      }

      struct gkyl_gyrokinetic_bc *bc_up =
        gk_fetch_bc_with_dir_edge(gks->info.bcs, 2 * app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        gks->upper_bc[d] = *bc_up;
      } else {
        gks->upper_bc[d].type = GKYL_BC_GK_SKIP;
      }
    } else {
      gks->lower_bc[d].type = GKYL_BC_GK_SPECIES_PERIODIC;
      gks->upper_bc[d].type = GKYL_BC_GK_SPECIES_PERIODIC;
    }
  }

  gk_species_init_io(app, gks);

  // Allocate distribution function arrays.
  gks->f = mkarr(app->use_gpu, gks->basis.num_basis, gks->local_ext.volume);

  gks->f_host = gks->f;
  if (app->use_gpu) {
    gks->f_host = mkarr(false, gks->basis.num_basis, gks->local_ext.volume);
  }

  // Allocate cflrate (scalar array).
  gks->cflrate = mkarr(app->use_gpu, 1, gks->local_ext.volume);

  // Allocate data for density (for charge density or upar calculation).
  gk_species_moment_init(app, gks, &gks->m0, GKYL_F_MOMENT_M0, false);

  if (gks->info.flr.type) {
    // Create operator needed for FLR effects.
    assert(app->cdim > 1);
    // Pointer to function performing the gyroaverage.
    gks->gyroaverage = gk_species_gyroaverage_enabled;
    // Gyroaveraged M0 and phi.
    gks->m0_gyroavg = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gks->gyro_phi = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);

    double gyroradius_bmag = gks->info.flr.bmag ? gks->info.flr.bmag : app->bmag_ref;

    double flr_weight =
      gks->info.flr.Tperp * gks->info.mass / (pow(gks->info.charge * gyroradius_bmag, 2.0));
    // Initialize the weight in the Laplacian operator.
    gks->flr_rhoSqD2 =
      mkarr(app->use_gpu, (2 * (app->cdim - 1) - 1) * app->basis.num_basis, app->local_ext.volume);
    gkyl_array_set_offset(
      gks->flr_rhoSqD2, flr_weight, app->gk_geom->geo_int.gxxj, 0 * app->basis.num_basis
    );
    if (app->cdim > 2) {
      gkyl_array_set_offset(
        gks->flr_rhoSqD2, flr_weight, app->gk_geom->geo_int.gxyj, 1 * app->basis.num_basis
      );
      gkyl_array_set_offset(
        gks->flr_rhoSqD2, flr_weight, app->gk_geom->geo_int.gyyj, 2 * app->basis.num_basis
      );
    }
    // Initialize the factor multiplying the field in the FLR operator.
    gks->flr_kSq = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_array_shiftc(gks->flr_kSq, -pow(sqrt(2.0), app->cdim), 0); // Sets kSq=-1.

    // Periodic along the periodic directions and Dirichlet otherwise, unless the field BCs say so.
    struct gkyl_poisson_bc flr_bc = {};
    for (int d = 0; d < app->cdim - 1; d++) {
      enum gkyl_poisson_bc_type bc_default =
        gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
      for (int i = 0; i < app->num_periodic_dir; i++) {
        if (app->periodic_dirs[i] == d) {
          bc_default = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
        }
      }
      flr_bc.lo_type[d] = bc_default;
      flr_bc.up_type[d] = bc_default;

      struct gkyl_gyrokinetic_bc *bc_lo =
        gk_fetch_bc_with_dir_edge(app->field->info.poisson_bcs, 2 * app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        if (bc_lo->type == GKYL_BC_GK_FIELD_PERIODIC) {
          flr_bc.lo_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
        } else {
          flr_bc.lo_type[d] =
            gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
        }
      }

      struct gkyl_gyrokinetic_bc *bc_up =
        gk_fetch_bc_with_dir_edge(app->field->info.poisson_bcs, 2 * app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        if (bc_up->type == GKYL_BC_GK_FIELD_PERIODIC) {
          flr_bc.up_type[d] = gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_PERIODIC);
        } else {
          flr_bc.up_type[d] =
            gkyl_gyrokinetic_translate_poisson_bc_type(GKYL_BC_GK_FIELD_DIRICHLET_VARYING);
        }
      }
    }
    // Deflated Poisson solve is performed on range assuming decomposition is *only* in z.
    gks->flr_op = gkyl_deflated_fem_poisson_new(
      app->grid, app->basis_on_dev, app->basis, app->local, app->local, gks->flr_rhoSqD2,
      gks->flr_kSq, flr_bc, NULL, app->use_gpu
    );
  } else {
    gks->gyroaverage = gk_species_gyroaverage_disabled;
    gks->m0_gyroavg = gkyl_array_acquire(gks->m0.marr);
    gks->gyro_phi = gkyl_array_acquire(app->field->phi_smooth);
  }

  // Initialize the collisionless solver.
  gks->collisionless = (struct gk_collisionless){};
  gk_species_collisionless_init(app, gks, &gks->collisionless);

  // Initialize an anomalous diffusion term.
  gks->anom_diff = (struct gk_anomalous_diff){};
  gk_species_anomalous_diff_init(app, gks, &gks->anom_diff);

  // Damping term -nu*f on RHS.
  gk_species_damping_init(app, gks, &gks->damping);

  // Function multiplying df/dt.
  gk_species_fdot_multiplier_init(app, gks, &gks->fdot_mult);

  // Allocate data for diagnostic moments.
  int ndm = gks->info.num_diag_moments;
  gks->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m = 0; m < ndm; ++m) {
    gk_species_moment_init(app, gks, &gks->moms[m], gks->info.diag_moments[m], false);
  }

  // Initialize projection routine for initial conditions.
  if (gks->info.init_from_file.type == 0) {
    gk_species_projection_init(app, gks, gks->info.projection, &gks->proj_init);
  } else {
    gk_species_file_import_init(app, gks, gks->info.init_from_file);
  }

  // Create skin/ghost ranges.
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &gks->local_lower_skin[dir], &gks->local_lower_ghost[dir], dir, GKYL_LOWER_EDGE,
      &gks->local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &gks->local_upper_skin[dir], &gks->local_upper_ghost[dir], dir, GKYL_UPPER_EDGE,
      &gks->local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &gks->global_lower_skin[dir], &gks->global_lower_ghost[dir], dir, GKYL_LOWER_EDGE,
      &gks->global_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &gks->global_upper_skin[dir], &gks->global_upper_ghost[dir], dir, GKYL_UPPER_EDGE,
      &gks->global_ext, ghost
    );
  }

  // Create a local range extended in the BC dir (for TS BCs).
  int par_dir = app->cdim - 1;
  int ndim = gks->local.ndim;
  int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
  for (int i = 0; i < ndim; i++) {
    lower_bcdir_ext[i] = gks->local.lower[i];
    upper_bcdir_ext[i] = gks->local.upper[i];
  }
  lower_bcdir_ext[par_dir] = gks->local_ext.lower[par_dir];
  upper_bcdir_ext[par_dir] = gks->local_ext.upper[par_dir];
  gkyl_sub_range_init(&gks->local_par_ext, &gks->local_ext, lower_bcdir_ext, upper_bcdir_ext);

  if (gk_app_inp->geometry.has_LCFS) {
    // Simulation spans the last-closed flux surface (LCFS). Create core and SOL global ranges.
    int idx_LCFS_lo = app->gk_geom->idx_LCFS_lo;
    // Length of lower and upper x ranges (one is core, the other SOL).
    int len_lo = idx_LCFS_lo;
    int len_up = gks->global.upper[0] - len_lo;
    // Lower and upper x ranges.
    struct gkyl_range *global_lo_r, *global_ext_lo_r, *global_up_r, *global_ext_up_r;
    struct gkyl_range *local_lo_r, *local_ext_lo_r, *local_up_r, *local_ext_up_r;
    struct gkyl_range *local_lower_skin_par_lo_r, *local_upper_skin_par_lo_r,
      *local_lower_ghost_par_lo_r, *local_upper_ghost_par_lo_r;
    struct gkyl_range *local_lower_skin_par_up_r, *local_upper_skin_par_up_r,
      *local_lower_ghost_par_up_r, *local_upper_ghost_par_up_r;
    if (app->gk_geom->geqdsk_sign_convention == 0) {
      // x increases towards SOL.
      global_lo_r = &gks->global_core;
      global_up_r = &gks->global_sol;
      global_ext_lo_r = &gks->global_ext_core;
      global_ext_up_r = &gks->global_ext_sol;
      local_lo_r = &gks->local_core;
      local_up_r = &gks->local_sol;
      local_ext_lo_r = &gks->local_ext_core;
      local_ext_up_r = &gks->local_ext_sol;
      local_lower_skin_par_lo_r = &gks->local_lower_skin_par_core;
      local_lower_ghost_par_lo_r = &gks->local_lower_ghost_par_core;
      local_lower_skin_par_up_r = &gks->local_lower_skin_par_sol;
      local_lower_ghost_par_up_r = &gks->local_lower_ghost_par_sol;
      local_upper_skin_par_lo_r = &gks->local_upper_skin_par_core;
      local_upper_ghost_par_lo_r = &gks->local_upper_ghost_par_core;
      local_upper_skin_par_up_r = &gks->local_upper_skin_par_sol;
      local_upper_ghost_par_up_r = &gks->local_upper_ghost_par_sol;
    } else {
      // x increases towards core.
      global_lo_r = &gks->global_sol;
      global_ext_lo_r = &gks->global_ext_sol;
      global_up_r = &gks->global_core;
      global_ext_up_r = &gks->global_ext_core;
      local_lo_r = &gks->local_sol;
      local_up_r = &gks->local_core;
      local_ext_lo_r = &gks->local_ext_sol;
      local_ext_up_r = &gks->local_ext_core;
      local_lower_skin_par_lo_r = &gks->local_lower_skin_par_sol;
      local_upper_skin_par_lo_r = &gks->local_upper_skin_par_sol;
      local_lower_ghost_par_lo_r = &gks->local_lower_ghost_par_sol;
      local_upper_ghost_par_lo_r = &gks->local_upper_ghost_par_sol;
      local_lower_skin_par_up_r = &gks->local_lower_skin_par_core;
      local_upper_skin_par_up_r = &gks->local_upper_skin_par_core;
      local_lower_ghost_par_up_r = &gks->local_lower_ghost_par_core;
      local_upper_ghost_par_up_r = &gks->local_upper_ghost_par_core;
    }

    // Lower and upper x ranges.
    gkyl_range_shorten_from_above(global_lo_r, &gks->global, 0, len_lo);
    gkyl_range_shorten_from_below(global_up_r, &gks->global, 0, len_up);
    gkyl_range_shorten_from_above(local_lo_r, &gks->local, 0, len_lo);
    gkyl_range_shorten_from_below(local_up_r, &gks->local, 0, len_up);

    // Extended lower and upper x ranges.
    int len_lo_ext = idx_LCFS_lo + 1;
    int len_up_ext = gks->global_ext.upper[0] - len_lo;
    gkyl_range_shorten_from_above(global_ext_lo_r, &gks->global_ext, 0, len_lo_ext);
    gkyl_range_shorten_from_below(global_ext_up_r, &gks->global_ext, 0, len_up_ext);
    gkyl_range_shorten_from_above(local_ext_lo_r, &gks->local_ext, 0, len_lo_ext);
    gkyl_range_shorten_from_below(local_ext_up_r, &gks->local_ext, 0, len_up_ext);

    // Parallel skin and ghost ranges, limited to the lower and upper x range.
    int par_dir = app->cdim - 1;
    for (int e = 0; e < 2; e++) {
      gkyl_range_shorten_from_above(
        e == 0 ? local_lower_skin_par_lo_r : local_upper_skin_par_lo_r,
        e == 0 ? &gks->local_lower_skin[par_dir] : &gks->local_upper_skin[par_dir], 0, len_lo
      );
      gkyl_range_shorten_from_above(
        e == 0 ? local_lower_ghost_par_lo_r : local_upper_ghost_par_lo_r,
        e == 0 ? &gks->local_lower_ghost[par_dir] : &gks->local_upper_ghost[par_dir], 0, len_lo
      );
      gkyl_range_shorten_from_below(
        e == 0 ? local_lower_skin_par_up_r : local_upper_skin_par_up_r,
        e == 0 ? &gks->local_lower_skin[par_dir] : &gks->local_upper_skin[par_dir], 0, len_up
      );
      gkyl_range_shorten_from_below(
        e == 0 ? local_lower_ghost_par_up_r : local_upper_ghost_par_up_r,
        e == 0 ? &gks->local_lower_ghost[par_dir] : &gks->local_upper_ghost[par_dir], 0, len_up
      );
    }

    // Core range extended in the parallel direction.
    int ndim = gks->local.ndim;
    int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
    for (int i = 0; i < ndim; i++) {
      lower_bcdir_ext[i] = gks->local_core.lower[i];
      upper_bcdir_ext[i] = gks->local_core.upper[i];
    }
    lower_bcdir_ext[par_dir] = gks->local_ext_core.lower[par_dir];
    upper_bcdir_ext[par_dir] = gks->local_ext_core.upper[par_dir];
    gkyl_sub_range_init(
      &gks->local_par_ext_core, &gks->local_ext_core, lower_bcdir_ext, upper_bcdir_ext
    );
  }

  // Initialize boundary fluxes.
  gks->bflux = (struct gk_boundary_fluxes){};
  // Additional bflux moments to step in time.
  struct gkyl_phase_diagnostics_inp add_bflux_moms_inp = (struct gkyl_phase_diagnostics_inp){};
  enum gkyl_species_bflux_type bflux_type = GK_SPECIES_BFLUX_NONE;
  // Check if using Boltzmann elc.
  bool ion_in_boltz_elc_field = (gks->info.charge > 0.0) && app->field->update_field &&
                                (app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN);
  // Check if sources are adaptive.
  bool adaptive_sources = gk_species_kinetic_do_I_adapt_src(app, gks);
  // Check if other species use the recycle_react_scale operation.
  bool recycle_react_scale = gk_species_kinetic_do_I_recycle_react_scale(app, gks);
  // Check if other species have recycling BCs.
  bool recycling_bcs = gk_species_kinetic_do_I_recycle(app, gks);
  if (gks->info.boundary_flux_diagnostics.num_diag_moments > 0 ||
      gks->info.boundary_flux_diagnostics.num_integrated_diag_moments > 0) {
    bflux_type = GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS_DIAGS;
  } else {
    // Set bflux_type to
    //   - GK_SPECIES_BFLUX_CALC_FLUX to only put bfluxes in ghost cells of rhs.
    //   - GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS to calc bfluxes and step its moments.
    // The latter also requires that you place the moment you desire in add_bflux_moms_inp below.

    if (recycling_bcs) {
      bflux_type = GK_SPECIES_BFLUX_CALC_FLUX;
    }
    if (ion_in_boltz_elc_field || adaptive_sources || recycle_react_scale) {
      // This is within an if-statement instead of an else if because it
      // superseeds (and is a superset) of GK_SPECIES_BFLUX_CALC_FLUX.
      bflux_type = GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS;
    }
  }
  int *nmom_extra = &add_bflux_moms_inp.num_diag_moments;
  if (ion_in_boltz_elc_field || adaptive_sources || recycle_react_scale) {
    add_bflux_moms_inp.diag_moments[nmom_extra[0]++] = GKYL_F_MOMENT_M0;
  }
  if (adaptive_sources) {
    add_bflux_moms_inp.diag_moments[nmom_extra[0]++] = GKYL_F_MOMENT_M2;
  }
  // Introduce new moments into moms_inp if needed.
  gk_species_bflux_init(app, gks, &gks->bflux, bflux_type, add_bflux_moms_inp);

  // Initialize empty structs. New methods will fill them if specified.
  gks->src = (struct gk_source){};
  gks->rad = (struct gk_rad_drag){};

  // Initialize a Maxwellian/LTE (local thermodynamic equilibrium) projection routine
  // Projection routine optionally corrects all the Maxwellian/LTE moments
  // This routine is utilized by both reactions and BGK collisions
  gks->lte = (struct gk_lte){};
  bool correct_all_moms = gks->info.correct.correct_all_moms;
  int max_iter = gks->info.correct.max_iter > 0 ? gks->info.correct.max_iter : 50;
  double iter_eps = gks->info.correct.iter_eps > 0 ? gks->info.correct.iter_eps : 1e-10;
  bool use_last_converged = gks->info.correct.use_last_converged;
  struct correct_all_moms_inp corr_inp = {
    .correct_all_moms = correct_all_moms,
    .max_iter = max_iter,
    .iter_eps = iter_eps,
    .use_last_converged = use_last_converged,
  };
  gk_species_lte_init(app, gks, &gks->lte, corr_inp);

  // Initialize the object that scales the species (meant for fluid neutrals for now).
  gks->sca = (struct gk_scaling){};
  gk_species_scaling_init(app, gks, &gks->sca);

  // Initialize reactions with charged species.
  gks->react = (struct gk_react){};
  gk_species_react_init(app, gks, gks->info.react, &gks->react, true);

  // Initialize reactions with neutral species.
  gks->react_neut = (struct gk_react){};
  gk_species_react_init(app, gks, gks->info.react_neut, &gks->react_neut, false);

  // Initialize LBO collisions.
  gks->lbo = (struct gk_lbo_collisions){};
  gk_species_lbo_init(app, gks, &gks->lbo);

  // Initialize BGK collisions.
  gks->bgk = (struct gk_bgk_collisions){};
  gk_species_bgk_init(app, gks, &gks->bgk);

  // Initialize a BGK source.
  gks->bgk_src = (struct gk_source_bgk){};
  gk_species_source_bgk_init(app, gks, &gks->bgk_src);

  // Initialize positivity enforcing operator.
  gks->positivity = (struct gk_positivity){};
  gk_species_positivity_init(app, gks, &gks->positivity);

  if (gks->info.is_static) {
    gk_species_init_static(gk_app_inp, app, gks);
  } else {
    gk_species_kinetic_init_dynamic(gk_app_inp, app, gks);
  }
}

void
gk_species_kinetic_release(const gkyl_gyrokinetic_app *app, const struct gk_species *gks)
{
  // Release resources for charged species.
  gkyl_array_release(gks->f);

  gkyl_array_release(gks->cflrate);

  if (gks->info.init_from_file.type == 0) {
    gk_species_projection_release(app, &gks->proj_init);
  }

  gkyl_comm_release(gks->comm);

  gkyl_msgpack_map_elem_release(gks->io_meta_basic_len, gks->io_meta_basic);
  gkyl_msgpack_map_elem_release(gks->io_meta_conf_len, gks->io_meta_conf);
  gkyl_msgpack_map_elem_release(gks->io_meta_phase_len, gks->io_meta_phase);

  if (app->use_gpu) {
    gkyl_array_release(gks->f_host);
    gkyl_cu_free(gks->basis_on_dev);
  }

  gkyl_velocity_map_release(gks->vel_map);

  gk_species_collisionless_release(app, &gks->collisionless);

  gk_species_scaling_release(app, &gks->sca);

  gk_species_anomalous_diff_release(app, &gks->anom_diff);

  // Release moment data.
  gk_species_moment_release(app, &gks->m0);
  for (int i = 0; i < gks->info.num_diag_moments; ++i) {
    gk_species_moment_release(app, &gks->moms[i]);
  }
  gkyl_free(gks->moms);

  gk_species_source_release(app, &gks->src);

  gk_species_damping_release(app, &gks->damping);

  gk_species_fdot_multiplier_release(app, &gks->fdot_mult);

  gk_species_lbo_release(app, &gks->lbo);

  gk_species_bgk_release(app, &gks->bgk);

  gk_species_source_bgk_release(app, &gks->bgk_src);

  gk_species_positivity_release(app, &gks->positivity);

  gk_species_radiation_release(app, &gks->rad);

  gk_species_react_release(app, &gks->react);

  gk_species_react_release(app, &gks->react_neut);

  // Free boundary flux memory.
  gk_species_bflux_release(app, gks, &gks->bflux);

  gk_species_lte_release(app, &gks->lte);

  gkyl_array_release(gks->m0_gyroavg);
  gkyl_array_release(gks->gyro_phi);
  if (gks->info.flr.type) {
    gkyl_array_release(gks->flr_rhoSqD2);
    gkyl_array_release(gks->flr_kSq);
    gkyl_deflated_fem_poisson_release(gks->flr_op);
  }

  gks->release_func(app, gks);
}
