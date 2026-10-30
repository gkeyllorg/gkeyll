#include <gkyl_gk_neut_species_priv.h>

static double
gk_neut_species_kinetic_rhs_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  double omega_cfl = 1 / DBL_MAX;
  gkyl_array_clear(gkns->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Collisionless terms.
  gk_neut_species_collisionless_rhs(app, gkns, &gkns->collisionless, fin, rhs);

  // Elastic collisions.
  gk_neut_species_bgk_rhs(app, gkns, &gkns->bgk, fin, rhs);

  // Reactions with charged species.
  gk_neut_species_react_rhs(app, gkns, &gkns->react_neut, fin, rhs);

  // Compute and store (in the ghost cell of rhs) the boundary fluxes.
  gk_neut_species_bflux_rhs(app, &gkns->bflux, fin, rhs);

  // Compute diagnostic moments of the boundary fluxes.
  gk_neut_species_bflux_calc_moms(app, &gkns->bflux, rhs, bflux_moms);

  app->stat.n_neut_species_omega_cfl += 1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(gkns->omega_cfl, gkns->cflrate, GKYL_MAX, &gkns->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, gkns->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    omega_cfl_ho[0] = gkns->omega_cfl[0];
  }
  omega_cfl = omega_cfl_ho[0];

  app->stat.neut_species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl / omega_cfl;
}

static double
gk_neut_species_kinetic_rhs_implicit_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt
)
{
  double omega_cfl = 1 / DBL_MAX;
  gkyl_array_clear(gkns->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Compute implicit update and update rhs to new time step.
  gk_neut_species_bgk_rhs_implicit(app, gkns, &gkns->bgk, fin, dt, rhs);

  gkyl_array_accumulate(gkyl_array_scale(rhs, dt), 1.0, fin);

  app->stat.n_neut_species_omega_cfl += 1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(gkns->omega_cfl, gkns->cflrate, GKYL_MAX, &gkns->local);

  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, gkns->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    omega_cfl_ho[0] = gkns->omega_cfl[0];
  }
  omega_cfl = omega_cfl_ho[0];

  app->stat.neut_species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl / omega_cfl;
}

static void
gk_neut_species_kinetic_apply_bc_dynamic(
  gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns, struct gkyl_array *f
)
{
  struct timespec wst = gkyl_wall_clock();

  int num_periodic_dir = app->num_periodic_dir, cdim = app->cdim;
  gkyl_comm_array_per_sync(
    gkns->comm, &gkns->local, &gkns->local_ext, num_periodic_dir, app->periodic_dirs, f
  );

  for (int d = 0; d < cdim; ++d) {
    if (gkns->bc_is_np[d]) {
      switch (gkns->lower_bc[d].type) {
        case GKYL_BC_GK_SPECIES_RECYCLE:
          gk_neut_species_recycle_apply_bc(app, &gkns->bc_recycle_lo, gkns, f);
          break;
        case GKYL_BC_GK_SPECIES_COPY:
        case GKYL_BC_GK_SPECIES_REFLECT:
        case GKYL_BC_GK_SPECIES_ABSORB:
          gkyl_bc_basic_advance(gkns->bc_lo[d], gkns->bc_buffer, f);
          break;
        case GKYL_BC_GK_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(gkns->bc_lo[d], gkns->bc_buffer_lo_fixed, f);
          break;
        case GKYL_BC_GK_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }

      switch (gkns->upper_bc[d].type) {
        case GKYL_BC_GK_SPECIES_RECYCLE:
          gk_neut_species_recycle_apply_bc(app, &gkns->bc_recycle_up, gkns, f);
          break;
        case GKYL_BC_GK_SPECIES_COPY:
        case GKYL_BC_GK_SPECIES_REFLECT:
        case GKYL_BC_GK_SPECIES_ABSORB:
          gkyl_bc_basic_advance(gkns->bc_up[d], gkns->bc_buffer, f);
          break;
        case GKYL_BC_GK_SPECIES_FIXED_FUNC:
          gkyl_bc_basic_advance(gkns->bc_up[d], gkns->bc_buffer_up_fixed, f);
          break;
        case GKYL_BC_GK_SPECIES_ZERO_FLUX:
          break; // do nothing, BCs already applied in hyper_dg loop by not updating flux
        default:
          break;
      }
    }
  }

  gkyl_comm_array_sync(gkns->comm, &gkns->local, &gkns->local_ext, f);

  app->stat.neut_species_bc_tm += gkyl_time_diff_now_sec(wst);
}

// release all resources for dynamic species
static void
gk_neut_species_kinetic_release_dynamic(
  const gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns
)
{
  // release various arrays
  gkyl_array_release(gkns->bc_buffer);
  gkyl_array_release(gkns->bc_buffer_lo_fixed);
  gkyl_array_release(gkns->bc_buffer_up_fixed);

  gkyl_array_release(gkns->f1);
  gkyl_array_release(gkns->fnew);
  gkyl_array_release(gkns->cflrate);

  if (app->use_gpu) {
    gkyl_cu_free(gkns->omega_cfl);
  } else {
    gkyl_free(gkns->omega_cfl);
  }

  gk_neut_species_source_release(app, &gkns->src);

  // Release integrated mom data.
  gk_neut_species_moment_release(app, &gkns->integ_moms);

  // Release integrated mom diag data.
  gkyl_dynvec_release(gkns->integ_diag);

  if (app->use_gpu) {
    gkyl_cu_free(gkns->red_integ_diag);
    gkyl_cu_free(gkns->red_integ_diag_global);
  } else {
    gkyl_free(gkns->red_integ_diag);
    gkyl_free(gkns->red_integ_diag_global);
  }

  // Copy BCs are allocated by default. Need to free.
  for (int d = 0; d < app->cdim; ++d) {
    if (gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_RECYCLE) {
      gk_neut_species_recycle_release(app, &gkns->bc_recycle_lo);
    } else if ((gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
               (gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
               (gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
               (gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
      gkyl_bc_basic_release(gkns->bc_lo[d]);
    }
    if (gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_RECYCLE) {
      gk_neut_species_recycle_release(app, &gkns->bc_recycle_up);
    } else if ((gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
               (gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
               (gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
               (gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
      gkyl_bc_basic_release(gkns->bc_up[d]);
    }
  }
}

static void
gk_neut_species_kinetic_release(const gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns)
{
  // Release resources for kinetic neutral species.
  gkyl_msgpack_map_elem_release(gkns->io_meta_basic_len, gkns->io_meta_basic);
  gkyl_msgpack_map_elem_release(gkns->io_meta_conf_len, gkns->io_meta_conf);
  gkyl_msgpack_map_elem_release(gkns->io_meta_phase_len, gkns->io_meta_phase);

  gkyl_array_release(gkns->f);
  if (gkns->info.init_from_file.type == 0) {
    gk_neut_species_projection_release(app, &gkns->proj_init);
  }
  gkyl_comm_release(gkns->comm);

  if (app->use_gpu) {
    gkyl_array_release(gkns->f_host);
    gkyl_cu_free(gkns->basis_on_dev);
  }

  gkyl_velocity_map_release(gkns->vel_map);
  gkyl_vlasov_velocity_map_release(gkns->vlasov_vel_map);
  gkyl_vlasov_position_map_release(gkns->vlasov_pos_map);

  // Release moment data.
  gk_neut_species_moment_release(app, &gkns->m0);
  for (int i = 0; i < gkns->info.num_diag_moments; ++i) {
    gk_neut_species_moment_release(app, &gkns->moms[i]);
  }
  gkyl_free(gkns->moms);

  gk_neut_species_bgk_release(app, &gkns->bgk);

  // Free boundary flux memory.
  gk_neut_species_bflux_release(app, gkns, &gkns->bflux);

  gk_neut_species_lte_release(app, &gkns->lte);

  gk_neut_species_react_release(app, &gkns->react_neut);

  gkyl_array_release(gkns->gij);
  gkyl_array_release(gkns->g_ij);

  gkyl_array_release(gkns->hamil);
  if (app->use_gpu) {
    gkyl_array_release(gkns->hamil_host);
  }
  gkyl_array_release(gkns->conf_poisson_tensor);
  gkyl_array_release(gkns->conf_poisson_tensor_surf);
  gkyl_array_release(gkns->f_no_J);

  gk_neut_species_collisionless_release(app, &gkns->collisionless);

  // Free memory for the object that scales the species.
  gk_neut_species_scaling_release(app, &gkns->sca);

  gkns->release_is_static_func(app, gkns);
}

void
gk_neut_species_kinetic_init_dynamic(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns
)
{
  int cdim = app->cdim, vdim = gkns->info.vdim;
  int pdim = cdim + vdim;

  // allocate additional distribution function arrays for time stepping
  gkns->f1 = mkarr(app->use_gpu, gkns->basis.num_basis, gkns->local_ext.volume);
  gkns->fnew = mkarr(app->use_gpu, gkns->basis.num_basis, gkns->local_ext.volume);

  // Allocate cflrate (scalar array).
  gkns->cflrate = mkarr(app->use_gpu, 1, gkns->local_ext.volume);

  if (app->use_gpu) {
    gkns->omega_cfl = gkyl_cu_malloc(sizeof(double));
  } else {
    gkns->omega_cfl = gkyl_malloc(sizeof(double));
  }

  // Allocate data for integrated moments.
  gk_neut_species_moment_init(app, gkns, &gkns->integ_moms, GKYL_F_MOMENT_M0M1M2, true);

  // Allocate data for integrated diagnostics.
  if (app->use_gpu) {
    gkns->red_integ_diag = gkyl_cu_malloc(sizeof(double[gkns->integ_moms.num_mom]));
    gkns->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[gkns->integ_moms.num_mom]));
  } else {
    gkns->red_integ_diag = gkyl_malloc(sizeof(double[gkns->integ_moms.num_mom]));
    gkns->red_integ_diag_global = gkyl_malloc(sizeof(double[gkns->integ_moms.num_mom]));
  }
  // Allocate dynamic-vector to store all-reduced integrated moments.
  gkns->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, gkns->integ_moms.num_mom);
  gkns->is_first_integ_write_call = true;

  // Allocate buffer needed for BCs.
  long buff_sz = 0;
  for (int dir = 0; dir < cdim; ++dir) {
    long vol = GKYL_MAX2(gkns->local_lower_skin[dir].volume, gkns->local_upper_skin[dir].volume);
    buff_sz = buff_sz > vol ? buff_sz : vol;
  }
  gkns->bc_buffer = mkarr(app->use_gpu, gkns->basis.num_basis, buff_sz);
  // Buffer arrays for fixed function boundary conditions on distribution function.
  gkns->bc_buffer_lo_fixed = mkarr(app->use_gpu, gkns->basis.num_basis, buff_sz);
  gkns->bc_buffer_up_fixed = mkarr(app->use_gpu, gkns->basis.num_basis, buff_sz);

  gkns->bc_recycle_lo = (struct gk_recycle_wall){};
  gkns->bc_recycle_up = (struct gk_recycle_wall){};
  for (int d = 0; d < cdim; ++d) {
    // Reflection flips velocity component d. In the triad frame (e1, e2, bhat) the
    // velocity along z is component 2 for any cdim, so in 1x/2x it is not component d.
    if (cdim < 3 && d == cdim - 1 &&
        (gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT ||
         gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT)) {
      gkyl_exit(
        "gk_neut_species: reflecting BCs along z are only supported for 3x kinetic neutrals (in 1x/2x the velocity along z is the third triad component, not the one the reflection flips)."
      );
    }

    // Copy BCs by default.
    if (gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_RECYCLE) {
      // Initialize recycling object.
      gk_neut_species_recycle_init(
        app, &gkns->bc_recycle_lo, d, GKYL_LOWER_EDGE, &gkns->lower_bc[d].emission, gkns,
        app->use_gpu
      );
    } else if ((gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
               (gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
               (gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
               (gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
      enum gkyl_bc_basic_type bctype =
        gkyl_gyrokinetic_translate_bc_basic_type(gkns->lower_bc[d].type);
      gkns->bc_lo[d] = gkyl_bc_basic_new(
        d, GKYL_LOWER_EDGE, bctype, gkns->basis_on_dev, &gkns->local_lower_skin[d],
        &gkns->local_lower_ghost[d], gkns->f->ncomp, app->cdim, app->use_gpu
      );

      if (gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) {
        // Fill the buffer used for BCs.
        struct gk_proj gk_proj_bc_lo;
        gk_neut_species_projection_init(app, gkns, gkns->lower_bc[d].projection, &gk_proj_bc_lo);
        gk_neut_species_projection_calc(
          app, gkns, &gk_proj_bc_lo, gkns->f1, 0.0
        ); // Temporarily use f1.
        gkyl_bc_basic_buffer_fixed_func(gkns->bc_lo[d], gkns->bc_buffer_lo_fixed, gkns->f1);
        gkyl_array_clear(gkns->f1, 0.0);
        gk_neut_species_projection_release(app, &gk_proj_bc_lo);
      }
    }

    if (gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_RECYCLE) {
      // Initialize recycling object.
      gk_neut_species_recycle_init(
        app, &gkns->bc_recycle_up, d, GKYL_UPPER_EDGE, &gkns->upper_bc[d].emission, gkns,
        app->use_gpu
      );
    } else if ((gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_COPY) ||
               (gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_ABSORB) ||
               (gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_REFLECT) ||
               (gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC)) {
      // Upper BC updater. Copy BCs by default.
      enum gkyl_bc_basic_type bctype =
        gkyl_gyrokinetic_translate_bc_basic_type(gkns->upper_bc[d].type);
      gkns->bc_up[d] = gkyl_bc_basic_new(
        d, GKYL_UPPER_EDGE, bctype, gkns->basis_on_dev, &gkns->local_upper_skin[d],
        &gkns->local_upper_ghost[d], gkns->f->ncomp, app->cdim, app->use_gpu
      );

      if (gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_FIXED_FUNC) {
        // Fill the buffer used for BCs.
        struct gk_proj gk_proj_bc_up;
        gk_neut_species_projection_init(app, gkns, gkns->upper_bc[d].projection, &gk_proj_bc_up);
        gk_neut_species_projection_calc(
          app, gkns, &gk_proj_bc_up, gkns->f1, 0.0
        ); // Temporarily use f1.
        gkyl_bc_basic_buffer_fixed_func(gkns->bc_up[d], gkns->bc_buffer_up_fixed, gkns->f1);
        gkyl_array_clear(gkns->f1, 0.0);
        gk_neut_species_projection_release(app, &gk_proj_bc_up);
      }
    }
  }

  // Set function pointers
  gkns->rhs_func = gk_neut_species_kinetic_rhs_dynamic;
  gkns->rhs_implicit_func = gk_neut_species_kinetic_rhs_implicit_dynamic;
  gkns->bc_func = gk_neut_species_kinetic_apply_bc_dynamic;
  gkns->release_func = gk_neut_species_kinetic_release;
  gkns->release_is_static_func = gk_neut_species_kinetic_release_dynamic;
  gkns->step_f_func = gk_neut_species_step_f_dynamic;
  gkns->combine_func = gk_neut_species_combine_dynamic;
  gkns->copy_func = gk_neut_species_copy_range_dynamic;
  gkns->write_func = gk_neut_species_write_dynamic;
  gkns->write_mom_func = gk_neut_species_write_mom_dynamic;
  gkns->calc_integrated_mom_func = gk_neut_species_calc_integrated_mom_dynamic;
  gkns->write_integrated_mom_func = gk_neut_species_write_integrated_mom_dynamic;
  gkns->report_n_iter_corr_func = gk_neut_species_n_iter_corr_enabled;
}

static void
gk_neut_species_kinetic_init_static(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns
)
{
  // Set pointers for RK methods.
  gkns->f1 = gkns->f;
  gkns->fnew = gkns->f;

  // Set function pointers
  gkns->rhs_func = gk_neut_species_rhs_static;
  gkns->rhs_implicit_func = gk_neut_species_rhs_implicit_static;
  gkns->bc_func = gk_neut_species_apply_bc_static;
  gkns->release_func = gk_neut_species_kinetic_release;
  gkns->release_is_static_func = gk_neut_species_release_static;
  gkns->step_f_func = gk_neut_species_step_f_static;
  gkns->combine_func = gk_neut_species_combine_static;
  gkns->copy_func = gk_neut_species_copy_range_static;
  gkns->write_func = gk_neut_species_write_init_only;
  gkns->write_mom_func = gk_neut_species_write_mom_static;
  gkns->calc_integrated_mom_func = gk_neut_species_calc_integrated_mom_static;
  gkns->write_integrated_mom_func = gk_neut_species_write_integrated_mom_static;
  gkns->report_n_iter_corr_func = gk_neut_species_n_iter_corr_disabled;
}

void
gk_neut_species_kinetic_read_distf_from_file(
  struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,

  struct gkyl_gyrokinetic_ic_import inp, struct gkyl_array *fout
)
{
  // Read a distribution function from a file onto fout. Intended options include importing:
  //   1) ICs with same grid.
  //   2) ICs one dimensionality lower (e.g. 2x2v for 3x2v sim).
  //   3) ICs with same grid extents but different resolution (NYI).

  struct gkyl_rect_grid grid = gkns->grid;
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
      //   - 1x3v for a 2x3v sim.
      //   - 2x3v for a 3x3v sim.
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
        gkns->basis.id, gkyl_msgpack_map_elem_get_string(elem_list_len, elem_list, "basis_type")
      ) == 0
    );
    assert(poly_order == gkyl_msgpack_map_elem_get_uint(elem_list_len, elem_list, "poly_order"));
    gkyl_msgpack_map_elem_release_string(elem_list_len, elem_list, "basis_type");

    gkyl_grid_sub_array_header_release(&hdr);
  }

  // Donor basis.
  struct gkyl_basis basis_do;
  // Basis is tensor for p=1 and ser for p>1
  if (poly_order > 1) {
    gkyl_cart_modal_serendip(&basis_do, pdim_do, poly_order);
  } else if (poly_order == 1) {
    gkyl_cart_modal_tensor(&basis_do, pdim_do, poly_order); // for canonical PB
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
  struct gkyl_comm *comm_do = gkyl_comm_split_comm(gkns->comm, 0, decomp_do);

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

  bool scale_by_jacobgeo = false;
  with_file(fp, inp.jacobtot_inv_file_name, "r")
  {
    // Set up configuration space donor grid and basis
    struct gkyl_rect_grid conf_grid_do;
    gkyl_rect_grid_init(&conf_grid_do, cdim_do, grid_do.lower, grid_do.upper, grid_do.cells);
    // Create configuration space global ranges
    struct gkyl_range conf_local_ext_do, conf_local_do, conf_global_ext_do, conf_global_do;
    gkyl_create_grid_ranges(&conf_grid_do, ghost_do, &conf_global_ext_do, &conf_global_do);
    // Create configuration space local ranges
    struct gkyl_rect_decomp *conf_decomp_do =
      gkyl_rect_decomp_new_from_cuts(cdim_do, cuts_do, &conf_global_do);
    gkyl_create_ranges(
      &conf_decomp_do->ranges[my_rank], ghost_do, &conf_local_ext_do, &conf_local_do
    );
    // Create a configuration space basis.
    struct gkyl_basis conf_basis_do;
    gkyl_cart_modal_serendip(&conf_basis_do, cdim_do, poly_order);
    // Array for Jacobian inverse
    struct gkyl_array *jacobgeo_inv_do_host =
      mkarr(false, conf_basis_do.num_basis, conf_local_ext_do.volume);
    rstat.io_status = gkyl_comm_array_read(
      comm_do, &conf_grid_do, &conf_local_do, jacobgeo_inv_do_host, inp.jacobtot_inv_file_name
    );
    gkyl_dg_mul_conf_phase_op_range(
      &conf_basis_do, &basis_do, fdo_host, jacobgeo_inv_do_host, fdo_host, &conf_local_ext_do,
      &local_ext_do
    );
    gkyl_array_release(jacobgeo_inv_do_host);
    gkyl_rect_decomp_release(conf_decomp_do);
    scale_by_jacobgeo = true;
  }

  if (app->use_gpu) {
    gkyl_array_copy(fdo, fdo_host);
  }

  if (pdim_do == pdim - 1) {
    struct gkyl_translate_dim *transdim =
      gkyl_translate_dim_new(cdim_do, basis_do, cdim, gkns->basis, -1, GKYL_NO_EDGE, app->use_gpu);
    gkyl_translate_dim_advance(transdim, &local_do, &gkns->local, fdo, 1, fout);
    gkyl_translate_dim_release(transdim);
  } else {
    if (same_res) {
      gkyl_array_copy(fout, fdo);
    } else {
      // Interpolate the donor distribution to the target grid.
      struct gkyl_dg_interpolate *interp = gkyl_dg_interpolate_new(
        app->cdim, &gkns->basis, &grid_do, &grid, &local_do, &gkns->local, ghost_do, app->use_gpu
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
                                   gkyl_array_acquire(xfac_ho);
    gkyl_proj_on_basis_advance(proj_conf_scale, 0.0, &app->local, xfac_ho);
    gkyl_array_copy(xfac, xfac_ho);
    gkyl_dg_mul_conf_phase_op_range(
      &app->basis, &gkns->basis, fout, xfac, fout, &app->local, &gkns->local
    );
    gkyl_proj_on_basis_release(proj_conf_scale);
    gkyl_array_release(xfac_ho);
    gkyl_array_release(xfac);
  }

  // Multiply fout by the Jacobian.
  if (scale_by_jacobgeo) {
    gkyl_dg_mul_conf_phase_op_range(
      &app->basis, &gkns->basis, fout, app->gk_geom->geo_int.jacobgeo, fout, &app->local,
      &gkns->local
    );
  }

  gkyl_rect_decomp_release(decomp_do);
  gkyl_comm_release(comm_do);
  gkyl_array_release(fdo);
  gkyl_array_release(fdo_host);
}

void
gk_neut_species_kinetic_file_import_init(
  struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns,
  struct gkyl_gyrokinetic_ic_import inp
)
{
  gk_neut_species_kinetic_read_distf_from_file(app, gkns, inp, gkns->f);
}

static void
gkyl_array_move_comp(
  struct gkyl_array *out, int cout, struct gkyl_array *in, int cin, struct gkyl_array *tmp
)
{
  // Move the 'cin' scalar field from the multi-component (multi-scalar field) array
  // 'in' to 'cin' component of the the multi-component 'out'. This requires a
  // temporary array 'tmp'. This is really a work around for the fact that some
  // of our array offset ops don't work as we'd like for multi-component arrays
  // with the same number of components.
  gkyl_array_set_offset(tmp, 1.0, in, cin);
  gkyl_array_set_offset(out, 1.0, tmp, cout);
}

void
gk_neut_species_kinetic_init(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns
)
{
  gkns->is_fluid = false; // Kinetic neutrals.

  int cdim = app->cdim, vdim = gkns->info.vdim;
  int pdim = cdim + vdim;

  assert(vdim > 0); // Ensure user provided vdim in input file.

  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  int cells_vel[GKYL_MAX_DIM], ghost_vel[GKYL_MAX_DIM];
  double lower_vel[GKYL_MAX_DIM], upper_vel[GKYL_MAX_DIM];

  for (int d = 0; d < cdim; ++d) {
    cells[d] = gk->cells[d];
    lower[d] = gk->lower[d];
    upper[d] = gk->upper[d];
    ghost[d] = 1;
  }
  for (int d = 0; d < vdim; ++d) {
    // full phase space grid
    cells[cdim + d] = gkns->info.cells[d];
    lower[cdim + d] = gkns->info.lower[d];
    upper[cdim + d] = gkns->info.upper[d];
    ghost[cdim + d] = 0; // no ghost-cells in velocity space

    // only velocity space
    cells_vel[d] = gkns->info.cells[d];
    lower_vel[d] = gkns->info.lower[d];
    upper_vel[d] = gkns->info.upper[d];
    ghost_vel[d] = 0; // no ghost-cells in velocity space
  }

  // Tensor p=1 hybrid basis: p=1 in configuration space tensored with p=2 in
  // velocity space. The p=2 velocity basis represents the quadratic H = v^2/2
  // exactly; the pure p=1 basis cannot, and its sparse-Hamiltonian moment and
  // energy diagnostics are wrong at the O(1) level.
  assert(app->poly_order == 1); // Neutral species are p=1 only.
  gkyl_cart_modal_hybrid(&gkns->basis, cdim, vdim);
  if (cdim < 3 && !app->gk_geom->has_dxdz_grad) {
    gkyl_exit(
      "gk_neut_species: kinetic neutrals on a 1x/2x geometry need the tangent gradients (geo_int_dxdz_grad), which these geometry files do not carry. Regenerate them with this version from a run that includes kinetic neutrals."
    );
  }

  if (app->use_gpu) {
    // Allocate device basis if we are using GPUs.
    gkns->basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    gkyl_cart_modal_hybrid_cu_dev(gkns->basis_on_dev, cdim, vdim);
  } else {
    gkns->basis_on_dev = &gkns->basis;
  }

  // Full phase space grid.
  gkyl_rect_grid_init(&gkns->grid, pdim, lower, upper, cells);
  gkyl_create_grid_ranges(&gkns->grid, ghost, &gkns->global_ext, &gkns->global);

  // Velocity space grid.
  gkyl_rect_grid_init(&gkns->grid_vel, vdim, lower_vel, upper_vel, cells_vel);
  gkyl_create_grid_ranges(&gkns->grid_vel, ghost_vel, &gkns->local_ext_vel, &gkns->local_vel);

  // Phase-space communicator.
  gkns->comm = gkyl_comm_extend_comm(app->comm, &gkns->local_vel);

  // Create local and local_ext from app local range.
  struct gkyl_range local;
  // local = conf-local X local_vel.
  gkyl_range_ten_prod(&local, &app->local, &gkns->local_vel);
  gkyl_create_ranges(&local, ghost, &gkns->local_ext, &gkns->local);

  // Model type is the triads (orthonormal-frame velocity) formulation to handle
  // curvilinear coordinates, with the separable H = v^2/2 (sparse velocity-space
  // Hamiltonian expansion). GK coordinates are too sheared/distorted for the
  // canonical-momentum formulation at practical resolutions.
  gkns->model_id = GKYL_MODEL_TRIAD;
  gkns->hamil_id = GKYL_HAMIL_VEL_SPARSE;
  // Field type is null as this species doesn't interact with fields.
  gkns->field_id = GKYL_FIELD_NULL;

  // Velocity space mapping.
  assert(gkns->info.mapc2p.mapping == 0); // mapped v-space not implemented for neutrals yet.
  gkns->vel_map = gkyl_velocity_map_new(
    gkns->info.mapc2p, gkns->grid, gkns->grid_vel, gkns->local, gkns->local_ext, gkns->local_vel,
    gkns->local_ext_vel, app->use_gpu
  );

  // Identity Vlasov velocity map: the Vlasov flux, moment and LTE updaters
  // always evaluate the velocity coordinate from the stored map. The velocity
  // basis is the p=2 tensor basis of the hybrid, which also selects the C^1
  // cubic velocity map and the p=2 velocity representation of the Hamiltonian.
  gkyl_cart_modal_tensor(&gkns->basis_vel, vdim, 2);
  struct gkyl_vlasov_velocity_map_inp inp_vvmap[GKYL_MAX_CDIM] = {0};
  gkns->vlasov_vel_map = gkyl_vlasov_velocity_map_new(
    &gkns->grid_vel, &gkns->local_vel, &gkns->basis_vel, inp_vvmap, false, app->use_gpu
  );

  // Identity Vlasov position map, required by the Vlasov flux updaters.
  struct gkyl_vlasov_position_map_inp inp_vpmap[GKYL_MAX_CDIM] = {0};
  gkns->vlasov_pos_map = gkyl_vlasov_position_map_new(
    &app->grid, &app->local, &app->local_ext, &app->basis, inp_vpmap, app->use_gpu
  );

  // Keep a copy of num_periodic_dir and periodic_dirs in species so we can
  // add the parallel direction in case TS BCs are needed.
  gkns->num_periodic_dir = app->num_periodic_dir;
  for (int d = 0; d < gkns->num_periodic_dir; ++d) {
    gkns->periodic_dirs[d] = app->periodic_dirs[d];
  }

  for (int d = 0; d < app->cdim; ++d) {
    gkns->bc_is_np[d] = true;
  }
  for (int d = 0; d < gkns->num_periodic_dir; ++d) {
    gkns->bc_is_np[gkns->periodic_dirs[d]] = false;
  }

  // Store the BCs from the input file.
  for (int d = 0; d < app->cdim; ++d) {
    if (gkns->bc_is_np[d]) {
      struct gkyl_gyrokinetic_bc *bc_lo =
        gk_fetch_bc_with_dir_edge(gkns->info.bcs, 2 * app->cdim, d, GKYL_LOWER_EDGE);
      if (bc_lo != 0) {
        gkns->lower_bc[d] = *bc_lo;
      } else {
        gkns->lower_bc[d].type = GKYL_BC_GK_SKIP;
      }

      struct gkyl_gyrokinetic_bc *bc_up =
        gk_fetch_bc_with_dir_edge(gkns->info.bcs, 2 * app->cdim, d, GKYL_UPPER_EDGE);
      if (bc_up != 0) {
        gkns->upper_bc[d] = *bc_up;
      } else {
        gkns->upper_bc[d].type = GKYL_BC_GK_SKIP;
      }
    } else {
      gkns->lower_bc[d].type = GKYL_BC_GK_SPECIES_PERIODIC;
      gkns->upper_bc[d].type = GKYL_BC_GK_SPECIES_PERIODIC;
    }
  }

  // Species properties metadata.
  struct gkyl_msgpack_map_elem io_meta_sprop[] = {
    {.key = "mass", .elem_type = GKYL_MP_DOUBLE, .dval = gkns->info.mass},
    {.key = "charge", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
    {.key = "vdim", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = gkns->info.vdim}
  };

  // Metadata for integrated quantities.
  const struct gkyl_msgpack_map_elem *io_meta_basic_union[] = {app->io_meta_basic, io_meta_sprop};
  int io_meta_basic_union_len[] = {
    app->io_meta_basic_len, sizeof(io_meta_sprop) / sizeof(io_meta_sprop[0])
  };
  gkns->io_meta_basic = gkyl_msgpack_map_elem_union(
    sizeof(io_meta_basic_union) / sizeof(io_meta_basic_union[0]), io_meta_basic_union_len,
    io_meta_basic_union, &gkns->io_meta_basic_len
  );

  // Metadata for conf-space quantities.
  struct gkyl_msgpack_map_elem io_meta_conf[] = {
    {.key = "value_form", .elem_type = GKYL_MP_STRING, .cval = "modal"},
    {.key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = app->basis.poly_order},
    {.key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = app->basis.id},
    {.key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
    {.key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0}
  };
  const struct gkyl_msgpack_map_elem *io_meta_conf_union[] = {
    app->io_meta_basic, io_meta_sprop, io_meta_conf
  };
  int io_meta_conf_union_len[] = {
    app->io_meta_basic_len, sizeof(io_meta_sprop) / sizeof(io_meta_sprop[0]),
    sizeof(io_meta_conf) / sizeof(io_meta_conf[0])
  };
  gkns->io_meta_conf = gkyl_msgpack_map_elem_union(
    sizeof(io_meta_conf_union) / sizeof(io_meta_conf_union[0]), io_meta_conf_union_len,
    io_meta_conf_union, &gkns->io_meta_conf_len
  );

  // Metadata for phase-space quantities.
  struct gkyl_msgpack_map_elem io_meta_phase[] = {
    {.key = "value_form", .elem_type = GKYL_MP_STRING, .cval = "modal"},
    {.key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = gkns->basis.poly_order},
    {.key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = gkns->basis.id},
    {.key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
    {.key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0}
  };
  const struct gkyl_msgpack_map_elem *io_meta_phase_union[] = {
    app->io_meta_basic, io_meta_sprop, io_meta_phase
  };
  int io_meta_phase_union_len[] = {
    app->io_meta_basic_len, sizeof(io_meta_sprop) / sizeof(io_meta_sprop[0]),
    sizeof(io_meta_phase) / sizeof(io_meta_phase[0])
  };
  gkns->io_meta_phase = gkyl_msgpack_map_elem_union(
    sizeof(io_meta_phase_union) / sizeof(io_meta_phase_union[0]), io_meta_phase_union_len,
    io_meta_phase_union, &gkns->io_meta_phase_len
  );

  // Allocate distribution function array for initialization and I/O.
  gkns->f = mkarr(app->use_gpu, gkns->basis.num_basis, gkns->local_ext.volume);

  gkns->f_host = gkns->f;
  if (app->use_gpu) {
    gkns->f_host = mkarr(false, gkns->basis.num_basis, gkns->local_ext.volume);
  }

  if (app->cdim < 3) {
    // Reorganize g_ij and gij as done in calculation of Hamiltonian to
    // compute momentum and temperature.
    gkns->gij = mkarr(
      app->use_gpu, app->gk_geom->geo_int.gij_neut->ncomp, app->gk_geom->geo_int.gij_neut->size
    );
    gkns->g_ij = mkarr(
      app->use_gpu, app->gk_geom->geo_int.g_ij_neut->ncomp, app->gk_geom->geo_int.g_ij_neut->size
    );

    // Reorganize the metric tensor so ignorable coordinates are last.
    int metric_reorg_idxs_1x[] = {5, 2, 4, 0, 1, 3};
    int metric_reorg_idxs_2x[] = {0, 2, 1, 5, 4, 3};
    int *metric_reorg_idxs = app->cdim == 1 ? metric_reorg_idxs_1x : metric_reorg_idxs_2x;
    int num_basis_conf = app->basis.num_basis;
    struct gkyl_array *tmp_arr =
      mkarr(app->use_gpu, num_basis_conf, app->gk_geom->geo_int.gij_neut->size);
    for (int i = 0; i < 6; i++) {
      gkyl_array_move_comp(
        gkns->gij, i * num_basis_conf, app->gk_geom->geo_int.gij_neut,
        metric_reorg_idxs[i] * num_basis_conf, tmp_arr
      );
      gkyl_array_move_comp(
        gkns->g_ij, i * num_basis_conf, app->gk_geom->geo_int.g_ij_neut,
        metric_reorg_idxs[i] * num_basis_conf, tmp_arr
      );
    }
    gkyl_array_release(tmp_arr);
  } else {
    gkns->gij = gkyl_array_acquire(app->gk_geom->geo_int.gij);
    gkns->g_ij = gkyl_array_acquire(app->gk_geom->geo_int.g_ij);
  }

  // Hamiltonian: in the orthonormal triad frame H = v^2/2, a function of
  // velocity space only (sparse expansion), built exactly as for
  // GKYL_MODEL_DEFAULT in the Vlasov app.
  gkns->hamil_range = gkns->local_vel;
  gkns->hamil = mkarr(app->use_gpu, gkns->basis_vel.num_basis, gkns->local_vel.volume);
  // The inverse Hamiltonian output is only meaningful relativistically; the
  // kernels write it unconditionally, so hand them a scratch array.
  struct gkyl_array *hamil_inv_scratch =
    mkarr(app->use_gpu, gkns->basis_vel.num_basis, gkns->local_vel.volume);
  gkyl_dg_vlasov_calc_hamil(
    &gkns->grid_vel, &gkns->basis_vel, &gkns->local_vel, GKYL_MODEL_DEFAULT, gkns->vlasov_vel_map,
    gkns->hamil, hamil_inv_scratch, app->use_gpu
  );
  gkyl_array_release(hamil_inv_scratch);

  gkns->hamil_host = gkns->hamil;
  if (app->use_gpu) {
    gkns->hamil_host = mkarr(false, gkns->basis_vel.num_basis, gkns->local_vel.volume);
    gkyl_array_copy(gkns->hamil_host, gkns->hamil);
  }

  // Poisson tensor of the b-aligned triad, stored at the configuration-space
  // Gauss-Legendre nodes (interior nodes and face nodes), built on the host from
  // the GK geometry's nodal Cartesian tangents. On GPUs the geometry's nodal
  // arrays live on the device, so take host copies first and copy the result back.
  int num_nodes = 1 << cdim, num_face_nodes = 1 << (cdim - 1);
  struct gkyl_array *pt_vol_host = mkarr(false, num_nodes * 18, app->local_ext.volume);
  struct gkyl_array *pt_surf_host =
    mkarr(false, cdim * 2 * num_face_nodes * 3, app->local_ext.volume);
  struct gk_geometry *geom = app->gk_geom;
  struct gkyl_array *dxdz_nodal_host = mkarr(false, 9, geom->nrange_int.volume);
  struct gkyl_array *dxdz_grad_nodal_host = cdim < 3 ? mkarr(false, 27, geom->nrange_int.volume) :
                                                       NULL;
  struct gkyl_array *bcart_nodal_host = mkarr(false, 3, geom->nrange_int.volume);
  struct gkyl_array *dxdz_surf_nodal_host[3] = {NULL};
  gkyl_array_copy(dxdz_nodal_host, geom->geo_int.dxdz_nodal);
  if (cdim < 3) {
    gkyl_array_copy(dxdz_grad_nodal_host, geom->geo_int.dxdz_grad_nodal);
  }
  gkyl_array_copy(bcart_nodal_host, geom->geo_int.bcart_nodal);
  for (int dir = 0; dir < cdim; ++dir) {
    dxdz_surf_nodal_host[dir] = mkarr(false, 9, geom->nrange_surf[dir].volume);
    gkyl_array_copy(dxdz_surf_nodal_host[dir], geom->geo_surf[dir].dxdz_nodal);
  }
  // Deflated geometry: Poisson-tensor coordinates are the grid directions followed
  // by the removed ones, (x, z, y) in 2x and (z, x, y) in 1x.
  int dir_map_2x[3] = {0, 2, 1}, dir_map_1x[3] = {2, 0, 1};
  const int *dir_map = cdim == 3 ? NULL : (cdim == 2 ? dir_map_2x : dir_map_1x);
  gkyl_vlasov_triad_geom_nodal_from_tangents(
    &app->grid, &app->local, &app->local_ext, &app->global, app->comm, app->use_gpu,
    &geom->nrange_int, dxdz_nodal_host, dxdz_grad_nodal_host, bcart_nodal_host, geom->nrange_surf,
    dxdz_surf_nodal_host, dir_map, true, pt_vol_host, pt_surf_host
  );
  gkns->conf_poisson_tensor = mkarr(app->use_gpu, pt_vol_host->ncomp, app->local_ext.volume);
  gkns->conf_poisson_tensor_surf = mkarr(app->use_gpu, pt_surf_host->ncomp, app->local_ext.volume);
  gkyl_array_copy(gkns->conf_poisson_tensor, pt_vol_host);
  gkyl_array_copy(gkns->conf_poisson_tensor_surf, pt_surf_host);
  gkyl_array_release(pt_vol_host);
  gkyl_array_release(pt_surf_host);
  gkyl_array_release(dxdz_nodal_host);
  gkyl_array_release(dxdz_grad_nodal_host);
  gkyl_array_release(bcart_nodal_host);
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_array_release(dxdz_surf_nodal_host[dir]);
  }

  // Distribution function with the velocity-space Jacobian divided out,
  // input to the velocity-space flux updater.
  gkns->f_no_J = mkarr(app->use_gpu, gkns->basis.num_basis, gkns->local_ext.volume);

  // Initialize the collisionless solver.
  gkns->collisionless = (struct gk_collisionless){};
  gk_neut_species_collisionless_init(app, gkns, &gkns->collisionless);

  // Create skin/ghost ranges.
  for (int dir = 0; dir < cdim; ++dir) {
    gkyl_skin_ghost_ranges(
      &gkns->local_lower_skin[dir], &gkns->local_lower_ghost[dir], dir, GKYL_LOWER_EDGE,
      &gkns->local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &gkns->local_upper_skin[dir], &gkns->local_upper_ghost[dir], dir, GKYL_UPPER_EDGE,
      &gkns->local_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &gkns->global_lower_skin[dir], &gkns->global_lower_ghost[dir], dir, GKYL_LOWER_EDGE,
      &gkns->global_ext, ghost
    );
    gkyl_skin_ghost_ranges(
      &gkns->global_upper_skin[dir], &gkns->global_upper_ghost[dir], dir, GKYL_UPPER_EDGE,
      &gkns->local_ext, ghost
    );
  }

  if (gkns->info.init_from_file.type == 0) {
    // Initialize projection routine for initial conditions.
    gk_neut_species_projection_init(app, gkns, gkns->info.projection, &gkns->proj_init);
  } else {
    // Read initial condition from file.
    gk_neut_species_kinetic_file_import_init(app, gkns, gkns->info.init_from_file);
  }

  // Allocate object for computing number .density
  gk_neut_species_moment_init(app, gkns, &gkns->m0, GKYL_F_MOMENT_M0, false);

  // Allocate objects for computing diagnostic moments.
  int ndm = gkns->info.num_diag_moments;
  gkns->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m = 0; m < ndm; ++m) {
    gk_neut_species_moment_init(app, gkns, &gkns->moms[m], gkns->info.diag_moments[m], false);
  }

  // Initialize boundary fluxes.
  gkns->bflux = (struct gk_boundary_fluxes){};
  // Additional bflux moments to step in time.
  struct gkyl_phase_diagnostics_inp add_bflux_moms_inp = (struct gkyl_phase_diagnostics_inp){};
  // Set the operation type for the bflux app.
  enum gkyl_species_bflux_type bflux_type = GK_SPECIES_BFLUX_NONE;
  if (gkns->info.boundary_flux_diagnostics.num_diag_moments > 0 ||
      gkns->info.boundary_flux_diagnostics.num_integrated_diag_moments > 0) {
    bflux_type = GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS_DIAGS;
  } else {
    // Set bflux_type to
    //   - GK_SPECIES_BFLUX_CALC_FLUX to only put bfluxes in ghost cells of rhs.
    //   - GK_SPECIES_BFLUX_CALC_FLUX_STEP_MOMS to calc bfluxes and step its moments.
    // The latter also requires that you place the moment you desire in add_bflux_moms_inp below.
  }
  // Introduce new moments into moms_inp if needed.
  gk_neut_species_bflux_init(app, gkns, &gkns->bflux, bflux_type, add_bflux_moms_inp);

  // Initialize a Maxwellian/LTE (local thermodynamic equilibrium) projection routine
  // Projection routine optionally corrects all the Maxwellian/LTE moments
  // This routine is utilized by both reactions and BGK collisions
  gkns->lte = (struct gk_lte){};
  bool correct_all_moms = gkns->info.correct.correct_all_moms;
  int max_iter = gkns->info.correct.max_iter > 0 ? gkns->info.correct.max_iter : 50;
  double iter_eps = gkns->info.correct.iter_eps > 0 ? gkns->info.correct.iter_eps : 1e-10;
  bool use_last_converged = gkns->info.correct.use_last_converged;
  struct correct_all_moms_inp corr_inp = {
    .correct_all_moms = correct_all_moms,
    .max_iter = max_iter,
    .iter_eps = iter_eps,
    .use_last_converged = use_last_converged,
  };
  gk_neut_species_lte_init(app, gkns, &gkns->lte, corr_inp);

  // Initialize elastic collisions.
  gkns->bgk = (struct gk_bgk_collisions){};
  gk_neut_species_bgk_init(app, gkns, &gkns->bgk);

  // Initialize the object that scales the species according to a balance
  // between recycling and reactions (meant for fluid neutrals for now).
  gkns->sca = (struct gk_scaling){};
  gk_neut_species_scaling_init(app, gkns, &gkns->sca);

  // Initialize reactions with charged species.
  gkns->react_neut = (struct gk_react){};
  gk_neut_species_react_init(app, gkns, gkns->info.react_neut, &gkns->react_neut);

  // Initialize empty structs. New methods will fill them if specified.
  gkns->src = (struct gk_source){};

  if (gkns->info.is_static) {
    gk_neut_species_kinetic_init_static(gk, app, gkns);
  } else {
    gk_neut_species_kinetic_init_dynamic(gk, app, gkns);
  }
}
