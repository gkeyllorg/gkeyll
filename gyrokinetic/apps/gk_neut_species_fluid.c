#include <gkyl_gk_neut_species_priv.h>

static double
gk_neut_species_fluid_rhs_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms
)
{
  double omega_cfl = 1 / DBL_MAX;
  gkyl_array_clear(gkns->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // Collisionless terms.
  struct timespec wst = gkyl_wall_clock();
  // Not ready.
  app->stat.neut_species_collisionless_tm += gkyl_time_diff_now_sec(wst);

  // Compute volume-integrated reactions in sca.
  gk_neut_species_scaling_rhs(app, gkns, &gkns->sca, fin, rhs);

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
gk_neut_species_fluid_rhs_implicit_dynamic(
  gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns, const struct gkyl_array *fin,
  struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt
)
{
  double omega_cfl = 1 / DBL_MAX;
  gkyl_array_clear(gkns->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // No implicit terms yet.

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
gk_neut_species_fluid_release_dynamic(
  const gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns
)
{
  // Release memory allocated for dynamic neutrals.
  gkyl_array_release(gkns->cflrate);

  if (app->use_gpu) {
    gkyl_cu_free(gkns->omega_cfl);
  } else {
    gkyl_free(gkns->omega_cfl);
  }

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
}

static void
gk_neut_species_fluid_release(const gkyl_gyrokinetic_app *app, const struct gk_neut_species *gkns)
{
  // Release resources for fluid neutral species.
  gkyl_msgpack_map_elem_release(gkns->io_meta_basic_len, gkns->io_meta_basic);
  gkyl_msgpack_map_elem_release(gkns->io_meta_conf_len, gkns->io_meta_conf);

  gkyl_array_release(gkns->f);
  gkyl_array_release(gkns->f1);
  gkyl_array_release(gkns->fnew);
  gkyl_array_release(gkns->f_host);

  if (gkns->info.init_from_file.type == 0) {
    gk_neut_species_projection_release(app, &gkns->proj_init);
  }
  gkyl_comm_release(gkns->comm);

  for (int i = 0; i < gkns->info.num_diag_moments; ++i) {
    gk_neut_species_moment_release(app, &gkns->moms[i]);
  }
  gkyl_free(gkns->moms);

  gk_neut_species_bgk_release(app, &gkns->bgk);

  gk_neut_species_positivity_release(app, &gkns->positivity);

  gk_neut_species_react_release(app, &gkns->react_neut);

  // Free boundary flux memory.
  gk_neut_species_bflux_release(app, gkns, &gkns->bflux);

  gk_neut_species_lte_release(app, &gkns->lte);

  // Free memory for the object that scales the species according to a balance
  // between recycling and reactions.
  gk_neut_species_scaling_release(app, &gkns->sca);

  gkns->release_is_static_func(app, gkns);
}

static void
gk_neut_species_fluid_init_dynamic(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns
)
{
  int cdim = app->cdim;

  // Allocate additional moment arrays for time stepping.
  gkns->f1 = mkarr(app->use_gpu, gkns->f->ncomp, gkns->f->size);
  gkns->fnew = mkarr(app->use_gpu, gkns->f->ncomp, gkns->f->size);

  // Allocate cflrate (scalar array).
  gkns->cflrate = mkarr(app->use_gpu, 1, gkns->local_ext.volume);

  gkns->omega_cfl = app->use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));

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

  // Set function pointers
  gkns->rhs_func = gk_neut_species_fluid_rhs_dynamic;
  gkns->rhs_implicit_func = gk_neut_species_fluid_rhs_implicit_dynamic;
  gkns->bc_func = gk_neut_species_apply_bc_static; // Not ready.
  gkns->release_func = gk_neut_species_fluid_release;
  gkns->release_is_static_func = gk_neut_species_fluid_release_dynamic;
  gkns->step_f_func = gk_neut_species_step_f_dynamic;
  gkns->combine_func = gk_neut_species_combine_dynamic;
  gkns->copy_func = gk_neut_species_copy_range_dynamic;
  gkns->write_func = gk_neut_species_write_dynamic;
  gkns->write_mom_func =
    gk_neut_species_write_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  gkns->calc_integrated_mom_func =
    gk_neut_species_calc_integrated_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  gkns->write_integrated_mom_func =
    gk_neut_species_write_integrated_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  gkns->report_n_iter_corr_func = gk_neut_species_n_iter_corr_disabled;
}

static void
gk_neut_species_fluid_init_static(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns
)
{
  // Set pointers for RK methods.
  gkns->f1 = gkyl_array_acquire(gkns->f);
  gkns->fnew = gkyl_array_acquire(gkns->f);

  // Set function pointers
  gkns->rhs_func = gk_neut_species_rhs_static;
  gkns->rhs_implicit_func = gk_neut_species_rhs_implicit_static;
  gkns->bc_func = gk_neut_species_apply_bc_static;
  gkns->release_func = gk_neut_species_fluid_release;
  gkns->release_is_static_func = gk_neut_species_release_static;
  gkns->step_f_func = gk_neut_species_step_f_static;
  gkns->combine_func = gk_neut_species_combine_static;
  gkns->copy_func = gk_neut_species_copy_range_static;
  gkns->write_func = gk_neut_species_write_init_only;
  gkns->write_mom_func = gk_neut_species_write_mom_init_only;
  gkns->calc_integrated_mom_func = gk_neut_species_calc_integrated_mom_static;
  gkns->write_integrated_mom_func = gk_neut_species_write_integrated_mom_static;
  gkns->report_n_iter_corr_func = gk_neut_species_n_iter_corr_disabled;
}

void
gk_neut_species_fluid_init(
  struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *gkns
)
{
  gkns->is_fluid = true; // Fluid neutrals.
  assert(
    gkns->info.vdim == 0
  ); // Ensure user provided vdim=0 in input file, or didn't provide it at all.

  gkns->model_id = GKYL_MODEL_DEFAULT;
  gkns->field_id = GKYL_FIELD_NULL;

  int cdim = app->cdim;

  // Number of moments depends on eqn_type.
  gkns->num_moments = 5; // rho, rho*ux, rho*uy, rho*uy, totE

  // Use the same basis as conf-space.
  gkns->basis = app->basis;
  gkns->basis_on_dev = app->basis_on_dev;

  // Use the same grid as conf-space.
  gkns->grid = app->grid;
  gkns->global_ext = app->global_ext;
  gkns->global = app->global;

  // Use the same communicator as conf-space.
  gkns->comm = gkyl_comm_acquire(app->comm);

  // Use the same local range as conf-space.
  gkns->local_ext = app->local_ext;
  gkns->local = app->local;

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
  }

  // Species properties metadata.
  struct gkyl_msgpack_map_elem io_meta_sprop[] = {
    {.key = "mass", .elem_type = GKYL_MP_DOUBLE, .dval = gkns->info.mass},
    {.key = "charge", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0},
    {.key = "gas_gamma", .elem_type = GKYL_MP_DOUBLE, .dval = gkns->info.gas_gamma},
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
  gkns->io_meta_phase = gkns->io_meta_conf;
  gkns->io_meta_phase_len = gkns->io_meta_conf_len;

  // Allocate distribution function array for initialization and I/O.
  gkns->f = mkarr(app->use_gpu, gkns->num_moments * gkns->basis.num_basis, gkns->local_ext.volume);
  gkns->f_host = app->use_gpu ? mkarr(false, gkns->f->ncomp, gkns->f->size) :
                                gkyl_array_acquire(gkns->f);

  // Create skin/ghost ranges fir applying BCs. Only used for dynamic neutrals but included here to avoid
  // code duplication since the "ghost" array is needed.
  int ghost[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    ghost[d] = 1;
  }

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

  // Initialize projection routine for initial conditions.
  if (gkns->info.init_from_file.type == 0) {
    gk_neut_species_projection_init(app, gkns, gkns->info.projection, &gkns->proj_init);
  }

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
  gk_neut_species_bflux_init(app, gkns, &gkns->bflux, bflux_type, add_bflux_moms_inp);

  // Initialize a Maxwellian/LTE (local thermodynamic equilibrium) projection routine
  // Projection routine optionally corrects all the Maxwellian/LTE moments
  // This routine is utilized by both reactions and BGK collisions
  gkns->lte = (struct gk_lte){};
  struct correct_all_moms_inp corr_inp = {
    .correct_all_moms = false,
    .max_iter = 0,
    .iter_eps = 10,
    .use_last_converged = false,
  };
  gk_neut_species_lte_init(app, gkns, &gkns->lte, corr_inp);

  // Initialize the object that scales the species according to a balance
  // between recycling and reactions.
  gkns->sca = (struct gk_scaling){};
  gk_neut_species_scaling_init(app, gkns, &gkns->sca);

  // Initialize BGK collisions with null type (not applicable to fluids).
  gkns->bgk = (struct gk_bgk_collisions){};
  gkns->info.collisions.collision_id = 0;
  gk_neut_species_bgk_init(app, gkns, &gkns->bgk);

  // Initialize positivity enforcing operator with null type (NYI for fluids).
  gkns->positivity = (struct gk_positivity){};
  gk_neut_species_positivity_init(app, gkns, &gkns->positivity);

  // Initialize reactions with charged species (NYI for fluids).
  gkns->react_neut = (struct gk_react){};
  gkns->info.react_neut.num_react = 0;
  gk_neut_species_react_init(app, gkns, gkns->info.react_neut, &gkns->react_neut);

  gkns->src = (struct gk_source){};
  if (!gkns->info.is_static) {
    gk_neut_species_fluid_init_dynamic(gk, app, gkns);
  } else {
    gk_neut_species_fluid_init_static(gk, app, gkns);
  }
}
