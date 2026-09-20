#include <gkyl_gk_neut_species_priv.h>

static double
gk_neut_species_fluid_rhs_dynamic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms)
{
  double omega_cfl = 1/DBL_MAX;   
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);
  
  // Collisionless terms.
  struct timespec wst = gkyl_wall_clock();
  // Not ready.
  app->stat.neut_species_collisionless_tm += gkyl_time_diff_now_sec(wst);

  // Compute volume-integrated reactions in sca.
  gk_neut_species_scaling_rhs(app, species, &species->sca, fin, rhs);

  app->stat.n_neut_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);
  
  double omega_cfl_ho[1];
  if (app->use_gpu)
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  else
    omega_cfl_ho[0] = species->omega_cfl[0];

  omega_cfl = omega_cfl_ho[0];
  
  app->stat.neut_species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl/omega_cfl;
}

static double
gk_neut_species_fluid_rhs_implicit_dynamic(gkyl_gyrokinetic_app *app, struct gk_neut_species *species,
  const struct gkyl_array *fin, struct gkyl_array *rhs, struct gkyl_array **bflux_moms, double dt)
{ 
  double omega_cfl = 1/DBL_MAX;
  gkyl_array_clear(species->cflrate, 0.0);
  gkyl_array_clear(rhs, 0.0);

  // No implicit terms yet.

  gkyl_array_accumulate(gkyl_array_scale(rhs, dt), 1.0, fin);
  
  app->stat.n_neut_species_omega_cfl +=1;
  struct timespec tm = gkyl_wall_clock();
  gkyl_array_reduce_range(species->omega_cfl, species->cflrate, GKYL_MAX, &species->local);
  
  double omega_cfl_ho[1];
  if (app->use_gpu) {
    gkyl_cu_memcpy(omega_cfl_ho, species->omega_cfl, sizeof(double), GKYL_CU_MEMCPY_D2H);
  }
  else {
    omega_cfl_ho[0] = species->omega_cfl[0];
  }
  omega_cfl = omega_cfl_ho[0];
  
  app->stat.neut_species_omega_cfl_tm += gkyl_time_diff_now_sec(tm);
  return app->cfl/omega_cfl;
}

static void
gk_neut_species_fluid_release_dynamic(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *ns)
{
  // Release memory allocated for dynamic neutrals.
  gkyl_array_release(ns->cflrate);
  
  if (app->use_gpu) {
    gkyl_cu_free(ns->omega_cfl);
  }
  else {
    gkyl_free(ns->omega_cfl);
  }

  // Release integrated mom data.
  gk_neut_species_moment_release(app, &ns->integ_moms); 

  // Release integrated mom diag data.
  gkyl_dynvec_release(ns->integ_diag);
  
  if (app->use_gpu) {
    gkyl_cu_free(ns->red_integ_diag);
    gkyl_cu_free(ns->red_integ_diag_global);
  }
  else {
    gkyl_free(ns->red_integ_diag);
    gkyl_free(ns->red_integ_diag_global);
  }
}

static void
gk_neut_species_fluid_release(const gkyl_gyrokinetic_app* app, const struct gk_neut_species *ns)
{
  // Release resources for fluid neutral species.
  gkyl_msgpack_map_elem_release(ns->io_meta_basic_len, ns->io_meta_basic);
  gkyl_msgpack_map_elem_release(ns->io_meta_conf_len , ns->io_meta_conf );

  gkyl_array_release(ns->f);
  gkyl_array_release(ns->f1);
  gkyl_array_release(ns->fnew);
  gkyl_array_release(ns->f_host);

  if (ns->info.init_from_file.type == 0) {
    gk_neut_species_projection_release(app, &ns->proj_init);
  }
  gkyl_comm_release(ns->comm);

  for (int i=0; i<ns->info.num_diag_moments; ++i)
    gk_neut_species_moment_release(app, &ns->moms[i]);
  gkyl_free(ns->moms);

  gk_neut_species_bgk_release(app, &ns->bgk);

  gk_neut_species_positivity_release(app, &ns->positivity);

  gk_neut_species_react_release(app, &ns->react_neut);

  // Free boundary flux memory.
  gk_neut_species_bflux_release(app, ns, &ns->bflux);

  gk_neut_species_lte_release(app, &ns->lte);

  // Free memory for the object that scales the species according to a balance
  // between recycling and reactions.
  gk_neut_species_scaling_release(app, &ns->sca);

  ns->release_is_static_func(app, ns);
}

static void
gk_neut_species_fluid_init_dynamic(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns)
{
  int cdim = app->cdim;
  
  // Allocate additional moment arrays for time stepping.
  ns->f1 = mkarr(app->use_gpu, ns->f->ncomp, ns->f->size);
  ns->fnew = mkarr(app->use_gpu, ns->f->ncomp, ns->f->size);
  
  // Allocate cflrate (scalar array).
  ns->cflrate = mkarr(app->use_gpu, 1, ns->local_ext.volume);

  ns->omega_cfl = app->use_gpu? gkyl_cu_malloc(sizeof(double))
                              : gkyl_malloc(sizeof(double));

  // Allocate data for integrated moments.
  gk_neut_species_moment_init(app, ns, &ns->integ_moms, GKYL_F_MOMENT_M0M1M2, true);

  // Allocate data for integrated diagnostics.
  if (app->use_gpu) {
    ns->red_integ_diag = gkyl_cu_malloc(sizeof(double[ns->integ_moms.num_mom]));
    ns->red_integ_diag_global = gkyl_cu_malloc(sizeof(double[ns->integ_moms.num_mom]));
  } else {
    ns->red_integ_diag = gkyl_malloc(sizeof(double[ns->integ_moms.num_mom]));
    ns->red_integ_diag_global = gkyl_malloc(sizeof(double[ns->integ_moms.num_mom]));
  }
  // Allocate dynamic-vector to store all-reduced integrated moments.
  ns->integ_diag = gkyl_dynvec_new(GKYL_DOUBLE, ns->integ_moms.num_mom);
  ns->is_first_integ_write_call = true;

  // Set function pointers
  ns->rhs_func = gk_neut_species_fluid_rhs_dynamic;
  ns->rhs_implicit_func = gk_neut_species_fluid_rhs_implicit_dynamic;
  ns->bc_func = gk_neut_species_apply_bc_static; // Not ready.
  ns->release_func = gk_neut_species_fluid_release;
  ns->release_is_static_func = gk_neut_species_fluid_release_dynamic;
  ns->step_f_func = gk_neut_species_step_f_dynamic;
  ns->combine_func = gk_neut_species_combine_dynamic;
  ns->copy_func = gk_neut_species_copy_range_dynamic;
  ns->write_func = gk_neut_species_write_dynamic;
  ns->write_mom_func = gk_neut_species_write_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  ns->calc_integrated_mom_func = gk_neut_species_calc_integrated_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  ns->write_integrated_mom_func = gk_neut_species_write_integrated_mom_dynamic; // MF 2025/07/18: currently works for fluid too.
  ns->report_n_iter_corr_func = gk_neut_species_n_iter_corr_disabled;
}

static void
gk_neut_species_fluid_init_static(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *s)
{
  // Set pointers for RK methods.
  s->f1 = gkyl_array_acquire(s->f);
  s->fnew = gkyl_array_acquire(s->f);

  // Set function pointers
  s->rhs_func = gk_neut_species_rhs_static;
  s->rhs_implicit_func = gk_neut_species_rhs_implicit_static;
  s->bc_func = gk_neut_species_apply_bc_static;
  s->release_func = gk_neut_species_fluid_release;
  s->release_is_static_func = gk_neut_species_release_static;
  s->step_f_func = gk_neut_species_step_f_static;
  s->combine_func = gk_neut_species_combine_static;
  s->copy_func = gk_neut_species_copy_range_static;
  s->write_func = gk_neut_species_write_init_only;
  s->write_mom_func = gk_neut_species_write_mom_init_only;
  s->calc_integrated_mom_func = gk_neut_species_calc_integrated_mom_static;
  s->write_integrated_mom_func = gk_neut_species_write_integrated_mom_static;
  s->report_n_iter_corr_func = gk_neut_species_n_iter_corr_disabled;
}

void
gk_neut_species_fluid_init(struct gkyl_gk *gk, struct gkyl_gyrokinetic_app *app, struct gk_neut_species *ns)
{
  ns->is_fluid = true; // Fluid neutrals.
  assert(ns->info.vdim == 0); // Ensure user provided vdim=0 in input file, or didn't provide it at all.

  ns->model_id = GKYL_MODEL_DEFAULT;
  ns->field_id = GKYL_FIELD_NULL;

  int cdim = app->cdim;

  // Number of moments depends on eqn_type.
  ns->num_moments = 5; // rho, rho*ux, rho*uy, rho*uy, totE

  // Use the same basis as conf-space.
  ns->basis = app->basis;
  ns->basis_on_dev = app->basis_on_dev;
  
  // Use the same grid as conf-space.
  ns->grid = app->grid;
  ns->global_ext = app->global_ext;
  ns->global = app->global;
  
  // Use the same communicator as conf-space.
  ns->comm = gkyl_comm_acquire(app->comm);

  // Use the same local range as conf-space.
  ns->local_ext = app->local_ext;
  ns->local = app->local;

  // Keep a copy of num_periodic_dir and periodic_dirs in species so we can
  // add the parallel direction in case TS BCs are needed.
  ns->num_periodic_dir = app->num_periodic_dir;
  for (int d=0; d<ns->num_periodic_dir; ++d)
    ns->periodic_dirs[d] = app->periodic_dirs[d];

  for (int d=0; d<app->cdim; ++d) ns->bc_is_np[d] = true;
  for (int d=0; d<ns->num_periodic_dir; ++d)
    ns->bc_is_np[ns->periodic_dirs[d]] = false;

  // Store the BCs from the input file.
  for (int d=0; d<app->cdim; ++d) {
    struct gkyl_gyrokinetic_bc *bc_lo = gk_fetch_bc_with_dir_edge(ns->info.bcs, 2*app->cdim, d, GKYL_LOWER_EDGE);
    if (bc_lo != 0)
      ns->lower_bc[d] = *bc_lo;
    else
      ns->lower_bc[d].type = GKYL_BC_GK_SKIP;

    struct gkyl_gyrokinetic_bc *bc_up = gk_fetch_bc_with_dir_edge(ns->info.bcs, 2*app->cdim, d, GKYL_UPPER_EDGE);
    if (bc_up != 0)
      ns->upper_bc[d] = *bc_up;
    else
      ns->upper_bc[d].type = GKYL_BC_GK_SKIP;
  }

  // Species properties metadata.
  struct gkyl_msgpack_map_elem io_meta_sprop[] = {
    { .key = "mass", .elem_type = GKYL_MP_DOUBLE, .dval = ns->info.mass },
    { .key = "charge", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0 },
    { .key = "gas_gamma", .elem_type = GKYL_MP_DOUBLE, .dval = ns->info.gas_gamma },
    { .key = "vdim", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = ns->info.vdim },
  };

  // Metadata for integrated quantities.
  const struct gkyl_msgpack_map_elem *io_meta_basic_union[] = {app->io_meta_basic, io_meta_sprop};
  int io_meta_basic_union_len[] = {app->io_meta_basic_len, sizeof(io_meta_sprop)/sizeof(io_meta_sprop[0])};
  ns->io_meta_basic = gkyl_msgpack_map_elem_union(sizeof(io_meta_basic_union)/sizeof(io_meta_basic_union[0]),
    io_meta_basic_union_len, io_meta_basic_union, &ns->io_meta_basic_len);

  // Metadata for conf-space quantities.
  struct gkyl_msgpack_map_elem io_meta_conf[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = app->basis.poly_order },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = app->basis.id },
    { .key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0 },
    { .key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
  };
  const struct gkyl_msgpack_map_elem *io_meta_conf_union[] = {app->io_meta_basic, io_meta_sprop, io_meta_conf};
  int io_meta_conf_union_len[] = {app->io_meta_basic_len, sizeof(io_meta_sprop)/sizeof(io_meta_sprop[0]),
    sizeof(io_meta_conf)/sizeof(io_meta_conf[0])};
  ns->io_meta_conf = gkyl_msgpack_map_elem_union(sizeof(io_meta_conf_union)/sizeof(io_meta_conf_union[0]),
    io_meta_conf_union_len, io_meta_conf_union, &ns->io_meta_conf_len);

  // Metadata for phase-space quantities.
  ns->io_meta_phase = ns->io_meta_conf;
  ns->io_meta_phase_len = ns->io_meta_conf_len;

  // Allocate distribution function array for initialization and I/O.
  ns->f = mkarr(app->use_gpu, ns->num_moments*ns->basis.num_basis, ns->local_ext.volume);
  ns->f_host = app->use_gpu? mkarr(false, ns->f->ncomp, ns->f->size)
                           : gkyl_array_acquire(ns->f);

  // Create skin/ghost ranges fir applying BCs. Only used for dynamic neutrals but included here to avoid
  // code duplication since the "ghost" array is needed.
  int ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d)
    ghost[d] = 1;

  // Create skin/ghost ranges.
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&ns->local_lower_skin[dir], &ns->local_lower_ghost[dir],
      dir, GKYL_LOWER_EDGE, &ns->local_ext, ghost);
    gkyl_skin_ghost_ranges(&ns->local_upper_skin[dir], &ns->local_upper_ghost[dir],
      dir, GKYL_UPPER_EDGE, &ns->local_ext, ghost);
    gkyl_skin_ghost_ranges(&ns->global_lower_skin[dir], &ns->global_lower_ghost[dir],
      dir, GKYL_LOWER_EDGE, &ns->global_ext, ghost); 
    gkyl_skin_ghost_ranges(&ns->global_upper_skin[dir], &ns->global_upper_ghost[dir],
      dir, GKYL_UPPER_EDGE, &ns->local_ext, ghost);
  }

  // Initialize projection routine for initial conditions.
  if (ns->info.init_from_file.type == 0) {
    gk_neut_species_projection_init(app, ns, ns->info.projection, &ns->proj_init);
  }

  // Allocate objects for computing diagnostic moments.
  int ndm = ns->info.num_diag_moments;
  ns->moms = gkyl_malloc(sizeof(struct gk_species_moment[ndm]));
  for (int m=0; m<ndm; ++m)
    gk_neut_species_moment_init(app, ns, &ns->moms[m], ns->info.diag_moments[m], false);

  // Initialize boundary fluxes.
  ns->bflux = (struct gk_boundary_fluxes) { };
  // Additional bflux moments to step in time.
  struct gkyl_phase_diagnostics_inp add_bflux_moms_inp = (struct gkyl_phase_diagnostics_inp) { };
  // Set the operation type for the bflux app.
  enum gkyl_species_bflux_type bflux_type = GK_SPECIES_BFLUX_NONE;
  gk_neut_species_bflux_init(app, ns, &ns->bflux, bflux_type, add_bflux_moms_inp);

  // Initialize a Maxwellian/LTE (local thermodynamic equilibrium) projection routine
  // Projection routine optionally corrects all the Maxwellian/LTE moments
  // This routine is utilized by both reactions and BGK collisions
  ns->lte = (struct gk_lte) { };
  struct correct_all_moms_inp corr_inp = { .correct_all_moms = false,
    .max_iter = 0, .iter_eps = 10, .use_last_converged = false };
  gk_neut_species_lte_init(app, ns, &ns->lte, corr_inp);

  // Initialize the object that scales the species according to a balance
  // between recycling and reactions.
  ns->sca = (struct gk_scaling) { };
  gk_neut_species_scaling_init(app, ns, &ns->sca);

  // Initialize BGK collisions with null type (not applicable to fluids).
  ns->bgk = (struct gk_bgk_collisions) { };
  ns->info.collisions.collision_id = 0;
  gk_neut_species_bgk_init(app, ns, &ns->bgk);

  // Initialize positivity enforcing operator with null type (NYI for fluids).
  ns->positivity = (struct gk_positivity) { };
  gk_neut_species_positivity_init(app, ns, &ns->positivity);

  // Initialize reactions with charged species (NYI for fluids).
  ns->react_neut = (struct gk_react) { };
  ns->info.react_neut.num_react = 0;
  gk_neut_species_react_init(app, ns, ns->info.react_neut, &ns->react_neut);

  ns->src = (struct gk_source) { };
  if (!ns->info.is_static) {
    gk_neut_species_fluid_init_dynamic(gk, app, ns);
  }
  else {
    gk_neut_species_fluid_init_static(gk, app, ns);
  }
}

