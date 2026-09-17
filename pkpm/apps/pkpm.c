#include <stdarg.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_basis.h>
#include <gkyl_comm_io.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_null_comm.h>
#include <gkyl_pkpm_priv.h>
#include <gkyl_app_priv.h>

#include <mpack.h>

// returned gkyl_array_meta must be freed using pkpm_array_meta_release
static struct gkyl_msgpack_data*
pkpm_array_meta_new(struct pkpm_output_meta meta)
{
  struct gkyl_msgpack_data *mt = gkyl_malloc(sizeof(*mt));

  mt->meta_sz = 0;
  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mt->meta, &mt->meta_sz);

  // add some data to mpack
  mpack_build_map(&writer);
  
  mpack_write_cstr(&writer, "time");
  mpack_write_double(&writer, meta.stime);

  mpack_write_cstr(&writer, "frame");
  mpack_write_i64(&writer, meta.frame);

  mpack_write_cstr(&writer, "polyOrder");
  mpack_write_i64(&writer, meta.poly_order);

  mpack_write_cstr(&writer, "basisType");
  mpack_write_cstr(&writer, meta.basis_type);

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    free(mt->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mt);
    mt = 0;
  }

  return mt;
}

static void
pkpm_array_meta_release(struct gkyl_msgpack_data *mt)
{
  if (!mt) return;
  MPACK_FREE(mt->meta);
  gkyl_free(mt);
}

static struct pkpm_output_meta
pkpm_meta_from_mpack(struct gkyl_msgpack_data *mt)
{
  struct pkpm_output_meta meta = { .frame = 0, .stime = 0.0 };

  if (mt->meta_sz > 0) {
    mpack_tree_t tree;
    mpack_tree_init_data(&tree, mt->meta, mt->meta_sz);
    mpack_tree_parse(&tree);
    mpack_node_t root = mpack_tree_root(&tree);

    mpack_node_t tm_node = mpack_node_map_cstr(root, "time");
    meta.stime = mpack_node_double(tm_node);

    mpack_node_t fr_node = mpack_node_map_cstr(root, "frame");
    meta.frame = mpack_node_i64(fr_node);

    mpack_node_t po_node = mpack_node_map_cstr(root, "polyOrder");
    meta.poly_order = mpack_node_i64(po_node);

    mpack_node_t bt_node = mpack_node_map_cstr(root, "basisType");
    char *basis_type = mpack_node_cstr_alloc(bt_node, 64);
    strcpy(meta.basis_type_nm, basis_type);
    meta.basis_type = meta.basis_type_nm;
    MPACK_FREE(basis_type);

    mpack_tree_destroy(&tree);
  }
  return meta;
}

gkyl_pkpm_app*
gkyl_pkpm_app_new(struct gkyl_pkpm *pkpm)
{
  disable_denorm_float();

  assert(pkpm->num_species <= GKYL_MAX_SPECIES);

  gkyl_pkpm_app *app = gkyl_malloc(sizeof(gkyl_pkpm_app));

  int cdim = app->cdim = pkpm->cdim;
  int vdim = app->vdim = pkpm->vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = pkpm->poly_order;
  int ns = app->num_species = pkpm->num_species;

  double cfl_frac = pkpm->cfl_frac == 0 ? 1.0 : pkpm->cfl_frac;
  app->cfl = cfl_frac;

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = pkpm->parallelism.use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif

  app->num_periodic_dir = pkpm->num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = pkpm->periodic_dirs[d];

  strcpy(app->name, pkpm->name);
  app->tcurr = 0.0; // reset on init

  if (app->use_gpu) {
    // allocate device basis if we are using GPUs
    app->basis_on_dev.basis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
    app->basis_on_dev.confBasis = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  }
  else {
    app->basis_on_dev.basis = &app->basis;
    app->basis_on_dev.confBasis = &app->confBasis;
  }

  // basis functions
  switch (pkpm->basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&app->confBasis, cdim, poly_order);
      if (poly_order > 1) {
        gkyl_cart_modal_serendip(&app->basis, pdim, poly_order);
        if (vdim > 0)
          gkyl_cart_modal_serendip(&app->velBasis, vdim, poly_order);
      } else if (poly_order == 1) {
        /* Force hybrid basis (p=2 in velocity space). */
        gkyl_cart_modal_hybrid(&app->basis, cdim, vdim);
        if (vdim > 0)
          gkyl_cart_modal_serendip(&app->velBasis, vdim, 2);
      }

      if (app->use_gpu) {
        gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.confBasis, cdim, poly_order);
        if (poly_order > 1) {
          gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev.basis, pdim, poly_order);
        } else if (poly_order == 1) {
          /* Force hybrid basis (p=2 in velocity space). */
          gkyl_cart_modal_hybrid_cu_dev(app->basis_on_dev.basis, cdim, vdim);
        }
      }
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      gkyl_cart_modal_tensor(&app->basis, pdim, poly_order);
      gkyl_cart_modal_tensor(&app->confBasis, cdim, poly_order);
      if (vdim > 0)
        gkyl_cart_modal_tensor(&app->velBasis, vdim, poly_order);
      break;

    default:
      assert(false);
      break;
  }

  gkyl_rect_grid_init(&app->grid, cdim, pkpm->lower, pkpm->upper, pkpm->cells);

  int ghost[] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&app->grid, ghost, &app->global_ext, &app->global);

  if (pkpm->parallelism.comm == 0) {
    int cuts[3] = { 1, 1, 1 };
    app->decomp = gkyl_rect_decomp_new_from_cuts(cdim, cuts, &app->global);
    
    app->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = app->decomp,
        .use_gpu = app->use_gpu
      }
    );
    
    // Global and local ranges are same, and so just copy them.
    memcpy(&app->local, &app->global, sizeof(struct gkyl_range));
    memcpy(&app->local_ext, &app->global_ext, sizeof(struct gkyl_range));
  }
  else {
    // Create decomp.
    app->decomp = gkyl_rect_decomp_new_from_cuts(app->cdim, pkpm->parallelism.cuts, &app->global);

    // Create a new communicator with the decomposition in it.
    app->comm = gkyl_comm_split_comm(pkpm->parallelism.comm, 0, app->decomp);

    // Create local and local_ext.
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    gkyl_create_ranges(&app->decomp->ranges[rank], ghost, &app->local_ext, &app->local);
  }

  // local skin and ghost ranges for configuration space fields
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&app->lower_skin[dir], &app->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->local_ext, ghost); 
    gkyl_skin_ghost_ranges(&app->upper_skin[dir], &app->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->local_ext, ghost);
  }

  // Configuration space geometry initialization
  // Note: *only* uses a p=1 DG representation of the geometry (JJ: 05/03/24)
  app->c2p_ctx = app->mapc2p = 0;  
  app->has_mapc2p = pkpm->mapc2p ? true : false;

  if (app->has_mapc2p) {
    // initialize computational to physical space mapping
    app->c2p_ctx = pkpm->c2p_ctx;
    app->mapc2p = pkpm->mapc2p;

    // we project mapc2p on p=1 basis functions
    struct gkyl_basis basis;
    gkyl_cart_modal_tensor(&basis, cdim, 1);

    // initialize DG field representing mapping
    struct gkyl_array *c2p = mkarr(false, cdim*basis.num_basis, app->local_ext.volume);
    gkyl_eval_on_nodes *ev_c2p = gkyl_eval_on_nodes_new(&app->grid, &basis, cdim, pkpm->mapc2p, pkpm->c2p_ctx);
    gkyl_eval_on_nodes_advance(ev_c2p, 0.0, &app->local_ext, c2p);

    // write DG projection of mapc2p to file
    cstr fileNm = cstr_from_fmt("%s-mapc2p.gkyl", app->name);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, 0, c2p, fileNm.str);
    cstr_drop(&fileNm);

    gkyl_array_release(c2p);
    gkyl_eval_on_nodes_release(ev_c2p);
  }

  // create geometry object
  app->geom = gkyl_wave_geom_new(&app->grid, &app->local_ext,
    app->mapc2p, app->c2p_ctx, app->use_gpu);

  // PKPM system *always* has a field to define the local magnetic field direction
  app->field = pkpm_field_new(pkpm, app);

  // allocate space to store species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct pkpm_species[ns])) : 0;

  // set info for each species: this needs to be done here as we need
  // to access species name from pkpm_species_init
  for (int i=0; i<ns; ++i)
    app->species[i].info = pkpm->species[i];

  // initialize each species
  for (int i=0; i<ns; ++i) 
    pkpm_species_init(pkpm, app, &app->species[i]);

  // initialize each species cross-species terms: this has to be done here
  // as need pointers to colliding species' collision objects
  // allocated in the previous step
  for (int i=0; i<ns; ++i) {
    if (app->species[i].collision_id == GKYL_LBO_COLLISIONS
      && app->species[i].lbo.num_cross_collisions) {
      pkpm_species_lbo_cross_init(app, &app->species[i], &app->species[i].lbo);
    }
  }

  // Set the appropriate update function for taking a single time step
  app->use_explicit_source = pkpm->use_explicit_source;
  if (app->use_explicit_source) {
    // If momentum-EM field coupling is explicit, 
    // we use a pure explicit SSP RK3 method.
    app->update_func = pkpm_update_explicit_ssp_rk3;
  }
  else {
    // By default: we perform a first-order operator split 
    // to implicitly treat the momentum-EM field coupling. 
    app->pkpm_em = pkpm_fluid_em_coupling_init(app);
    app->update_func = pkpm_update_op_split;
  }
  
  // initialize stat object
  app->stat = (struct gkyl_pkpm_stat) {
    .use_gpu = app->use_gpu,
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };

  return app;
}

struct pkpm_species *
pkpm_find_species(const gkyl_pkpm_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return &app->species[i];
  return 0;
}

int
pkpm_find_species_idx(const gkyl_pkpm_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return i;
  return -1;
}

void
gkyl_pkpm_app_apply_ic(gkyl_pkpm_app* app, double t0)
{
  app->tcurr = t0;
  gkyl_pkpm_app_apply_ic_field(app, t0);
  for (int i=0; i<app->num_species; ++i) {
    gkyl_pkpm_app_apply_ic_species(app, i, t0);
  }
}

void
gkyl_pkpm_app_apply_ic_field(gkyl_pkpm_app* app, double t0)
{
  app->tcurr = t0;

  struct timespec wtm = gkyl_wall_clock();
  pkpm_field_apply_ic(app, app->field, t0);
  app->stat.init_field_tm += gkyl_time_diff_now_sec(wtm);

  pkpm_field_apply_bc(app, app->field, app->field->em);
  // bvar is computed over extended range, so calculate it after
  // we apply BCs in the initialization step. Note that the apply_bc call
  // may not apply BCs in the corner cells and the corner values will just
  // be what comes from initializing the EM field over the extended range
  pkpm_field_calc_bvar(app, app->field, app->field->em); 
}

void
gkyl_pkpm_app_apply_ic_species(gkyl_pkpm_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  pkpm_species_apply_ic(app, &app->species[sidx], t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);

  pkpm_species_apply_bc(app, &app->species[sidx], app->species[sidx].f);
  pkpm_fluid_species_apply_bc(app, &app->species[sidx], app->species[sidx].fluid);
}

void
gkyl_pkpm_app_calc_integrated_mom(gkyl_pkpm_app* app, double tm)
{
  double avals[9], avals_global[9];

  struct timespec wst = gkyl_wall_clock();

  for (int i=0; i<app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    gkyl_array_clear(s->integ_pkpm_mom, 0.0);

    pkpm_species_calc_pkpm_vars(app, s, s->f, s->fluid);
    gkyl_dg_calc_pkpm_integrated_vars(s->calc_pkpm_vars, &app->local, 
      s->pkpm_moms.marr, s->fluid, 
      s->pkpm_prim, s->integ_pkpm_mom);
    gkyl_array_scale_range(s->integ_pkpm_mom, app->grid.cellVolume, &(app->local));
    if (app->use_gpu) {
      gkyl_array_reduce_range(s->red_integ_diag, s->integ_pkpm_mom, GKYL_SUM, &(app->local));
      gkyl_cu_memcpy(avals, s->red_integ_diag, sizeof(double[9]), GKYL_CU_MEMCPY_D2H);
    }
    else { 
      gkyl_array_reduce_range(avals, s->integ_pkpm_mom, GKYL_SUM, &(app->local));
    }

    gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_SUM, 9, avals, avals_global);
    gkyl_dynvec_append(s->integ_diag, tm, avals_global);
  }

  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += app->num_species;
}

void
gkyl_pkpm_app_calc_integrated_L2_f(gkyl_pkpm_app* app, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  for (int i=0; i<app->num_species; ++i) {
    struct pkpm_species *s = &app->species[i];
    pkpm_species_calc_L2(app, tm, s);
  }
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += app->num_species;
}

void
gkyl_pkpm_app_calc_field_energy(gkyl_pkpm_app* app, double tm)
{
  struct timespec wst = gkyl_wall_clock();
  pkpm_field_calc_energy(app, tm, app->field);
  app->stat.diag_tm += gkyl_time_diff_now_sec(wst);
  app->stat.n_diag += 1;
}

void
gkyl_pkpm_app_write(gkyl_pkpm_app* app, double tm, int frame)
{
  app->stat.n_io += 1;
  struct timespec wtm = gkyl_wall_clock();
  
  gkyl_pkpm_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i) {
    gkyl_pkpm_app_write_species(app, i, tm, frame);
    gkyl_pkpm_app_write_mom(app, i, tm, frame);
  }

  app->stat.io_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_pkpm_app_write_field(gkyl_pkpm_app* app, double tm, int frame)
{
  struct gkyl_msgpack_data *mt = pkpm_array_meta_new( (struct pkpm_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  const char *fmt = "%s-field_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(app->field->em_host, app->field->em);
  }
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
    mt, app->field->em_host, fileNm);

  if (app->field->has_ext_em) {
    // Only write out external fields at t=0 or if they are time-dependent
    if (frame == 0 || app->field->ext_em_evolve) {
      const char *fmt_ext_em = "%s-field_ext_em_%d.gkyl";
      int sz_ext_em = gkyl_calc_strlen(fmt_ext_em, app->name, frame);
      char fileNm_ext_em[sz_ext_em+1]; // ensures no buffer overflow
      snprintf(fileNm_ext_em, sizeof fileNm_ext_em, fmt_ext_em, app->name, frame);

      // External EM field computed with project on basis, so just use host copy 
      pkpm_field_calc_ext_em(app, app->field, tm);
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
        mt, app->field->ext_em_host, fileNm_ext_em);
    }
  }

  if (app->field->has_app_current) {
    // Only write out applied currents at t=0 or if they are time-dependent
    if (frame == 0 || app->field->app_current_evolve) {
      const char *fmt_app_current = "%s-field_app_current_%d.gkyl";
      int sz_app_current = gkyl_calc_strlen(fmt_app_current, app->name, frame);
      char fileNm_app_current[sz_app_current+1]; // ensures no buffer overflow
      snprintf(fileNm_app_current, sizeof fileNm_app_current, fmt_app_current, app->name, frame);

      // Applied currents computed with project on basis, so just use host copy 
      pkpm_field_calc_app_current(app, app->field, tm);
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
        mt, app->field->app_current_host, fileNm_app_current);
    }
  }  

  pkpm_array_meta_release(mt);    
}

void
gkyl_pkpm_app_write_species(gkyl_pkpm_app* app, int sidx, double tm, int frame)
{
  struct gkyl_msgpack_data *mt = pkpm_array_meta_new( (struct pkpm_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->basis.id
    }
  );

  struct pkpm_species *s = &app->species[sidx];

  const char *fmt = "%s-%s_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->info.name, frame);

  if (app->use_gpu) {
    // copy data from device to host before writing it out
    gkyl_array_copy(s->f_host, s->f);
  }
  gkyl_comm_array_write(s->comm, &s->grid, &s->local,
    mt, s->f_host, fileNm);  

  pkpm_array_meta_release(mt);    
}

void
gkyl_pkpm_app_write_mom(gkyl_pkpm_app* app, int sidx, double tm, int frame)
{
  struct gkyl_msgpack_data *mt = pkpm_array_meta_new( (struct pkpm_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  struct pkpm_species *s = &app->species[sidx];

  // Construct the file handles for the three quantities (PKPM moments, PKPM fluid variables, PKPM update variables)
  const char *fmt = "%s-%s_pkpm_moms_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->info.name, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->info.name, frame);

  const char *fmt_fluid = "%s-%s_pkpm_fluid_%d.gkyl";
  int sz_fluid = gkyl_calc_strlen(fmt_fluid, app->name, s->info.name, frame);
  char fileNm_fluid[sz_fluid+1]; // ensures no buffer overflow
  snprintf(fileNm_fluid, sizeof fileNm_fluid, fmt_fluid, app->name, s->info.name, frame);

  const char *fmt_pkpm_vars = "%s-%s_pkpm_vars_%d.gkyl";
  int sz_pkpm_vars = gkyl_calc_strlen(fmt_pkpm_vars, app->name, s->info.name, frame);
  char fileNm_pkpm_vars[sz_pkpm_vars+1]; // ensures no buffer overflow
  snprintf(fileNm_pkpm_vars, sizeof fileNm_pkpm_vars, fmt_pkpm_vars, app->name, s->info.name, frame);

  // Compute the PKPM variables including moments and primitive variables 
  // and construct arrays for writing out fluid and other pkpm variables.
  pkpm_species_moment_calc(&s->pkpm_moms_diag, s->local, app->local, s->f);
  pkpm_species_calc_pkpm_vars(app, s, s->f, s->fluid);
  pkpm_species_calc_pkpm_update_vars(app, s, s->f); 
  gkyl_dg_calc_pkpm_vars_io(s->calc_pkpm_vars, &app->local, 
    s->pkpm_moms.marr, s->fluid, 
    s->pkpm_p_ij, s->pkpm_prim, 
    s->pkpm_accel, s->fluid_io, s->pkpm_vars_io);

  // copy data from device to host before writing it out
  if (app->use_gpu) {
    gkyl_array_copy(s->pkpm_moms_diag.marr_host, s->pkpm_moms_diag.marr);
    gkyl_array_copy(s->fluid_io_host, s->fluid_io);
    gkyl_array_copy(s->pkpm_vars_io_host, s->pkpm_vars_io);
  }

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
    mt, s->pkpm_moms_diag.marr_host, fileNm);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
    mt, s->fluid_io_host, fileNm_fluid);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, 
    mt, s->pkpm_vars_io_host, fileNm_pkpm_vars);

  pkpm_array_meta_release(mt);     
}

void
gkyl_pkpm_app_write_integrated_mom(gkyl_pkpm_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated diagnostic moments
      const char *fmt = "%s-%s-%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, app->species[i].info.name,
        "imom");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[i].info.name,
        "imom");

      if (app->species[i].is_first_integ_write_call) {
        gkyl_dynvec_write(app->species[i].integ_diag, fileNm);
        app->species[i].is_first_integ_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(app->species[i].integ_diag, fileNm);
      }
    }
    gkyl_dynvec_clear(app->species[i].integ_diag);
  }
}

void
gkyl_pkpm_app_write_integrated_L2_f(gkyl_pkpm_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out integrated L^2
      const char *fmt = "%s-%s-%s.gkyl";
      int sz = gkyl_calc_strlen(fmt, app->name, app->species[i].info.name,
        "L2");
      char fileNm[sz+1]; // ensures no buffer overflow
      snprintf(fileNm, sizeof fileNm, fmt, app->name, app->species[i].info.name,
        "L2");

      if (app->species[i].is_first_integ_L2_write_call) {
        // write to a new file (this ensure previous output is removed)
        gkyl_dynvec_write(app->species[i].integ_L2_f, fileNm);
        app->species[i].is_first_integ_L2_write_call = false;
      }
      else {
        // append to existing file
        gkyl_dynvec_awrite(app->species[i].integ_L2_f, fileNm);
      }
    }
    gkyl_dynvec_clear(app->species[i].integ_L2_f);
  }
}

void
gkyl_pkpm_app_write_field_energy(gkyl_pkpm_app* app)
{
  // write out diagnostic moments
  const char *fmt = "%s-field-energy.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    if (app->field->is_first_energy_write_call) {
      // write to a new file (this ensure previous output is removed)
      gkyl_dynvec_write(app->field->integ_energy, fileNm);
      app->field->is_first_energy_write_call = false;
    }
    else {
      // append to existing file
      gkyl_dynvec_awrite(app->field->integ_energy, fileNm);
    }
  }
  gkyl_dynvec_clear(app->field->integ_energy);
}

void
gkyl_pkpm_app_train(gkyl_pkpm_app* app, double tm, int frame, struct gkyl_kann_net** ann, int num_input_moms, int* input_moms, int num_output_moms, int* output_moms,
  struct gkyl_kn_vec* input_data, struct gkyl_kn_vec* output_data)
{
  for (int i = 0; i < app->num_species; i++) {
    gkyl_pkpm_app_train_mom(app, i, tm, frame, ann, num_input_moms, input_moms, num_output_moms, output_moms, input_data, output_data);
  }
}

void
gkyl_pkpm_app_train_mom(gkyl_pkpm_app* app, int sidx, double tm, int frame, struct gkyl_kann_net** ann, int num_input_moms, int* input_moms, int num_output_moms, int* output_moms,
  struct gkyl_kn_vec* input_data, struct gkyl_kn_vec* output_data)
{
  struct pkpm_species *s = &app->species[sidx];
  pkpm_species_moment_calc(&s->pkpm_moms_diag, s->local, app->local, s->f);

  if (app->use_gpu) {
    gkyl_array_copy(s->pkpm_moms_diag.marr_host, s->pkpm_moms_diag.marr);
  }

  int cell_count = 0;
  if (app->cdim == 1) {
    cell_count = app->grid.cells[0];
  }
  else if (app->cdim == 2) {
    cell_count = app->grid.cells[0] * app->grid.cells[1];
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local, iter.idx);

    const double *pkpm_moms_diag_d = gkyl_array_cfetch(s->pkpm_moms_diag.marr_host, loc);

    for (int i = 0; i < num_input_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          input_data->vals[count][i * 2] = (float)pkpm_moms_diag_d[input_moms[i] * 2];
          input_data->vals[count][(i * 2) + 1] = (float)pkpm_moms_diag_d[(input_moms[i] * 2) + 1];
        }
        else if (app->cdim == 2) {
          input_data->vals[count][i * 4] = (float)pkpm_moms_diag_d[input_moms[i] * 4];
          input_data->vals[count][(i * 4) + 1] = (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 1];
          input_data->vals[count][(i * 4) + 2] = (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 2];
          input_data->vals[count][(i * 4) + 3] = (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 3];
        }
      }
      else if (app->poly_order == 2) {
        input_data->vals[count][i * 3] = (float)pkpm_moms_diag_d[input_moms[i] * 3];
        input_data->vals[count][(i * 3) + 1] = (float)pkpm_moms_diag_d[(input_moms[i] * 3) + 1];
        input_data->vals[count][(i * 3) + 2] = (float)pkpm_moms_diag_d[(input_moms[i] * 3) + 2];
      }
    }

    for (int i = 0; i < num_output_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          output_data->vals[count][i * 2] = (float)pkpm_moms_diag_d[output_moms[i] * 2];
          output_data->vals[count][(i * 2) + 1] = (float)pkpm_moms_diag_d[(output_moms[i] * 2) + 1];
        }
        else if (app->cdim == 2) {
          output_data->vals[count][i * 4] = (float)pkpm_moms_diag_d[output_moms[i] * 4];
          output_data->vals[count][(i * 4) + 1] = (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 1];
          output_data->vals[count][(i * 4) + 2] = (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 2];
          output_data->vals[count][(i * 4) + 3] = (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 3];
        }
      }
      else if (app->poly_order == 2) {
        output_data->vals[count][i * 3] = (float)pkpm_moms_diag_d[output_moms[i] * 3];
        output_data->vals[count][(i * 3) + 1] = (float)pkpm_moms_diag_d[(output_moms[i] * 3) + 1];
        output_data->vals[count][(i * 3) + 2] = (float)pkpm_moms_diag_d[(output_moms[i] * 3) + 2];
      }
    }

    count += 1;
  }

  struct gkyl_kann_train_params params = {
    .learning_rate = 0.0001f,
    .mini_size = 64,
    .max_epoch = 50,
    .max_drop_streak = 10,
    .frac_val = 0.1f,
  };

  if (gkyl_kann_net_is_cu_dev(ann[sidx])) {
    // Mirror the host training data onto the device and train there.
    struct gkyl_kn_vec *input_data_cu = gkyl_kn_vec_cu_dev_new(input_data->nvec, input_data->N);
    struct gkyl_kn_vec *output_data_cu = gkyl_kn_vec_cu_dev_new(output_data->nvec, output_data->N);
    gkyl_kn_vec_copy(input_data_cu, input_data);
    gkyl_kn_vec_copy(output_data_cu, output_data);

    gkyl_kann_net_train_fnn1(ann[sidx], &params, input_data_cu, output_data_cu);

    gkyl_kn_vec_release(input_data_cu);
    gkyl_kn_vec_release(output_data_cu);
  }
  else {
    gkyl_kann_net_train_fnn1(ann[sidx], &params, input_data, output_data);
  }
}

void
gkyl_pkpm_app_write_nn(gkyl_pkpm_app* app, double tm, int frame, struct gkyl_kann_net** ann)
{
  for (int i = 0; i < app->num_species; i++) {
    gkyl_pkpm_app_write_nn_mom(app, i, tm, frame, ann);
  }
}

void
gkyl_pkpm_app_write_nn_mom(gkyl_pkpm_app* app, int sidx, double tm, int frame, struct gkyl_kann_net** ann)
{
  struct pkpm_species *s = &app->species[sidx];

  const char *fmt = "%s_moms_nn_%d-%s.dat";
  int sz = gkyl_calc_strlen(fmt, app->name, frame, s->info.name);
  char fileNm[sz + 1];
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame, s->info.name);

  gkyl_kann_net_save(ann[sidx], fileNm);
}

void
gkyl_pkpm_app_test(gkyl_pkpm_app* app, double tm, int frame, struct gkyl_kann_net** ann, int num_input_moms, int* input_moms, int num_output_moms, int* output_moms,
  struct gkyl_kn_vec* input_data_real, struct gkyl_kn_vec* output_data_real, struct gkyl_kn_vec* output_data_predicted)
{
  for (int i = 0; i < app->num_species; i++) {
    gkyl_pkpm_app_test_mom(app, i, tm, frame, ann, num_input_moms, input_moms, num_output_moms, output_moms,
      input_data_real, output_data_real, output_data_predicted);
  }
}

void
gkyl_pkpm_app_test_mom(gkyl_pkpm_app* app, int sidx, double tm, int frame, struct gkyl_kann_net** ann, int num_input_moms, int* input_moms, int num_output_moms, int* output_moms,
  struct gkyl_kn_vec* input_data_real, struct gkyl_kn_vec* output_data_real, struct gkyl_kn_vec* output_data_predicted)
{
  struct gkyl_msgpack_data *mt = pkpm_array_meta_new( (struct pkpm_output_meta) {
      .frame = frame,
      .stime = tm,
      .poly_order = app->poly_order,
      .basis_type = app->confBasis.id
    }
  );

  struct pkpm_species *s = &app->species[sidx];
  pkpm_species_moment_calc(&s->pkpm_moms_diag, s->local, app->local, s->f);

  if (app->use_gpu) {
    gkyl_array_copy(s->pkpm_moms_diag.marr_host, s->pkpm_moms_diag.marr);
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local, iter.idx);

    const double *pkpm_moms_diag_d = gkyl_array_cfetch(s->pkpm_moms_diag.marr_host, loc);

    for (int i = 0; i < num_input_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          input_data_real->vals[count][i * 2] = (float)pkpm_moms_diag_d[input_moms[i] * 2];
          input_data_real->vals[count][(i * 2) + 1] = (float)pkpm_moms_diag_d[(input_moms[i] * 2) + 1];
        }
        else if (app->cdim == 2) {
          input_data_real->vals[count][i * 4] = (float)pkpm_moms_diag_d[input_moms[i] * 4];
          input_data_real->vals[count][(i * 4) + 1] = (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 1];
          input_data_real->vals[count][(i * 4) + 2] = (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 2];
          input_data_real->vals[count][(i * 4) + 3] = (float)pkpm_moms_diag_d[(input_moms[i] * 4) + 3];
        }
      }
      else if (app->poly_order == 2) {
        input_data_real->vals[count][i * 3] = (float)pkpm_moms_diag_d[input_moms[i] * 3];
        input_data_real->vals[count][(i * 3) + 1] = (float)pkpm_moms_diag_d[(input_moms[i] * 3) + 1];
        input_data_real->vals[count][(i * 3) + 2] = (float)pkpm_moms_diag_d[(input_moms[i] * 3) + 2];
      }
    }

    for (int i = 0; i < num_output_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          output_data_real->vals[count][i * 2] = (float)pkpm_moms_diag_d[output_moms[i] * 2];
          output_data_real->vals[count][(i * 2) + 1] = (float)pkpm_moms_diag_d[(output_moms[i] * 2) + 1];
        }
        else if (app->cdim == 2) {
          output_data_real->vals[count][i * 4] = (float)pkpm_moms_diag_d[output_moms[i] * 4];
          output_data_real->vals[count][(i * 4) + 1] = (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 1];
          output_data_real->vals[count][(i * 4) + 2] = (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 2];
          output_data_real->vals[count][(i * 4) + 3] = (float)pkpm_moms_diag_d[(output_moms[i] * 4) + 3];
        }
      }
      else if (app->poly_order == 2) {
        output_data_real->vals[count][i * 3] = (float)pkpm_moms_diag_d[output_moms[i] * 3];
        output_data_real->vals[count][(i * 3) + 1] = (float)pkpm_moms_diag_d[(output_moms[i] * 3) + 1];
        output_data_real->vals[count][(i * 3) + 2] = (float)pkpm_moms_diag_d[(output_moms[i] * 3) + 2];
      }
    }

    count += 1;
  }

  // Batched inference over all sampled cells. On GPU, stage the input and
  // predicted output through device kn_vecs.
  if (gkyl_kann_net_is_cu_dev(ann[sidx])) {
    struct gkyl_kn_vec *input_data_cu = gkyl_kn_vec_cu_dev_new(input_data_real->nvec, input_data_real->N);
    struct gkyl_kn_vec *output_data_cu = gkyl_kn_vec_cu_dev_new(output_data_predicted->nvec, output_data_predicted->N);
    gkyl_kn_vec_copy(input_data_cu, input_data_real);

    gkyl_kann_net_apply(ann[sidx], input_data_cu, output_data_cu);

    gkyl_kn_vec_copy(output_data_predicted, output_data_cu);
    gkyl_kn_vec_release(input_data_cu);
    gkyl_kn_vec_release(output_data_cu);
  }
  else {
    gkyl_kann_net_apply(ann[sidx], input_data_real, output_data_predicted);
  }

  struct gkyl_range_iter iter_new;
  gkyl_range_iter_init(&iter_new, &app->local);
  long count_new = 0;
  while (gkyl_range_iter_next(&iter_new)) {
    long loc_new = gkyl_range_idx(&app->local, iter_new.idx);

    double *pkpm_moms_diag_d_new = gkyl_array_fetch(s->pkpm_moms_diag.marr_host, loc_new);

    for (int i = 0; i < num_output_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          pkpm_moms_diag_d_new[output_moms[i] * 2] = (double)output_data_predicted->vals[count_new][i * 2];
          pkpm_moms_diag_d_new[(output_moms[i] * 2) + 1] = (double)output_data_predicted->vals[count_new][(i * 2) + 1];
        }
        else if (app->cdim == 2) {
          pkpm_moms_diag_d_new[output_moms[i] * 4] = (double)output_data_predicted->vals[count_new][i * 4];
          pkpm_moms_diag_d_new[(output_moms[i] * 4) + 1] = (double)output_data_predicted->vals[count_new][(i * 4) + 1];
          pkpm_moms_diag_d_new[(output_moms[i] * 4) + 2] = (double)output_data_predicted->vals[count_new][(i * 4) + 2];
          pkpm_moms_diag_d_new[(output_moms[i] * 4) + 3] = (double)output_data_predicted->vals[count_new][(i * 4) + 3];
        }
      }
      else if (app->poly_order == 2) {
        pkpm_moms_diag_d_new[output_moms[i] * 3] = (double)output_data_predicted->vals[count_new][i * 3];
        pkpm_moms_diag_d_new[(output_moms[i] * 3) + 1] = (double)output_data_predicted->vals[count_new][(i * 3) + 1];
        pkpm_moms_diag_d_new[(output_moms[i] * 3) + 2] = (double)output_data_predicted->vals[count_new][(i * 3) + 2];
      }
    }

    count_new += 1;
  }

  const char *fmt = "%s-%s_pkpm_moms_predicted_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, s->info.name, frame);
  char fileNm[sz + 1];
  snprintf(fileNm, sizeof fileNm, fmt, app->name, s->info.name, frame);

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, s->pkpm_moms_diag.marr_host, fileNm);

  struct gkyl_range_iter iter_old;
  gkyl_range_iter_init(&iter_old, &app->local);
  long count_old = 0;
  while (gkyl_range_iter_next(&iter_old)) {
    long loc_old = gkyl_range_idx(&app->local, iter_old.idx);

    double *pkpm_moms_diag_d_old = gkyl_array_fetch(s->pkpm_moms_diag.marr_host, loc_old);

    for (int i = 0; i < num_output_moms; i++) {
      if (app->poly_order == 1) {
        if (app->cdim == 1) {
          pkpm_moms_diag_d_old[output_moms[i] * 2] = (double)output_data_real->vals[count_old][i * 2];
          pkpm_moms_diag_d_old[(output_moms[i] * 2) + 1] = (double)output_data_real->vals[count_old][(i * 2) + 1];
        }
        else if (app->cdim == 2) {
          pkpm_moms_diag_d_old[output_moms[i] * 4] = (double)output_data_real->vals[count_old][i * 4];
          pkpm_moms_diag_d_old[(output_moms[i] * 4) + 1] = (double)output_data_real->vals[count_old][(i * 4) + 1];
          pkpm_moms_diag_d_old[(output_moms[i] * 4) + 2] = (double)output_data_real->vals[count_old][(i * 4) + 2];
          pkpm_moms_diag_d_old[(output_moms[i] * 4) + 3] = (double)output_data_real->vals[count_old][(i * 4) + 3];
        }
      }
      else if (app->poly_order == 2) {
        pkpm_moms_diag_d_old[output_moms[i] * 3] = (double)output_data_real->vals[count_old][i * 3];
        pkpm_moms_diag_d_old[(output_moms[i] * 3) + 1] = (double)output_data_real->vals[count_old][(i * 3) + 1];
        pkpm_moms_diag_d_old[(output_moms[i] * 3) + 2] = (double)output_data_real->vals[count_old][(i * 3) + 2];
      }
    }

    count_old += 1;
  }
}

struct gkyl_update_status
gkyl_pkpm_update(gkyl_pkpm_app* app, double dt)
{
  app->stat.nup += 1;

  struct timespec wst = gkyl_wall_clock();
  struct gkyl_update_status status = app->update_func(app, dt);
  app->tcurr += status.dt_actual;

  app->stat.total_tm += gkyl_time_diff_now_sec(wst);

  // Check for any CUDA errors during time step
  if (app->use_gpu)
    checkCuda(cudaGetLastError());
  return status;
}

struct gkyl_pkpm_stat
gkyl_pkpm_app_stat(gkyl_pkpm_app* app)
{
  pkpm_species_tm(app);
  pkpm_species_coll_tm(app);
  return app->stat;
}

static void
range_stat_write(gkyl_pkpm_app* app, const char *nm, const struct gkyl_range *r, FILE *fp)
{
  gkyl_pkpm_app_cout(app, fp, " %s_cells : [ ", nm);
  for (int i=0; i<r->ndim; ++i)
    gkyl_pkpm_app_cout(app, fp, " %d, ", gkyl_range_shape(r, i));
  gkyl_pkpm_app_cout(app, fp, " ],\n");
}

// ensure stats across processors are made consistent
static void
comm_reduce_app_stat(const gkyl_pkpm_app* app,
  const struct gkyl_pkpm_stat *local, struct gkyl_pkpm_stat *global)
{
  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz == 1) {
    memcpy(global, local, sizeof(struct gkyl_pkpm_stat));
    return;
  }

  global->use_gpu = local->use_gpu;

  enum { NUP, NFEULER, NSTAGE_2_FAIL, NSTAGE_3_FAIL, L_END };
  int64_t l_red[] = {
    [NUP] = local->nup,
    [NFEULER] = local->nfeuler,
    [NSTAGE_2_FAIL] = local->nstage_2_fail,
    [NSTAGE_3_FAIL] = local->nstage_3_fail
  };

  int64_t l_red_global[L_END];
  gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, L_END, l_red, l_red_global);

  global->nup = l_red_global[NUP];
  global->nfeuler = l_red_global[NFEULER];
  global->nstage_2_fail = l_red_global[NSTAGE_2_FAIL];
  global->nstage_3_fail = l_red_global[NSTAGE_3_FAIL];  

  enum {
    TOTAL_TM, RK3_TM, PKPM_EM_TM, INIT_SPECIES_TM, INIT_FLUID_SPECIES_TM, INIT_FIELD_TM,
    SPECIES_RHS_TM, FLUID_SPECIES_RHS_TM, SPECIES_COLL_MOM_TM,
    SPECIES_COL_TM, SPECIES_PKPM_VARS_TM, FIELD_RHS_TM, FIELD_EM_VARS_TM, CURRENT_TM,
    SPECIES_OMEGA_CFL_TM, FIELD_OMEGA_CFL_TM, DIAG_TM, IO_TM,
    SPECIES_BC_TM, FLUID_SPECIES_BC_TM, FIELD_BC_TM,
    D_END
  };

  double d_red[D_END] = {
    [TOTAL_TM] = local->total_tm,
    [RK3_TM] = local->rk3_tm,
    [PKPM_EM_TM] = local->pkpm_em_tm,
    [INIT_SPECIES_TM] = local->init_species_tm,
    [INIT_FLUID_SPECIES_TM] = local->init_fluid_species_tm,
    [INIT_FIELD_TM] = local->field_rhs_tm,
    [SPECIES_RHS_TM] = local->species_rhs_tm,
    [FLUID_SPECIES_RHS_TM] = local->fluid_species_rhs_tm,
    [SPECIES_COLL_MOM_TM] = local->species_coll_mom_tm,
    [SPECIES_COL_TM] = local->species_coll_tm,
    [SPECIES_PKPM_VARS_TM] = local->species_pkpm_vars_tm,
    [FIELD_RHS_TM] = local->field_rhs_tm,
    [FIELD_EM_VARS_TM] = local->field_em_vars_tm,
    [CURRENT_TM] = local->current_tm,
    [SPECIES_OMEGA_CFL_TM] = local->species_omega_cfl_tm,
    [FIELD_OMEGA_CFL_TM] = local->field_omega_cfl_tm,
    [DIAG_TM] = local->diag_tm,
    [IO_TM] = local->io_tm,
    [SPECIES_BC_TM] = local->species_bc_tm,
    [FLUID_SPECIES_BC_TM] = local->fluid_species_bc_tm,
    [FIELD_BC_TM] = local->field_bc_tm
  };

  double d_red_global[D_END];
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, D_END, d_red, d_red_global);
  
  global->total_tm = d_red_global[TOTAL_TM];
  global->rk3_tm = d_red_global[RK3_TM];
  global->pkpm_em_tm = d_red_global[PKPM_EM_TM];
  global->init_species_tm = d_red_global[INIT_SPECIES_TM];
  global->init_fluid_species_tm = d_red_global[INIT_FLUID_SPECIES_TM];
  global->field_rhs_tm = d_red_global[INIT_FIELD_TM];
  global->species_rhs_tm = d_red_global[SPECIES_RHS_TM];
  global->fluid_species_rhs_tm = d_red_global[FLUID_SPECIES_RHS_TM];
  global->species_coll_mom_tm = d_red_global[SPECIES_COLL_MOM_TM];
  global->species_coll_tm = d_red_global[SPECIES_COL_TM];
  global->species_pkpm_vars_tm = d_red_global[SPECIES_PKPM_VARS_TM];
  global->field_rhs_tm = d_red_global[FIELD_RHS_TM];
  global->field_em_vars_tm = d_red_global[FIELD_EM_VARS_TM];
  global->current_tm = d_red_global[CURRENT_TM];
  global->species_omega_cfl_tm = d_red_global[SPECIES_OMEGA_CFL_TM];
  global->field_omega_cfl_tm = d_red_global[FIELD_OMEGA_CFL_TM];
  global->diag_tm = d_red_global[DIAG_TM];
  global->io_tm = d_red_global[IO_TM];
  global->species_bc_tm = d_red_global[SPECIES_BC_TM];
  global->fluid_species_bc_tm = d_red_global[FLUID_SPECIES_BC_TM];
  global->field_bc_tm = d_red_global[FIELD_BC_TM];

  // misc data needing reduction

  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_2_dt_diff,
    global->stage_2_dt_diff);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_3_dt_diff,
    global->stage_3_dt_diff);

  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_drag_tm,
    global->species_lbo_coll_drag_tm);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, GKYL_MAX_SPECIES, local->species_lbo_coll_diff_tm,
    global->species_lbo_coll_diff_tm);
}

void
gkyl_pkpm_app_stat_write(gkyl_pkpm_app* app)
{
  const char *fmt = "%s-%s";
  int sz = gkyl_calc_strlen(fmt, app->name, "stat.json");
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, "stat.json");

  int num_ranks;
  gkyl_comm_get_size(app->comm, &num_ranks);

  char buff[70];
  time_t t = time(NULL);
  struct tm curr_tm = *localtime(&t);

  pkpm_species_coll_tm(app);
  pkpm_species_tm(app);

  struct gkyl_pkpm_stat stat = { };
  comm_reduce_app_stat(app, &app->stat, &stat);
  
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  if (rank == 0) fp = fopen(fileNm, "a");

  gkyl_pkpm_app_cout(app, fp, "{\n");

  if (strftime(buff, sizeof buff, "%c", &curr_tm))
    gkyl_pkpm_app_cout(app, fp, " date : %s,\n", buff);

  gkyl_pkpm_app_cout(app, fp, " use_gpu : %d,\n", stat.use_gpu);
  gkyl_pkpm_app_cout(app, fp, " num_ranks : %d,\n", num_ranks); 
  
  for (int s=0; s<app->num_species; ++s)
    range_stat_write(app, app->species[s].info.name, &app->species[s].global, fp);
  
  gkyl_pkpm_app_cout(app, fp, " nup : %ld,\n", stat.nup);
  gkyl_pkpm_app_cout(app, fp, " nfeuler : %ld,\n", stat.nfeuler);
  gkyl_pkpm_app_cout(app, fp, " nstage_2_fail : %ld,\n", stat.nstage_2_fail);
  gkyl_pkpm_app_cout(app, fp, " nstage_3_fail : %ld,\n", stat.nstage_3_fail);

  gkyl_pkpm_app_cout(app, fp, " stage_2_dt_diff : [ %lg, %lg ],\n",
    stat.stage_2_dt_diff[0], stat.stage_2_dt_diff[1]);
  gkyl_pkpm_app_cout(app, fp, " stage_3_dt_diff : [ %lg, %lg ],\n",
    stat.stage_3_dt_diff[0], stat.stage_3_dt_diff[1]);

  gkyl_pkpm_app_cout(app, fp, " total_tm : %lg,\n", stat.total_tm);
  gkyl_pkpm_app_cout(app, fp, " rk3_tm : %lg,\n", stat.rk3_tm);
  if (!app->use_explicit_source) {
    gkyl_pkpm_app_cout(app, fp, " fluid_em_coupling_tm : %lg,\n", stat.pkpm_em_tm);
  }
  gkyl_pkpm_app_cout(app, fp, " init_species_tm : %lg,\n", stat.init_species_tm);
  gkyl_pkpm_app_cout(app, fp, " init_field_tm : %lg,\n", stat.init_field_tm);
  
  gkyl_pkpm_app_cout(app, fp, " species_rhs_tm : %lg,\n", stat.species_rhs_tm);
  gkyl_pkpm_app_cout(app, fp, " species_bc_tm : %lg,\n", stat.species_bc_tm);
  gkyl_pkpm_app_cout(app, fp, " species_pkpm_vars_tm : %lg,\n", stat.species_pkpm_vars_tm);

  gkyl_pkpm_app_cout(app, fp, " species_coll_mom_tm : %lg,\n", stat.species_coll_mom_tm);
  gkyl_pkpm_app_cout(app, fp, " species_coll_tm : %lg,\n", stat.species_coll_tm);
  for (int s=0; s<app->num_species; ++s) {
    gkyl_pkpm_app_cout(app, fp, " species_coll_drag_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_drag_tm[s]);
    gkyl_pkpm_app_cout(app, fp, " species_coll_diff_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_diff_tm[s]);
  }
  
  gkyl_pkpm_app_cout(app, fp, " fluid_species_rhs_tm : %lg,\n", stat.fluid_species_rhs_tm);
  gkyl_pkpm_app_cout(app, fp, " fluid_species_bc_tm : %lg,\n", stat.fluid_species_bc_tm);

  gkyl_pkpm_app_cout(app, fp, " field_rhs_tm : %lg,\n", stat.field_rhs_tm);
  gkyl_pkpm_app_cout(app, fp, " field_bc_tm : %lg,\n", stat.field_bc_tm);
  gkyl_pkpm_app_cout(app, fp, " field_em_vars_tm : %lg,\n", stat.field_em_vars_tm);
  if (app->use_explicit_source) {
    gkyl_pkpm_app_cout(app, fp, " current_tm : %lg,\n", stat.current_tm);
  }

  gkyl_pkpm_app_cout(app, fp, " ndiag : %ld,\n", stat.n_diag);
  gkyl_pkpm_app_cout(app, fp, " diag_tm : %lg\n", stat.diag_tm);
  
  gkyl_pkpm_app_cout(app, fp, " nspecies_omega_cfl : %ld,\n", stat.n_species_omega_cfl);
  gkyl_pkpm_app_cout(app, fp, " species_omega_cfl_tm : %lg\n", stat.species_omega_cfl_tm);

  gkyl_pkpm_app_cout(app, fp, " nfield_omega_cfl : %ld,\n", stat.n_field_omega_cfl);
  gkyl_pkpm_app_cout(app, fp, " field_omega_cfl_tm : %lg\n", stat.field_omega_cfl_tm);

  gkyl_pkpm_app_cout(app, fp, " nio : %ld,\n", stat.n_io);
  gkyl_pkpm_app_cout(app, fp, " io_tm : %lg\n", stat.io_tm);
  
  gkyl_pkpm_app_cout(app, fp, "}\n");

  if (rank == 0)
    fclose(fp);  

}

static struct gkyl_app_restart_status
header_from_file(gkyl_pkpm_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = { .io_status = 0 };
  
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (hdr.etype != GKYL_DOUBLE)
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
    }

    struct pkpm_output_meta meta =
      pkpm_meta_from_mpack( &(struct gkyl_msgpack_data) {
          .meta = hdr.meta,
          .meta_sz = hdr.meta_size
        }
      );

    rstat.frame = meta.frame;
    rstat.stime = meta.stime;

    gkyl_grid_sub_array_header_release(&hdr);
  }
  
  return rstat;
}

struct gkyl_app_restart_status
gkyl_pkpm_app_from_file_field(gkyl_pkpm_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->field->em_host, fname);
    if (app->use_gpu) {
      gkyl_array_copy(app->field->em, app->field->em_host);
    }
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      pkpm_field_apply_bc(app, app->field, app->field->em);
    }
  }

  // Compute external EM field and applied current if present
  // Computation necessary in case external EM field or applied current
  // are time-independent and not computed in the time-stepping loop
  // since they are not read-in as part of restarts. 
  pkpm_field_calc_ext_em(app, app->field, rstat.stime);
  pkpm_field_calc_app_current(app, app->field, rstat.stime);

  return rstat;
}

struct gkyl_app_restart_status 
gkyl_pkpm_app_from_file_species(gkyl_pkpm_app *app, int sidx,
  const char *fname)
{
  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  struct pkpm_species *s = &app->species[sidx];
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(s->comm, &s->grid, &s->local, s->f_host, fname);
    if (app->use_gpu) {
      gkyl_array_copy(s->f, s->f_host);
    }
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      pkpm_species_apply_bc(app, s, s->f);
    }
  }

  // Compute applied acceleration if present.
  // Computation necessary in case applied acceleration
  // is time-independent and not computed in the time-stepping loop
  // since it is not read-in as part of restarts. 
  pkpm_species_calc_app_accel(app, s, rstat.stime);

  return rstat;
}

struct gkyl_app_restart_status 
gkyl_pkpm_app_from_file_fluid_species(gkyl_pkpm_app *app, int sidx,
  const char *fname)
{
  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  struct pkpm_species *s = &app->species[sidx];
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    // Read in the full 10 component fluid array
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, s->fluid_io_host, fname);
    if (app->use_gpu) {
      gkyl_array_copy(s->fluid_io, s->fluid_io_host);
    }
    // Copy the relevant components of 10 component fluid array (components 1-4)
    gkyl_array_set_offset(s->fluid, 1.0, s->fluid_io, app->confBasis.num_basis);
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      pkpm_fluid_species_apply_bc(app, s, s->fluid);
    }
  }

  return rstat;
}

struct gkyl_app_restart_status
gkyl_pkpm_app_from_frame_field(gkyl_pkpm_app *app, int frame)
{
  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
  struct gkyl_app_restart_status rstat = gkyl_pkpm_app_from_file_field(app, fileNm.str);
  app->field->is_first_energy_write_call = false; // append to existing diagnostic
  cstr_drop(&fileNm);
  
  return rstat;
}

struct gkyl_app_restart_status
gkyl_pkpm_app_from_frame_species(gkyl_pkpm_app *app, int sidx, int frame)
{
  struct pkpm_species *s = &app->species[sidx];

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, s->info.name, frame);
  struct gkyl_app_restart_status rstat = gkyl_pkpm_app_from_file_species(app, sidx, fileNm.str);

  cstr fileNm_fluid = cstr_from_fmt("%s-%s_pkpm_fluid_%d.gkyl", app->name, s->info.name, frame);
  struct gkyl_app_restart_status rstat_fluid = gkyl_pkpm_app_from_file_fluid_species(app, sidx, fileNm_fluid.str);

  app->species[sidx].is_first_integ_write_call = false; // append to existing diagnostic
  app->species[sidx].is_first_integ_L2_write_call = false; // append to existing diagnostic
  cstr_drop(&fileNm);
  
  return rstat;
}

struct gkyl_app_restart_status
gkyl_pkpm_app_read_from_frame(gkyl_pkpm_app *app, int frame)
{
  struct gkyl_app_restart_status rstat;
  
  rstat = gkyl_pkpm_app_from_frame_field(app, frame);

  for (int i=0; i<app->num_species; i++) {
    rstat = gkyl_pkpm_app_from_frame_species(app, i, frame);
  }

  return rstat;
}

// private function to handle variable argument list for printing
static void
v_pkpm_app_cout(const gkyl_pkpm_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp) {
    vfprintf(fp, fmt, argp);
    fflush(fp);
  }
}

void
gkyl_pkpm_app_cout(const gkyl_pkpm_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_pkpm_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_pkpm_app_release(gkyl_pkpm_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    pkpm_species_release(app, &app->species[i]);
  }
  if (app->num_species > 0) {
    gkyl_free(app->species);
  }

  pkpm_field_release(app, app->field);
  if (!app->use_explicit_source) {
    pkpm_fluid_em_coupling_release(app, app->pkpm_em);
  }

  gkyl_comm_release(app->comm);
  gkyl_rect_decomp_release(app->decomp);
  
  gkyl_wave_geom_release(app->geom);

  if (app->use_gpu) {
    gkyl_cu_free(app->basis_on_dev.basis);
    gkyl_cu_free(app->basis_on_dev.confBasis);
  }

  gkyl_free(app);
}
