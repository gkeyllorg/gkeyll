#include <gkyl_array_rio_priv.h>
#include <gkyl_moment_priv.h>
#include <gkyl_null_comm.h>
#include <gkyl_util.h>

#include <mpack.h>

static inline int int_max(int a, int b) { return a > b ? a : b; }

struct gkyl_msgpack_data*
moment_array_meta_new(struct moment_output_meta meta)
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

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    free(mt->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mt);
    mt = 0;
  }

  return mt;
}

void
moment_array_meta_release(struct gkyl_msgpack_data *mt)
{
  if (!mt) return;
  MPACK_FREE(mt->meta);
  gkyl_free(mt);
}

struct moment_output_meta
moment_meta_from_mpack(struct gkyl_msgpack_data *mt)
{
  struct moment_output_meta meta = { .frame = 0, .stime = 0.0 };

  if (mt->meta_sz > 0) {
    mpack_tree_t tree;
    mpack_tree_init_data(&tree, mt->meta, mt->meta_sz);
    mpack_tree_parse(&tree);
    mpack_node_t root = mpack_tree_root(&tree);
    mpack_node_t tm_node = mpack_node_map_cstr(root, "time");
    meta.stime = mpack_node_double(tm_node);
    mpack_node_t fr_node = mpack_node_map_cstr(root, "frame");
    meta.frame = mpack_node_i64(fr_node);
    mpack_tree_destroy(&tree);
  }
  return meta;
}

gkyl_moment_app*
gkyl_moment_app_new(struct gkyl_moment *mom)
{
  disable_denorm_float();
  
  struct gkyl_moment_app *app = gkyl_calloc(1, sizeof(*app));

  int ndim = app->ndim = mom->ndim;
  strcpy(app->name, mom->name);
  app->tcurr = 0.0; // reset on init

  app->scheme_type = mom->scheme_type;
  
  app->mp_recon = mom->mp_recon;
  app->use_hybrid_flux_kep = mom->use_hybrid_flux_kep;
  
  if (app->scheme_type == GKYL_MOMENT_WAVE_PROP)
    app->update_func = moment_update_one_step;
  else if (app->scheme_type == GKYL_MOMENT_MP) 
    app->update_func = moment_update_ssp_rk3;
  else if (app->scheme_type == GKYL_MOMENT_KEP)
    app->update_func = moment_update_ssp_rk3;

  int ghost[3] = { 2, 2, 2 }; // 2 ghost-cells for wave
  if (mom->scheme_type != GKYL_MOMENT_WAVE_PROP)
    for (int d=0; d<3; ++d) ghost[d] = 3; // 3 for MP scheme and KEP

  for (int d=0; d<3; ++d) app->nghost[d] = ghost[d];
  
  gkyl_rect_grid_init(&app->grid, ndim, mom->lower, mom->upper, mom->cells);
  gkyl_create_grid_ranges(&app->grid, ghost, &app->global_ext, &app->global);

  if (mom->parallelism.comm == 0) {
    int cuts[3] = { 1, 1, 1 };
    app->decomp = gkyl_rect_decomp_new_from_cuts(app->ndim, cuts, &app->global);
    
    app->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .decomp = app->decomp,
        .sync_corners = true, // If no communicator, since some moment apps need corner syncs, turn on corner syncs. 
      }
    );
    
    // Global and local ranges are same, and so just copy them.
    memcpy(&app->local, &app->global, sizeof(struct gkyl_range));
    memcpy(&app->local_ext, &app->global_ext, sizeof(struct gkyl_range));
  }
  else {
    // Create decomp.
    app->decomp = gkyl_rect_decomp_new_from_cuts(app->ndim, mom->parallelism.cuts, &app->global);

    // Create a new communicator with the decomposition in it.
    app->comm = gkyl_comm_split_comm(mom->parallelism.comm, 0, app->decomp);

    // Create local and local_ext.
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    gkyl_create_ranges(&app->decomp->ranges[rank], ghost, &app->local_ext, &app->local);
  }

  skin_ghost_ranges_init(&app->skin_ghost, &app->local_ext, ghost);  
  
  app->c2p_ctx = app->mapc2p = 0;  
  app->has_mapc2p = mom->mapc2p ? true : false;

  if (app->has_mapc2p) {
    // initialize computational to physical space mapping
    app->c2p_ctx = mom->c2p_ctx;
    app->mapc2p = mom->mapc2p;

    // we project mapc2p on p=1 basis functions
    struct gkyl_basis basis;
    gkyl_cart_modal_tensor(&basis, ndim, 1);

    // initialize DG field representing mapping
    struct gkyl_array *c2p = mkarr(false, ndim*basis.num_basis, app->local_ext.volume);
    gkyl_eval_on_nodes *ev_c2p = gkyl_eval_on_nodes_new(&app->grid, &basis, ndim, mom->mapc2p, mom->c2p_ctx);
    gkyl_eval_on_nodes_advance(ev_c2p, 0.0, &app->local_ext, c2p);

    // write DG projection of mapc2p to file
    cstr fileNm = cstr_from_fmt("%s-mapc2p.gkyl", app->name);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, 0, c2p, fileNm.str);
    cstr_drop(&fileNm);

    gkyl_array_release(c2p);
    gkyl_eval_on_nodes_release(ev_c2p);
  }

  // create geometry object (no GPU support in fluids right now JJ: 11/26/23)
  app->geom = gkyl_wave_geom_new(&app->grid, &app->local_ext,
    app->mapc2p, app->c2p_ctx, false);

  double cfl_frac = mom->cfl_frac == 0 ? 0.95 : mom->cfl_frac;
  app->cfl = 1.0*cfl_frac;
  if (app->scheme_type == GKYL_MOMENT_MP)
    app->cfl = 0.4*cfl_frac; // this should be 1/(1+alpha) = 0.2 but is set to a larger value

  app->num_periodic_dir = mom->num_periodic_dir;
  for (int d=0; d<ndim; ++d)
    app->periodic_dirs[d] = mom->periodic_dirs[d];

  // construct list of directions to skip
  for (int d=0; d<3; ++d)
    app->is_dir_skipped[d] = 0;
  for (int i=0; i<mom->num_skip_dirs; ++i)
    app->is_dir_skipped[mom->skip_dirs[i]] = 1;

  app->has_field = 0;
  // Are we running with a field?
  if (mom->field.init) {
    app->has_field = 1;
  }
  // Initialize a (potentially null) field object for safety.
  moment_field_init(mom, &mom->field, app, &app->field);

  // Are we running with Braginskii transport?
  app->has_braginskii = mom->has_braginskii;
  app->coll_fac = mom->coll_fac;

  int ns = app->num_species = mom->num_species;
  // allocate space to store species objects
  app->species = ns>0 ? gkyl_calloc(ns, sizeof(struct moment_species)) : 0;
  // create species
  for (int i=0; i<ns; ++i) {
    moment_species_init(mom, &mom->species[i], app, &app->species[i]);
  }

  // specify collision parameters in the exposed app
  app->has_collision = mom->has_collision;
  int num_entries = app->num_species * (app->num_species-1) / 2;
  for (int s=0; s<app->num_species; ++s)
    for (int r=0; r<app->num_species; ++r)
      app->nu_base[s][r] = mom->nu_base[s][r];

  // There are a significant number of options which necessitate the source solve in fluids
  // (e.g. applied acceleration, geometric sources, multi-species transport, electromagnetic coupling, etc.).
  // To facilitate these options, each fluid species stores a boolean for whether or not it will require a 
  // a source update and sets the value to true if any input parameter would necessitate a source, from
  // app->has_field (because the field couples to the fluid via sources) to whether an applied accerleation
  // pointer is set to various closure-related input parameters. 
  // We also check if there are sources in EM-only simulations, such as volume expansion or applied current
  // simulations. 
  app->update_sources = false;
  if (app->field.has_volume_sources || app->field.has_app_current) {
    app->update_sources = true; 
  }
  for (int s=0; s<app->num_species; ++s) {
    if (app->species[s].update_sources) {
      app->update_sources = true; 
    }
  }
  if (app->update_sources) {
    moment_coupling_init(app, &app->sources);
  }

  app->update_mhd_source = false;
  if (ns==1 && mom->species[0].equation->type==GKYL_EQN_MHD) {
    app->update_mhd_source = true;
    mhd_src_init(app, &mom->species[0], &app->mhd_source);
  }

  // allocate work array for use in MP scheme
  if (app->scheme_type == GKYL_MOMENT_MP || app->scheme_type == GKYL_MOMENT_KEP) {
    int max_eqn = 0;
    for (int i=0; i<ns; ++i)
      max_eqn = int_max(max_eqn, app->species[i].num_equations);
    if (app->has_field)
      max_eqn = int_max(max_eqn, 8); // maxwell equations have 8 components
    app->ql = mkarr(false, max_eqn, app->local_ext.volume);
    app->qr = mkarr(false, max_eqn, app->local_ext.volume);

    app->amdq = mkarr(false, max_eqn, app->local_ext.volume);
    app->apdq = mkarr(false, max_eqn, app->local_ext.volume);
  }

  // initialize stat object to all zeros
  app->stat = (struct gkyl_moment_stat) {
  };

  return app;
}

double
gkyl_moment_app_max_dt(gkyl_moment_app* app)
{
  double max_dt = DBL_MAX;
  for (int i=0;  i<app->num_species; ++i) 
    max_dt = fmin(max_dt, moment_species_max_dt(app, &app->species[i]));

  if (app->has_field)
    max_dt = fmin(max_dt, moment_field_max_dt(app, &app->field));

  double max_dt_global;
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &max_dt, &max_dt_global);

  return max_dt_global;
}

void
gkyl_moment_app_apply_ic(gkyl_moment_app* app, double t0)
{
  app->tcurr = t0;
  gkyl_moment_app_apply_ic_field(app, t0);
  for (int i=0;  i<app->num_species; ++i)
    gkyl_moment_app_apply_ic_species(app, i, t0);
}

void
gkyl_moment_app_apply_ic_field(gkyl_moment_app* app, double t0)
{
  if (app->has_field != 1) return;
  
  app->tcurr = t0;
  int num_quad = app->scheme_type == GKYL_MOMENT_MP ? 4 : 2;
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, num_quad, 8, app->field.init, app->field.ctx);
  
  gkyl_fv_proj_advance(proj, t0, &app->local, app->field.fcurr);
  gkyl_fv_proj_release(proj);

  if (app->field.has_ext_em) {
    gkyl_fv_proj_advance(app->field.ext_em_proj, t0, &app->local, app->field.ext_em);
  }
  if (app->field.has_app_current) {
    gkyl_fv_proj_advance(app->field.app_current_proj, t0, &app->local, app->field.app_current);
  }

  moment_field_apply_bc(app, t0, &app->field, app->field.fcurr);
}

void
gkyl_moment_app_apply_ic_species(gkyl_moment_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  app->tcurr = t0;
  int num_quad = app->scheme_type == GKYL_MOMENT_MP ? 4 : 2;  
  gkyl_fv_proj *proj = gkyl_fv_proj_new(&app->grid, num_quad, app->species[sidx].num_equations,
    app->species[sidx].init, app->species[sidx].ctx);
  
  gkyl_fv_proj_advance(proj, t0, &app->local, app->species[sidx].fcurr);
  gkyl_fv_proj_release(proj);

  if (app->species[sidx].has_app_accel) {
    gkyl_fv_proj_advance(app->species[sidx].app_accel_proj, t0, &app->local, app->species[sidx].app_accel);
  }

  moment_species_apply_bc(app, t0, &app->species[sidx], app->species[sidx].fcurr);
}

void
gkyl_moment_app_write(const gkyl_moment_app* app, double tm, int frame)
{
  gkyl_moment_app_write_field(app, tm, frame);
  for (int i=0; i<app->num_species; ++i)
    gkyl_moment_app_write_species(app, i, tm, frame);
}

void
gkyl_moment_app_write_field(const gkyl_moment_app* app, double tm, int frame)
{
  if (app->has_field != 1) return;

  struct gkyl_msgpack_data *mt = moment_array_meta_new( (struct moment_output_meta) {
      .frame = frame,
      .stime= tm
    }
  );

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, app->field.fcurr, fileNm.str);
  cstr_drop(&fileNm);

  // write external EM field if it is present
  if (app->field.ext_em_proj) {
    if (app->field.ext_em_evolve || frame == 0) {
      cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "ext_em", frame);
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, app->field.ext_em, fileNm.str);
      cstr_drop(&fileNm);
    }
  }

  // write applied currents if it is present
  if (app->field.app_current_proj) {
    if (app->field.app_current_evolve || frame == 0) {
      cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "app_current", frame);
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, app->field.app_current, fileNm.str);
      cstr_drop(&fileNm);
    }
  }

  moment_array_meta_release(mt);
}

void
gkyl_moment_app_write_field_energy(gkyl_moment_app *app)
{
  if (app->has_field) {
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out field energy
      cstr fileNm = cstr_from_fmt("%s-field-energy.gkyl", app->name);
      
      if (app->field.is_first_energy_write_call) {
        // write to a new file (this ensure previous output is removed)
        gkyl_dynvec_write(app->field.integ_energy, fileNm.str);
        app->field.is_first_energy_write_call = false;
      }
      else {
        // append to existing file
        gkyl_dynvec_awrite(app->field.integ_energy, fileNm.str);
      }
      cstr_drop(&fileNm);
    }
    gkyl_dynvec_clear(app->field.integ_energy);
  }
}

void
gkyl_moment_app_write_integrated_mom(gkyl_moment_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    if (rank == 0) {
      // write out diagnostic moments
      cstr fileNm = cstr_from_fmt("%s-%s-%s.gkyl", app->name, app->species[i].name,
        "imom");
      
      if (app->species[i].is_first_q_write_call) {
        gkyl_dynvec_write(app->species[i].integ_q, fileNm.str);
        app->species[i].is_first_q_write_call = false;
      }
      else {
        gkyl_dynvec_awrite(app->species[i].integ_q, fileNm.str);
      }
      cstr_drop(&fileNm);
    }
    gkyl_dynvec_clear(app->species[i].integ_q);
  }
}

void
gkyl_moment_app_write_species(const gkyl_moment_app* app, int sidx, double tm, int frame)
{
  struct gkyl_msgpack_data *mt = moment_array_meta_new( (struct moment_output_meta) {
      .frame = frame,
      .stime = tm
    }
  );
  
  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, app->species[sidx].name, frame);
  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, app->species[sidx].fcurr, fileNm.str);
  cstr_drop(&fileNm);

  if (app->scheme_type == GKYL_MOMENT_KEP) {
    cstr fileNm = cstr_from_fmt("%s-%s-alpha_%d.gkyl", app->name, app->species[sidx].name, frame);
    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, app->species[sidx].alpha, fileNm.str);
    cstr_drop(&fileNm);
  }

  if (app->species[sidx].has_app_accel) {
    if (app->species[sidx].app_accel_evolve || frame == 0) {
      cstr fileNm = cstr_from_fmt("%s-%s-app_accel_%d.gkyl", app->name, app->species[sidx].name, frame);
      gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, app->species[sidx].app_accel, fileNm.str);
      cstr_drop(&fileNm);      
    }
  }

  moment_array_meta_release(mt);
}

struct gkyl_update_status
gkyl_moment_update(gkyl_moment_app* app, double dt)
{
  gkyl_comm_barrier(app->comm);
  
  app->stat.nup += 1;
  
  struct timespec wst = gkyl_wall_clock();
  struct gkyl_update_status status = app->update_func(app, dt);
  app->tcurr += status.dt_actual;
  
  app->stat.total_tm += gkyl_time_diff_now_sec(wst);
  
  return status;
}

int
gkyl_moment_app_field_energy_ndiag(gkyl_moment_app *app)
{
  return 6;
}

void
gkyl_moment_app_get_field_energy(gkyl_moment_app *app, double *vals)
{
  double energy_global[6] = { 0.0 };
  if (app->has_field) {
    double energy[6];
    calc_integ_quant(app->field.maxwell, app->grid.cellVolume, app->field.fcurr, app->geom,
      app->local, energy);
    
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, 6, energy, energy_global);
  }
  for (int i=0; i<6; ++i)
    vals[i] = energy_global[i];
}

void
gkyl_moment_app_calc_field_energy(gkyl_moment_app* app, double tm)
{
  if (app->has_field) {
    double energy[6] = { 0.0 };
    gkyl_moment_app_get_field_energy(app, energy);
    gkyl_dynvec_append(app->field.integ_energy, tm, energy);

  }
}

void
gkyl_moment_app_calc_integrated_mom(gkyl_moment_app *app, double tm)
{
  for (int sidx=0; sidx<app->num_species; ++sidx) {

    int num_diag = app->species[sidx].equation->num_diag;
    double q_integ[num_diag];

    calc_integ_quant(app->species[sidx].equation, app->grid.cellVolume, app->species[sidx].fcurr, app->geom,
      app->local, q_integ);

    double q_integ_global[num_diag];
    gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_SUM, num_diag, q_integ, q_integ_global);
    gkyl_dynvec_append(app->species[sidx].integ_q, tm, q_integ_global);
  }
}

void
gkyl_moment_app_nghost(gkyl_moment_app *app, int nghost[3])
{
  for (int i=0; i<app->ndim; ++i)
    nghost[i] = app->nghost[i];
}

struct gkyl_array*
gkyl_moment_app_get_write_array_species(const gkyl_moment_app* app, int sidx)
{
  // this needs to be consistent with the write_species method
  return app->species[sidx].fcurr;
}

struct gkyl_array*
gkyl_moment_app_get_write_array_field(const gkyl_moment_app* app)
{
  if (app->has_field != 1) return 0;
// this needs to be consistent with the write_field method
  return app->field.fcurr;
}

struct gkyl_moment_stat
gkyl_moment_app_stat(gkyl_moment_app* app)
{
  return app->stat;
}

// ensure stats across processors are made consistent
static void
comm_reduce_app_stat(const gkyl_moment_app* app, const struct gkyl_moment_stat *local,
  struct gkyl_moment_stat *global)
{
  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz == 1) {
    memcpy(global, local, sizeof(struct gkyl_moment_stat));
    return;
  }

  enum { NUP, NFAIL, NFEULER, NSTAGE_2_FAIL, NSTAGE_3_FAIL, L_END };
  int64_t l_red[] = {
    [NUP] = local->nup,
    [NFAIL] = local->nfail,
    [NFEULER] = local->nfeuler,
    [NSTAGE_2_FAIL] = local->nstage_2_fail,
    [NSTAGE_3_FAIL] = local->nstage_3_fail
  };

  int64_t l_red_global[L_END];
  gkyl_comm_allreduce(app->comm, GKYL_INT_64, GKYL_MAX, L_END, l_red, l_red_global);

  global->nup = l_red_global[NUP];
  global->nfail = l_red_global[NFAIL];
  global->nfeuler = l_red_global[NFEULER];
  global->nstage_2_fail = l_red_global[NSTAGE_2_FAIL];
  global->nstage_3_fail = l_red_global[NSTAGE_3_FAIL];

  enum { TOTAL_TM, SPECIES_TM, FIELD_TM, SOURCES_TM, INIT_SPECIES_TM, INIT_FIELD_TM,
    SPECIES_RHS_TM, FIELD_RHS_TM, SPECIES_BC_TM, FIELD_BC_TM,
    D_END
  };

  double d_red[] = {
    [TOTAL_TM] = local->total_tm,
    [SPECIES_TM] = local->species_tm,
    [FIELD_TM] = local->field_tm,
    [SOURCES_TM] = local->sources_tm,
    [INIT_SPECIES_TM] = local->init_species_tm,
    [INIT_FIELD_TM] = local->init_field_tm,
    [SPECIES_RHS_TM] = local->species_rhs_tm,
    [FIELD_RHS_TM] = local->field_rhs_tm,
    [SPECIES_BC_TM] = local->species_bc_tm,
    [FIELD_BC_TM] = local->field_bc_tm
  };

  double_t d_red_global[D_END];
  gkyl_comm_allreduce(app->comm, GKYL_DOUBLE, GKYL_MAX, D_END, d_red, d_red_global);

  global->total_tm = d_red_global[TOTAL_TM];
  global->species_tm = d_red_global[SPECIES_TM];
  global->field_tm = d_red_global[FIELD_TM];
  global->sources_tm = d_red_global[SOURCES_TM];
  global->init_species_tm = d_red_global[INIT_SPECIES_TM];
  global->init_field_tm = d_red_global[INIT_FIELD_TM];
  global->species_rhs_tm = d_red_global[SPECIES_RHS_TM];
  global->field_rhs_tm = d_red_global[FIELD_RHS_TM];
  global->species_bc_tm = d_red_global[SPECIES_BC_TM];
  global->field_bc_tm = d_red_global[FIELD_BC_TM];
}

static void
comm_reduce_wave_prop_stats(const gkyl_moment_app* app, const struct gkyl_wave_prop_stats *local,
  struct gkyl_wave_prop_stats *global)
{
  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz == 1) {
    memcpy(global, local, sizeof(struct gkyl_wave_prop_stats));
    return;
  }

  enum { N_CALLS, N_BAD_ADVANCE_CALLS, N_MAX_BAD_CELLS, L_END };
  int64_t l_red[] = {
    [N_CALLS] = local->n_calls,
    [N_BAD_ADVANCE_CALLS] = local->n_bad_advance_calls,
    [N_MAX_BAD_CELLS] = local->n_max_bad_cells
  };

  int64_t l_red_global[L_END];
  gkyl_comm_allreduce(app->comm, GKYL_INT_64, GKYL_MAX, L_END, l_red, l_red_global); 

  global->n_calls = l_red_global[N_CALLS];
  global->n_bad_advance_calls = l_red_global[N_BAD_ADVANCE_CALLS];
  global->n_max_bad_cells = l_red_global[N_MAX_BAD_CELLS];

  int64_t n_bad_cells_local = local->n_bad_cells, n_bad_cells = 0;
  
  gkyl_comm_allreduce(app->comm, GKYL_INT_64, GKYL_SUM, 1, &n_bad_cells_local, &n_bad_cells);
  global->n_bad_cells = n_bad_cells;
}

void
gkyl_moment_app_stat_write(const gkyl_moment_app* app)
{
  cstr fileNm = cstr_from_fmt("%s-%s", app->name, "stat.json");

  char buff[70];
  time_t t = time(NULL);
  struct tm curr_tm = *localtime(&t);

  // total number of cells updated in simulation
  long tot_cells_up = app->global.volume*app->num_species*app->ndim*app->stat.nup;

  struct gkyl_moment_stat stat = { };
  comm_reduce_app_stat(app, &app->stat, &stat);

  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  int num_ranks;
  gkyl_comm_get_size(app->comm, &num_ranks);

  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  if (rank == 0) fp = fopen(fileNm.str, "a");

  gkyl_moment_app_cout(app, fp, "{\n");

  if (strftime(buff, sizeof buff, "%c", &curr_tm))
    gkyl_moment_app_cout(app, fp, " date : %s\n", buff);

  gkyl_moment_app_cout(app, fp, " num_ranks : %d,\n", num_ranks);
  
  gkyl_moment_app_cout(app, fp, " nup : %ld,\n", stat.nup);
  gkyl_moment_app_cout(app, fp, " nfail : %ld,\n", stat.nfail);
  gkyl_moment_app_cout(app, fp, " total_tm : %lg,\n", stat.total_tm);

  if (app->scheme_type == GKYL_MOMENT_WAVE_PROP) {
    gkyl_moment_app_cout(app, fp, " species_tm : %lg,\n", stat.species_tm);
    gkyl_moment_app_cout(app, fp, " field_tm : %lg,\n", stat.field_tm);
    gkyl_moment_app_cout(app, fp, " sources_tm : %lg\n", stat.sources_tm);
  }
  else if (app->scheme_type == GKYL_MOMENT_MP || app->scheme_type == GKYL_MOMENT_KEP) {
    
    gkyl_moment_app_cout(app, fp, " nfeuler : %ld,\n", stat.nfeuler);
    gkyl_moment_app_cout(app, fp, " nstage_2_fail : %ld,\n", stat.nstage_2_fail);
    gkyl_moment_app_cout(app, fp, " nstage_3_fail : %ld,\n", stat.nstage_3_fail);

    gkyl_moment_app_cout(app, fp, " stage_2_dt_diff : [ %lg, %lg ],\n",
      stat.stage_2_dt_diff[0], stat.stage_2_dt_diff[1]);
    gkyl_moment_app_cout(app, fp, " stage_3_dt_diff : [ %lg, %lg ],\n",
      stat.stage_3_dt_diff[0], stat.stage_3_dt_diff[1]);

    gkyl_moment_app_cout(app, fp, " total_tm : %lg,\n", stat.total_tm);
    gkyl_moment_app_cout(app, fp, " init_species_tm : %lg,\n", stat.init_species_tm);
    if (app->has_field)
      gkyl_moment_app_cout(app, fp, " init_field_tm : %lg,\n", stat.init_field_tm);

    gkyl_moment_app_cout(app, fp, " species_rhs_tm : %lg,\n", stat.species_rhs_tm);

    if (app->has_field)
      gkyl_moment_app_cout(app, fp, " field_rhs_tm : %lg,\n", stat.field_rhs_tm);
  }

  gkyl_moment_app_cout(app, fp, " species_bc_tm : %lg,\n", stat.species_bc_tm);
  gkyl_moment_app_cout(app, fp, " field_bc_tm : %lg,\n", stat.field_bc_tm);    

  for (int i = 0; i < app->num_species; ++i) {
    long tot_bad_cells = 0L;
    
    if (app->scheme_type == GKYL_MOMENT_WAVE_PROP) {
      for (int d = 0; d < app->ndim; ++d) {
        struct gkyl_wave_prop_stats wvs_local = gkyl_wave_prop_stats(app->species[i].slvr[d]);
        struct gkyl_wave_prop_stats wvs = {};
        comm_reduce_wave_prop_stats(app, &wvs_local, &wvs);

        gkyl_moment_app_cout(app, fp, " %s_n_bad_1D_sweeps[%d] = %ld\n",
          app->species[i].name, d, wvs.n_bad_advance_calls);
        gkyl_moment_app_cout(app, fp, " %s_n_bad_cells[%d] = %ld\n",
          app->species[i].name, d, wvs.n_bad_cells);
        gkyl_moment_app_cout(app, fp, " %s_n_max_bad_cells[%d] = %ld\n",
          app->species[i].name, d, wvs.n_max_bad_cells);

        tot_bad_cells += wvs.n_bad_cells;
      }
    }
    gkyl_moment_app_cout(app, fp, " %s_bad_cell_frac = %lg\n",
      app->species[i].name, (double)tot_bad_cells/tot_cells_up);
  }

  gkyl_moment_app_cout(app, fp, "}\n");

  if (rank == 0)
    fclose(fp);

  cstr_drop(&fileNm);
}

static struct gkyl_app_restart_status
header_from_file(gkyl_moment_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = { .io_status = 0 };
  
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (!gkyl_rect_grid_cmp(&app->grid, &grid))
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
      if (hdr.etype != GKYL_DOUBLE)
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
    }

    struct moment_output_meta meta =
      moment_meta_from_mpack( &(struct gkyl_msgpack_data) {
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
gkyl_moment_app_from_file_field(gkyl_moment_app *app, const char *fname)
{
  if (app->has_field != 1)
    return (struct gkyl_app_restart_status) {
      .io_status = GKYL_ARRAY_RIO_SUCCESS,
      .frame = 0,
      .stime = 0.0
    };

  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->field.fcurr, fname);
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      moment_field_apply_bc(app, rstat.stime, &app->field, app->field.fcurr);
    }
  }

  // Compute external EM field and applied current if present
  // Computation necessary in case external EM field or applied current
  // are time-independent and not computed in the time-stepping loop
  // since they are not read-in as part of restarts. 
  if (app->field.has_ext_em) {
    gkyl_fv_proj_advance(app->field.ext_em_proj, rstat.stime, &app->local, app->field.ext_em);
  }
  if (app->field.has_app_current) {
    gkyl_fv_proj_advance(app->field.app_current_proj, rstat.stime, &app->local, app->field.app_current);
  }  

  return rstat;
}

struct gkyl_app_restart_status 
gkyl_moment_app_from_file_species(gkyl_moment_app *app, int sidx,
  const char *fname)
{
  struct gkyl_app_restart_status rstat = header_from_file(app, fname);
  
  if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->species[sidx].fcurr, fname);
    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      moment_species_apply_bc(app, rstat.stime, &app->species[sidx], app->species[sidx].fcurr);
    }
  }

  // Compute applied acceleration if present.
  // Computation necessary in case applied acceleration
  // is time-independent and not computed in the time-stepping loop
  // since it is not read-in as part of restarts. 
  if (app->species[sidx].has_app_accel) {
    gkyl_fv_proj_advance(app->species[sidx].app_accel_proj, rstat.stime, &app->local, app->species[sidx].app_accel);
  }  

  return rstat;
}

struct gkyl_app_restart_status
gkyl_moment_app_from_frame_field(gkyl_moment_app *app, int frame)
{
  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
  struct gkyl_app_restart_status rstat = gkyl_moment_app_from_file_field(app, fileNm.str);
  app->field.is_first_energy_write_call = false; // append to existing diagnostic
  cstr_drop(&fileNm);
  
  return rstat;
}

struct gkyl_app_restart_status
gkyl_moment_app_from_frame_species(gkyl_moment_app *app, int sidx, int frame)
{
  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, app->species[sidx].name, frame);
  struct gkyl_app_restart_status rstat = gkyl_moment_app_from_file_species(app, sidx, fileNm.str);
  app->species[sidx].is_first_q_write_call = false; // append to existing diagnostic
  cstr_drop(&fileNm);
  
  return rstat;
}

struct gkyl_app_restart_status
gkyl_moment_app_read_from_frame(gkyl_moment_app *app, int frame)
{
  struct gkyl_app_restart_status rstat;

  rstat = gkyl_moment_app_from_frame_field(app, frame);

  for (int i = 0; i < app->num_species; i++) {
    rstat = gkyl_moment_app_from_frame_species(app, i, frame);
  }

  return rstat;
}

// private function to handle variable argument list for printing
static void
v_moment_app_cout(const gkyl_moment_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp) {
    vfprintf(fp, fmt, argp);
    fflush(fp);
  }
}

void
gkyl_moment_app_cout(const gkyl_moment_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_moment_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_moment_app_release(gkyl_moment_app* app)
{
  if (app->update_sources)
    moment_coupling_release(app, &app->sources);

  gkyl_comm_release(app->comm);
  gkyl_rect_decomp_release(app->decomp);
  
  for (int i=0; i<app->num_species; ++i)
    moment_species_release(&app->species[i]);
  gkyl_free(app->species);

  moment_field_release(&app->field);

  if (app->update_mhd_source)
    mhd_src_release(&app->mhd_source);

  gkyl_wave_geom_release(app->geom);

  if (app->scheme_type == GKYL_MOMENT_MP || app->scheme_type == GKYL_MOMENT_KEP) {
    gkyl_array_release(app->ql);
    gkyl_array_release(app->qr);
    gkyl_array_release(app->amdq);
    gkyl_array_release(app->apdq);
  }

  gkyl_free(app);
}

