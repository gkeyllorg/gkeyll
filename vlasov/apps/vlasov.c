#include <stdarg.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_basis.h>
#include <gkyl_comm_io.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_null_comm.h>

#include <gkyl_vlasov_priv.h>
#include <gkyl_app_priv.h>

#include <mpack.h>

// returned gkyl_array_meta must be freed using vlasov_array_meta_release
struct gkyl_msgpack_data*
vlasov_array_meta_new(struct vlasov_output_meta meta)
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

void
vlasov_array_meta_release(struct gkyl_msgpack_data *mt)
{
  if (!mt) return;
  MPACK_FREE(mt->meta);
  gkyl_free(mt);
}

struct vlasov_output_meta
vlasov_meta_from_mpack(struct gkyl_msgpack_data *mt)
{
  struct vlasov_output_meta meta = { .frame = 0, .stime = 0.0 };

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

gkyl_vlasov_app*
gkyl_vlasov_app_new(struct gkyl_vm *vm)
{
  disable_denorm_float();

  assert(vm->num_species <= GKYL_MAX_SPECIES);

  gkyl_vlasov_app *app = gkyl_malloc(sizeof(gkyl_vlasov_app));

  int ntot = vm->num_species, ns = 0, nsf = 0;
  for (int i=0; i<ntot; ++i) {
    if (vm->species[i].type != GKYL_SPECIES_FLUID) ns += 1;
    if (vm->species[i].type != GKYL_SPECIES_VLASOV) nsf += 1;
  }

  int cdim = app->cdim = vm->cdim;
  int vdim = app->vdim = vm->vdim;
  int pdim = cdim+vdim;
  int poly_order = app->poly_order = vm->poly_order;
  app->num_species = ntot;
  app->num_kinetic_species = ns;
  app->num_fluid_species = nsf;

  double cfl_frac = vm->cfl_frac == 0 ? 1.0 : vm->cfl_frac;
  app->cfl = cfl_frac;

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = vm->parallelism.use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif

  app->num_periodic_dir = vm->num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = vm->periodic_dirs[d];

  strcpy(app->name, vm->name);
  app->tcurr = 0.0; // reset on init

  if (app->use_gpu) {
    // allocate device basis if we are using GPUs
    app->basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  }
  else {
    app->basis_on_dev = &app->basis;
  }

  // basis functions
  switch (vm->basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&app->basis, cdim, poly_order);
      if (app->use_gpu) {
        gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev, cdim, poly_order);
      }
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      gkyl_cart_modal_tensor(&app->basis, cdim, poly_order);
      if (app->use_gpu) {
        gkyl_cart_modal_tensor_cu_dev(app->basis_on_dev, cdim, poly_order);
      }
      break;

    default:
      assert(false);
      break;
  }

  gkyl_rect_grid_init(&app->grid, cdim, vm->lower, vm->upper, vm->cells);

  int ghost[] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&app->grid, ghost, &app->global_ext, &app->global);

  if (vm->parallelism.comm == 0) {
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
    app->decomp = gkyl_rect_decomp_new_from_cuts(app->cdim, vm->parallelism.cuts, &app->global);

    // Create a new communicator with the decomposition in it.
    app->comm = gkyl_comm_split_comm(vm->parallelism.comm, 0, app->decomp);

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
  // Note: *only* uses a p=1 DG representation of the geometry (JJ: 11/24/23)
  app->c2p_ctx = app->mapc2p = 0;  
  app->has_mapc2p = vm->mapc2p ? true : false;

  if (app->has_mapc2p) {
    // initialize computational to physical space mapping
    app->c2p_ctx = vm->c2p_ctx;
    app->mapc2p = vm->mapc2p;

    // we project mapc2p on p=1 basis functions
    struct gkyl_basis basis;
    gkyl_cart_modal_tensor(&basis, cdim, 1);

    // initialize DG field representing mapping
    struct gkyl_array *c2p = mkarr(false, cdim*basis.num_basis, app->local_ext.volume);
    gkyl_eval_on_nodes *ev_c2p = gkyl_eval_on_nodes_new(&app->grid, &basis, cdim, vm->mapc2p, vm->c2p_ctx);
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

  // Create the configuration-space position map (C^0 piecewise linear,
  // diagonal). Always created; directions without a user map are the identity.
  // Shared by all species, which acquire a reference.
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = { 0 };
  for (int d=0; d<cdim; ++d) {
    inp_pmap[d].eval_pmap = vm->mapc2p_pos[d].mapc2p_pos_func;
    inp_pmap[d].ctx = vm->mapc2p_pos[d].mapc2p_pos_ctx;
  }
  app->pos_map = gkyl_vlasov_position_map_new(&app->grid, &app->local,
    &app->local_ext, &app->basis, inp_pmap, app->use_gpu);

  // Finalize the position-map ghost cells. The constructor filled every ghost
  // with the adjacent local-skin value (a copy BC, valid at physical
  // boundaries); now overwrite the periodic and inter-rank ghosts with the
  // wrapped/neighbor values, mirroring the conf-array sync in vm_*_apply_bc.
  struct gkyl_array *pos_map_arrs[] = {
    app->pos_map->pmap, app->pos_map->jacob_pos,
    app->pos_map->jacob_pos_surf, app->pos_map->jacob_pos_gauss
  };
  for (int ia=0; ia<4; ++ia) {
    gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
      app->num_periodic_dir, app->periodic_dirs, pos_map_arrs[ia]);
    gkyl_comm_array_sync(app->comm, &app->local, &app->local_ext, pos_map_arrs[ia]);
  }

  // The position map is static in time; write it once here (uniform/identity
  // grids included), mirroring the mapc2p write above. Both the p=3 map and its
  // p=0 cell average are written, with metadata built from the map's I/O basis.
  gkyl_vlasov_position_map_write(app->pos_map, app->comm, app->name,
    "position-map");

  // allocate space to store vlasov-maxwell geometry objects
  app->vm_geom = gkyl_malloc(sizeof(struct vm_geom));
  vm_geom_init(vm, app, app->vm_geom);

  app->has_field = !vm->skip_field; // note inversion of truth value
  // A field object is always created: the null field when no field is present.
  app->field = vlasov_field_new(vm, app);

  app->species = ntot>0 ? gkyl_malloc(sizeof(struct vlasov_species[ntot])) : 0;

  // Construct each species container (allocates its aspect, stores the input,
  // and sets its methods); the aspects are initialized below, once all
  // containers exist, since the initializations look species up by name.
  for (int i=0; i<ntot; ++i)
    vlasov_species_init(app, &vm->species[i], &app->species[i]);

  // initialize each kinetic species
  for (int i=0; i<ntot; ++i)
    if (app->species[i].kinetic)
      vm_species_init(vm, app, app->species[i].kinetic);

  // initialize species wall emission terms: these rely
  // on other species which must be allocated in the previous step
  for (int i=0; i<ntot; ++i) {
    struct vm_species *vms = app->species[i].kinetic;
    if (!vms) continue;
    if (vms->emit_lo)
      vm_species_emission_cross_init(app, vms, &vms->bc_emission_lo);
    if (vms->emit_up)
      vm_species_emission_cross_init(app, vms, &vms->bc_emission_up);
  }
  
  // initialize each species cross-species terms: this has to be done here
  // as need pointers to colliding species' collision objects
  // allocated in the previous step
  for (int i=0; i<ntot; ++i) {
    struct vm_species *vms = app->species[i].kinetic;
    if (!vms) continue;
    vm_species_lbo_cross_init(app, vms, &vms->lbo);
    vm_species_bgk_cross_init(app, vms, &vms->bgk);
  }

  // initialize each species source terms: this has to be done here
  // as they may initialize a bflux updater for their source species.
  for (int i=0; i<ntot; ++i) {
    struct vm_species *vms = app->species[i].kinetic;
    if (vms && vms->source_id)
      vm_species_source_init(app, vms, &vms->src);
  }

  // initialize each fluid species
  // Fluid species must be initialized after kinetic species, as some fluid species couple
  // to kinetic species and pointers are allocated by the kinetic species objects
  for (int i=0; i<ntot; ++i)
    if (app->species[i].fluid)
      vm_fluid_species_init(vm, app, app->species[i].fluid);

  for (int i=0; i<ntot; ++i) {
    struct vm_fluid_species *vmf = app->species[i].fluid;
    if (vmf && vmf->source_id)
      vm_fluid_species_source_init(app, vmf, &vmf->src);
  }

  // Check if there are both any fluid species and an EM field. 
  // If there are, initialize the implicit fluid-EM coupling solver.
  app->has_fluid_em_coupling = false;
  if (nsf > 0 && app->has_field) {
    app->has_fluid_em_coupling = true;
    app->fl_em = vm_fluid_em_coupling_new(app);
  }

  // Use implicit BGK collisions if any species requests them.
  app->has_implicit_coll_scheme = false;
  for (int i=0; i<ntot; ++i){
    if (app->species[i].kinetic && app->species[i].kinetic->info.collisions.is_implicit){
      app->has_implicit_coll_scheme = true;
    }
  }

  // Set the appropriate update function for taking a single time step
  // If we have implicit fluid-EM coupling or implicit BGK collisions, 
  // we perform a first-order operator split and treat those terms implicitly.
  // Otherwise, we default to an SSP-RK3 method. 
  if (app->has_implicit_coll_scheme || app->has_fluid_em_coupling) {
    app->update_func = vlasov_update_op_split;
  }
  else {
    app->update_func = vlasov_update_ssp_rk3;
  }

  // initialize stat object
  app->stat = (struct gkyl_vlasov_stat) {
    .use_gpu = app->use_gpu,
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };

  return app;
}

int
vlasov_find_species_idx(const gkyl_vlasov_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].name) == 0)
      return i;
  return -1;
}

struct vm_species *
vm_find_species(const gkyl_vlasov_app *app, const char *nm)
{
  int i = vm_find_species_idx(app, nm);
  return i >= 0 ? app->species[i].kinetic : 0;
}

int
vm_find_species_idx(const gkyl_vlasov_app *app, const char *nm)
{
  int i = vlasov_find_species_idx(app, nm);
  return (i >= 0 && app->species[i].kinetic) ? i : -1;
}

void
vm_apply_bc(gkyl_vlasov_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array *fluid[], struct gkyl_array *emfield)
{
  // distf[] and fluid[] are indexed over the overall species count.
  int num_species = app->num_species;
  for (int i=0; i<num_species; ++i)
    vlasov_species_apply_bc(app, &app->species[i], distf[i], fluid[i], tcurr);
  // No-op for Vlasov-Poisson and the null field (no EM boundary conditions).
  vlasov_field_apply_bc(app, emfield);
}

void
gkyl_vlasov_app_apply_ic(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_species; ++i)
    vlasov_species_apply_ic(app, &app->species[i], t0);

  gkyl_vlasov_app_apply_ic_field(app, t0); // no-op for the null field

  // BCs must be done after all species initialize for emission BCs to work.
  int num_species = app->num_species;
  struct gkyl_array *distf[num_species];
  struct gkyl_array *fluid[num_species];
  vlasov_species_gather_state(app, distf, fluid);
  vm_apply_bc(app, t0, distf, fluid, app->field->em);
}

void
gkyl_vlasov_app_apply_ic_field(gkyl_vlasov_app* app, double t0)
{
  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();

  // Dispatches to the Maxwell or Poisson IC (no-op for the null field).
  const struct gkyl_array *fin[app->num_species];
  vlasov_species_gather_dist(app, fin);
  vlasov_field_apply_ic(app, fin, t0);

  app->stat.init_field_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_vlasov_app_apply_ic_species(gkyl_vlasov_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);
  vlasov_species_apply_ic(app, &app->species[sidx], t0);
}

int
gkyl_vlasov_app_find_species(const gkyl_vlasov_app* app, const char *nm)
{
  return vlasov_find_species_idx(app, nm);
}

void
gkyl_vlasov_app_calc_integrated_mom(gkyl_vlasov_app* app, double tm)
{
  for (int i=0; i<app->num_species; ++i)
    vlasov_species_calc_integrated_mom(app, &app->species[i], tm);
}

void
gkyl_vlasov_app_calc_integrated_L2_f(gkyl_vlasov_app* app, double tm)
{
  for (int i=0; i<app->num_species; ++i)
    vlasov_species_calc_integrated_L2_f(app, &app->species[i], tm);
}

void
gkyl_vlasov_app_calc_field_energy(gkyl_vlasov_app* app, double tm)
{
  // The distributions let Vlasov-Poisson solve for the potential at time tm.
  const struct gkyl_array *fin[app->num_species];
  vlasov_species_gather_dist(app, fin);
  vlasov_field_calc_energy(app, tm, fin); // no-op for the null field
}

void
gkyl_vlasov_app_write(gkyl_vlasov_app* app, double tm, int frame)
{
  gkyl_vlasov_app_write_field(app, tm, frame); // no-op for the null field
  for (int i=0; i<app->num_species; ++i)
    vlasov_species_write(app, &app->species[i], tm, frame);
}

void
gkyl_vlasov_app_write_field(gkyl_vlasov_app* app, double tm, int frame)
{
  // The distributions let Vlasov-Poisson solve for the potential at time tm.
  const struct gkyl_array *fin[app->num_species];
  vlasov_species_gather_dist(app, fin);
  vlasov_field_write(app, tm, frame, fin);
}

void
gkyl_vlasov_app_write_species(gkyl_vlasov_app* app, int sidx, double tm, int frame)
{
  vlasov_species_write(app, &app->species[sidx], tm, frame);
}

void
gkyl_vlasov_app_write_mom(gkyl_vlasov_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i)
    vlasov_species_write_mom(app, &app->species[i], tm, frame);
}

void
gkyl_vlasov_app_write_integrated_mom(gkyl_vlasov_app *app)
{
  for (int i=0; i<app->num_species; ++i)
    vlasov_species_write_integrated_mom(app, &app->species[i]);
}

void
gkyl_vlasov_app_write_integrated_L2_f(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    vlasov_species_write_integrated_L2_f(app, &app->species[i]);
}

void
gkyl_vlasov_app_write_field_energy(gkyl_vlasov_app* app)
{
  vlasov_field_write_energy(app); // no-op for the null field
}

void
gkyl_vlasov_app_write_lte_corr_status(gkyl_vlasov_app* app)
{
  for (int i=0; i<app->num_species; ++i)
    vlasov_species_write_lte_corr_status(app, &app->species[i]);
}

struct gkyl_update_status
gkyl_vlasov_update(gkyl_vlasov_app* app, double dt)
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

struct gkyl_vlasov_stat
gkyl_vlasov_app_stat(gkyl_vlasov_app* app)
{
  vm_species_n_iter_corr(app);
  return app->stat;
}

static void
range_stat_write(gkyl_vlasov_app* app, const char *nm, const struct gkyl_range *r, FILE *fp)
{
  gkyl_vlasov_app_cout(app, fp, " %s_cells : [ ", nm);
  for (int i=0; i<r->ndim; ++i)
    gkyl_vlasov_app_cout(app, fp, " %d, ", gkyl_range_shape(r, i));
  gkyl_vlasov_app_cout(app, fp, " ],\n");
}

// ensure stats across processors are made consistent
static void
comm_reduce_app_stat(const gkyl_vlasov_app* app,
  const struct gkyl_vlasov_stat *local, struct gkyl_vlasov_stat *global)
{
  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz == 1) {
    memcpy(global, local, sizeof(struct gkyl_vlasov_stat));
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

  int64_t l_red_n_iter_corr[app->num_species];
  int64_t l_red_num_corr[app->num_species];
  for (int s=0; s<app->num_species; ++s) {
    l_red_n_iter_corr[s] = local->n_iter_corr[s];
    l_red_num_corr[s] = local->num_corr[s];
  }

  int64_t l_red_global_n_iter_corr[app->num_species];
  int64_t l_red_global_num_corr[app->num_species];
  gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, app->num_species, 
    l_red_n_iter_corr, l_red_global_n_iter_corr);
  gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, app->num_species, 
    l_red_num_corr, l_red_global_num_corr);

  for (int s=0; s<app->num_species; ++s) {
    global->n_iter_corr[s] = l_red_global_n_iter_corr[s];
    global->num_corr[s] = l_red_global_num_corr[s];
  }

  enum {
    TOTAL_TM, RK3_TM, FL_EM_TM, 
    INIT_SPECIES_TM, INIT_FLUID_SPECIES_TM, INIT_FIELD_TM, 
    SPECIES_RHS_TM, FLUID_SPECIES_RHS_TM, FLUID_SPECIES_VARS_TM, 
    SPECIES_COLL_MOM_TM, SPECIES_COL_TM, SPECIES_RAD_TM, SPECIES_LTE_TM, 
    FIELD_RHS_TM, CURRENT_TM,
    SPECIES_OMEGA_CFL_TM, FIELD_OMEGA_CFL_TM, 
    SPECIES_BC_TM, FLUID_SPECIES_BC_TM, FIELD_BC_TM,
    SPECIES_IO_TM, SPECIES_DIAG_CALC_TM, SPECIES_DIAG_IO_TM, 
    FIELD_IO_TM, FIELD_DIAG_CALC_TM, FIELD_DIAG_IO_TM, 
    D_END
  };

  double d_red[D_END] = {
    [TOTAL_TM] = local->total_tm,
    [RK3_TM] = local->rk3_tm,
    [FL_EM_TM] = local->fl_em_tm,
    [INIT_SPECIES_TM] = local->init_species_tm,
    [INIT_FLUID_SPECIES_TM] = local->init_fluid_species_tm,
    [INIT_FIELD_TM] = local->init_field_tm,
    [SPECIES_RHS_TM] = local->species_rhs_tm,
    [FLUID_SPECIES_RHS_TM] = local->fluid_species_rhs_tm,
    [FLUID_SPECIES_VARS_TM] = local->fluid_species_vars_tm,
    [SPECIES_COLL_MOM_TM] = local->species_coll_mom_tm,
    [SPECIES_COL_TM] = local->species_coll_tm,
    [SPECIES_RAD_TM] = local->species_rad_tm,
    [SPECIES_LTE_TM] = local->species_lte_tm,
    [FIELD_RHS_TM] = local->field_rhs_tm,
    [CURRENT_TM] = local->current_tm,
    [SPECIES_OMEGA_CFL_TM] = local->species_omega_cfl_tm,
    [FIELD_OMEGA_CFL_TM] = local->field_omega_cfl_tm,
    [SPECIES_BC_TM] = local->species_bc_tm,
    [FLUID_SPECIES_BC_TM] = local->fluid_species_bc_tm,
    [FIELD_BC_TM] = local->field_bc_tm, 
    [SPECIES_IO_TM] = local->species_io_tm,
    [SPECIES_DIAG_CALC_TM] = local->species_diag_calc_tm,
    [SPECIES_DIAG_IO_TM] = local->species_diag_io_tm,
    [FIELD_IO_TM] = local->field_io_tm,
    [FIELD_DIAG_CALC_TM] = local->field_diag_calc_tm,
    [FIELD_DIAG_IO_TM] = local->field_diag_io_tm,
  };

  double d_red_global[D_END];
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, D_END, d_red, d_red_global);
  
  global->total_tm = d_red_global[TOTAL_TM];
  global->rk3_tm = d_red_global[RK3_TM];
  global->fl_em_tm = d_red_global[FL_EM_TM];
  global->init_species_tm = d_red_global[INIT_SPECIES_TM];
  global->init_fluid_species_tm = d_red_global[INIT_FLUID_SPECIES_TM];
  global->init_field_tm = d_red_global[INIT_FIELD_TM];
  global->species_rhs_tm = d_red_global[SPECIES_RHS_TM];
  global->fluid_species_rhs_tm = d_red_global[FLUID_SPECIES_RHS_TM];
  global->fluid_species_vars_tm = d_red_global[FLUID_SPECIES_VARS_TM];
  global->species_coll_mom_tm = d_red_global[SPECIES_COLL_MOM_TM];
  global->species_coll_tm = d_red_global[SPECIES_COL_TM];
  global->species_rad_tm = d_red_global[SPECIES_RAD_TM];
  global->species_lte_tm = d_red_global[SPECIES_LTE_TM];
  global->field_rhs_tm = d_red_global[FIELD_RHS_TM];
  global->current_tm = d_red_global[CURRENT_TM];
  global->species_omega_cfl_tm = d_red_global[SPECIES_OMEGA_CFL_TM];
  global->field_omega_cfl_tm = d_red_global[FIELD_OMEGA_CFL_TM];
  global->species_bc_tm = d_red_global[SPECIES_BC_TM];
  global->fluid_species_bc_tm = d_red_global[FLUID_SPECIES_BC_TM];
  global->field_bc_tm = d_red_global[FIELD_BC_TM];
  global->species_io_tm = d_red_global[SPECIES_IO_TM];
  global->species_diag_calc_tm = d_red_global[SPECIES_DIAG_CALC_TM];
  global->species_diag_io_tm = d_red_global[SPECIES_DIAG_IO_TM];
  global->field_io_tm = d_red_global[FIELD_IO_TM];
  global->field_diag_calc_tm = d_red_global[FIELD_DIAG_CALC_TM];
  global->field_diag_io_tm = d_red_global[FIELD_DIAG_IO_TM];

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
gkyl_vlasov_app_stat_write(gkyl_vlasov_app* app)
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

  vm_species_n_iter_corr(app);

  struct gkyl_vlasov_stat stat = { };
  comm_reduce_app_stat(app, &app->stat, &stat);
  
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  if (rank == 0) fp = fopen(fileNm, "a");

  gkyl_vlasov_app_cout(app, fp, "{\n");

  if (strftime(buff, sizeof buff, "%c", &curr_tm))
    gkyl_vlasov_app_cout(app, fp, " date : %s,\n", buff);

  gkyl_vlasov_app_cout(app, fp, " use_gpu : %d,\n", stat.use_gpu);
  gkyl_vlasov_app_cout(app, fp, " num_ranks : %d,\n", num_ranks); 
  
  for (int s=0; s<app->num_species; ++s)
    if (app->species[s].kinetic)
      range_stat_write(app, app->species[s].kinetic->name, &app->species[s].kinetic->global, fp);
  
  gkyl_vlasov_app_cout(app, fp, " nup : %ld,\n", stat.nup);
  gkyl_vlasov_app_cout(app, fp, " nfeuler : %ld,\n", stat.nfeuler);
  gkyl_vlasov_app_cout(app, fp, " nstage_2_fail : %ld,\n", stat.nstage_2_fail);
  gkyl_vlasov_app_cout(app, fp, " nstage_3_fail : %ld,\n", stat.nstage_3_fail);

  gkyl_vlasov_app_cout(app, fp, " stage_2_dt_diff : [ %lg, %lg ],\n",
    stat.stage_2_dt_diff[0], stat.stage_2_dt_diff[1]);
  gkyl_vlasov_app_cout(app, fp, " stage_3_dt_diff : [ %lg, %lg ],\n",
    stat.stage_3_dt_diff[0], stat.stage_3_dt_diff[1]);

  gkyl_vlasov_app_cout(app, fp, " total_tm : %lg,\n", stat.total_tm);
  gkyl_vlasov_app_cout(app, fp, " rk3_tm : %lg,\n", stat.rk3_tm);
  gkyl_vlasov_app_cout(app, fp, " fluid_em_coupling_tm : %lg,\n", stat.fl_em_tm);
  gkyl_vlasov_app_cout(app, fp, " init_species_tm : %lg,\n", stat.init_species_tm);
  if (app->has_field)
    gkyl_vlasov_app_cout(app, fp, " init_field_tm : %lg,\n", stat.init_field_tm);
  
  gkyl_vlasov_app_cout(app, fp, " species_rhs_tm : %lg,\n", stat.species_rhs_tm);

  for (int s=0; s<app->num_species; ++s) {
    gkyl_vlasov_app_cout(app, fp, " n_iter_corr[%d] : %ld,\n", s, 
      stat.n_iter_corr[s]);
    gkyl_vlasov_app_cout(app, fp, " num_corr[%d] : %ld,\n", s, 
      stat.num_corr[s]);   
    gkyl_vlasov_app_cout(app, fp, " species_coll_drag_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_drag_tm[s]);
    gkyl_vlasov_app_cout(app, fp, " species_coll_diff_tm[%d] : %lg,\n", s,
      stat.species_lbo_coll_diff_tm[s]);
  }

  gkyl_vlasov_app_cout(app, fp, " species_coll_mom_tm : %lg,\n", stat.species_coll_mom_tm);
  gkyl_vlasov_app_cout(app, fp, " species_coll_tm : %lg,\n", stat.species_coll_tm);
  gkyl_vlasov_app_cout(app, fp, " species_rad_tm : %lg,\n", stat.species_rad_tm);

  gkyl_vlasov_app_cout(app, fp, " species_lte_tm : %lg,\n", stat.species_lte_tm);

  gkyl_vlasov_app_cout(app, fp, " species_bc_tm : %lg,\n", stat.species_bc_tm);
  
  gkyl_vlasov_app_cout(app, fp, " fluid_species_rhs_tm : %lg,\n", stat.fluid_species_rhs_tm);

  gkyl_vlasov_app_cout(app, fp, " fluid_species_bc_tm : %lg,\n", stat.fluid_species_bc_tm);

  if (app->has_field) {
    gkyl_vlasov_app_cout(app, fp, " field_rhs_tm : %lg,\n", stat.field_rhs_tm);
    gkyl_vlasov_app_cout(app, fp, " field_bc_tm : %lg,\n", stat.field_bc_tm);
    
    gkyl_vlasov_app_cout(app, fp, " current_tm : %lg,\n", stat.current_tm);
  }

  gkyl_vlasov_app_cout(app, fp, " species_omega_cfl_tm : %lg\n", stat.species_omega_cfl_tm);
  gkyl_vlasov_app_cout(app, fp, " field_omega_cfl_tm : %lg\n", stat.field_omega_cfl_tm);  
  gkyl_vlasov_app_cout(app, fp, " species_io_tm : %lg\n", stat.species_io_tm);
  gkyl_vlasov_app_cout(app, fp, " species_diag_tm : %lg\n", stat.species_diag_calc_tm);
  gkyl_vlasov_app_cout(app, fp, " species_diag_io_tm : %lg\n", stat.species_diag_io_tm);
  gkyl_vlasov_app_cout(app, fp, " field_io_tm : %lg\n", stat.field_io_tm);
  gkyl_vlasov_app_cout(app, fp, " field_diag_tm : %lg\n", stat.field_diag_calc_tm);
  gkyl_vlasov_app_cout(app, fp, " field_diag_io_tm : %lg\n", stat.field_diag_io_tm);

  gkyl_vlasov_app_cout(app, fp, " n_species_omega_cfl : %ld,\n", stat.n_species_omega_cfl);
  gkyl_vlasov_app_cout(app, fp, " n_field_omega_cfl : %ld,\n", stat.n_field_omega_cfl);
  gkyl_vlasov_app_cout(app, fp, " n_mom : %ld,\n", stat.n_mom);
  gkyl_vlasov_app_cout(app, fp, " n_diag : %ld,\n", stat.n_diag);
  gkyl_vlasov_app_cout(app, fp, " n_io : %ld,\n", stat.n_io);
  gkyl_vlasov_app_cout(app, fp, " n_diag_io : %ld,\n", stat.n_diag_io);  
  gkyl_vlasov_app_cout(app, fp, " n_field_diag : %ld,\n", stat.n_field_diag);
  gkyl_vlasov_app_cout(app, fp, " n_field_io : %ld,\n", stat.n_field_io);
  gkyl_vlasov_app_cout(app, fp, " n_field_diag_io : %ld,\n", stat.n_field_diag_io);  

  gkyl_vlasov_app_cout(app, fp, "}\n");

  if (rank == 0)
    fclose(fp);  

}

struct gkyl_app_restart_status
vlasov_header_from_file(gkyl_vlasov_app *app, const char *fname)
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

    struct vlasov_output_meta meta =
      vlasov_meta_from_mpack( &(struct gkyl_msgpack_data) {
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
gkyl_vlasov_app_from_file_field(gkyl_vlasov_app *app, const char *fname)
{
  return vlasov_field_from_file(app, fname);
}

struct gkyl_app_restart_status 
gkyl_vlasov_app_from_file_species(gkyl_vlasov_app *app, int sidx,
  const char *fname)
{
  return vlasov_species_from_file(app, &app->species[sidx], fname);
}

struct gkyl_app_restart_status
gkyl_vlasov_app_from_frame_field(gkyl_vlasov_app *app, int frame)
{
  return vlasov_field_read_from_frame(app, frame);
}

struct gkyl_app_restart_status
gkyl_vlasov_app_from_frame_species(gkyl_vlasov_app *app, int sidx, int frame)
{
  return vlasov_species_read_from_frame(app, &app->species[sidx], frame);
}

struct gkyl_app_restart_status
gkyl_vlasov_app_read_from_frame(gkyl_vlasov_app *app, int frame)
{
  struct gkyl_app_restart_status rstat;

  // Field and species restarts are independent: the Poisson potential is solved
  // from the distribution whenever it is needed, so no re-solve happens here.
  rstat = gkyl_vlasov_app_from_frame_field(app, frame);

  for (int i = 0; i < app->num_species; i++) {
    rstat = vlasov_species_read_from_frame(app, &app->species[i], frame);
  }

  return rstat;
}

// private function to handle variable argument list for printing
static void
v_vlasov_app_cout(const gkyl_vlasov_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp) {
    vfprintf(fp, fmt, argp);
    fflush(fp);
  }
}

void
gkyl_vlasov_app_cout(const gkyl_vlasov_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_vlasov_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_vlasov_app_release(gkyl_vlasov_app* app)
{
  vm_geom_release(app, app->vm_geom);
  gkyl_free(app->vm_geom);
  for (int i=0; i<app->num_species; ++i)
    vlasov_species_release(app, &app->species[i]);
  // 'fluid_species' is a view into 'species', which owns the backing array.
  if (app->species)
    gkyl_free(app->species);
  // A field object always exists (the null field when no field is present).
  vlasov_field_release(app);
  if (app->has_fluid_em_coupling)
    vm_fluid_em_coupling_release(app, app->fl_em);

  gkyl_comm_release(app->comm);
  gkyl_rect_decomp_release(app->decomp);

  gkyl_wave_geom_release(app->geom);
  gkyl_vlasov_position_map_release(app->pos_map);

  if (app->use_gpu) {
    gkyl_cu_free(app->basis_on_dev);
  }

  gkyl_free(app);
}
