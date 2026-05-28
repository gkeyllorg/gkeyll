#include <stdarg.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_basis.h>
#include <gkyl_comm_io.h>
#include <gkyl_dflt.h>
#include <gkyl_dynvec.h>
#include <gkyl_null_comm.h>
#include <gkyl_nodal_ops.h>

#include <gkyl_gyrokinetic_priv.h>
#include <gkyl_app_priv.h>

#include <mpack.h>

void
gyrokinetic_cuts_check(struct gkyl_gyrokinetic_app* app, struct gkyl_comm *comm, const int *cuts, FILE *iostream)
{
  // A temporary function that checks the consistency of the communicator and
  // cuts provided for a simulation (e.g at the moment we only decompose along
  // z, and we need to check that cuts meets that requirement).
  int cdim = app->cdim;

  // Create decomposition.
  int cuts_used[cdim];
#ifdef GKYL_HAVE_MPI
  for (int d = 0; d < cdim; d++)
    cuts_used[d] = cuts[d];
#else
  for (int d = 0; d < cdim; d++) cuts_used[d] = 1;
#endif

  int comm_rank, comm_size;
  gkyl_comm_get_rank(comm, &comm_rank);
  gkyl_comm_get_size(comm, &comm_size);

  int ncuts = 1;
  for (int d = 0; d < cdim; d++) ncuts *= cuts_used[d];

  if (ncuts != comm_size) {
    if (comm_rank == 0)
      fprintf(iostream, "\n*** Number of ranks, %d, does not match total cuts, %d!\n\n", comm_size, ncuts);
    assert(false);
  }

  for (int d = 0; d < cdim - 1; d++) {
    if (cuts_used[d] > 1) {
      if (comm_rank == 0)
        fprintf(iostream,
          "\n*** Parallelization only allowed in z. Number of ranks, %d, in direction %d cannot be > 1!\n\n", cuts_used[d], d);
      assert(false);
    }
  }

  if (app->grid.cells[cdim - 1] % cuts_used[cdim - 1] != 0) {
    if (comm_rank == 0)
      fprintf(iostream,
        "\n*** Number of cells in z, %d, not divisible by number of cuts, %d!\n\n",
        app->grid.cells[cdim - 1], cuts_used[cdim - 1]);
    assert(false);
  }
}

static bool gyrokinetic_str_ends_in_b0(char *name){
  size_t len = strlen(name);
  int i = len - 1;
  int digit_count = 0;
  while (i >= 0 && isdigit((unsigned char)name[i])) {
    i--;
    digit_count++;
  }
  if (digit_count > 0 && i >= 1 && name[i] == 'b' && name[i-1] == '_') {
    const char *num_str = &name[i + 1];
    int num = atoi(num_str);
    if ( num == 0)
      return true;
    else
      return false;
  }
  else {
    return true;
  }
}

gkyl_gyrokinetic_app*
gkyl_gyrokinetic_app_new_geom(struct gkyl_gk *gk)
{
  disable_denorm_float();

  assert(gk->num_species <= GKYL_MAX_SPECIES);

  gkyl_gyrokinetic_app *app = gkyl_malloc(sizeof(gkyl_gyrokinetic_app));

  int cdim = app->cdim = gk->cdim;
  int poly_order = app->poly_order = gk->poly_order;
  int ns = app->num_species = gk->num_species;
  int neuts = app->num_neut_species = gk->num_neut_species;

  double cfl_frac = gk->cfl_frac == 0 ? 1.0 : gk->cfl_frac;
  app->cfl = cfl_frac;

  // The value 1.7 here is based on figure 2.4a in Durran's "Numerical methods
  // for fluid dynamics" textbook for a purely oscillatory mode and RK3.
  double cfl_frac_omegaH = fabs(gk->cfl_frac_omegaH) < 1e-16 ? 1.7 : gk->cfl_frac_omegaH;
  app->cfl_omegaH = cfl_frac_omegaH;

#ifdef GKYL_HAVE_CUDA
  app->use_gpu = gk->parallelism.use_gpu;
#else
  app->use_gpu = false; // can't use GPUs if we don't have them!
#endif

  app->num_periodic_dir = gk->num_periodic_dir;
  for (int d=0; d<cdim; ++d)
    app->periodic_dirs[d] = gk->periodic_dirs[d];

  strcpy(app->name, gk->name);
  app->tcurr = 0.0; // reset on init

  if (app->use_gpu) {
    // allocate device basis if we are using GPUs
    app->basis_on_dev = gkyl_cu_malloc(sizeof(struct gkyl_basis));
  }
  else {
    app->basis_on_dev = &app->basis;
  }

  // basis functions
  switch (gk->basis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      gkyl_cart_modal_serendip(&app->basis, cdim, poly_order);
      if (app->use_gpu) {
        gkyl_cart_modal_serendip_cu_dev(app->basis_on_dev, cdim, poly_order);
      }
      break;
    default:
      assert(false);
      break;
  }

  gkyl_rect_grid_init(&app->grid, cdim, gk->lower, gk->upper, gk->cells);

  int ghost[] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&app->grid, ghost, &app->global_ext, &app->global);

  if (gk->parallelism.comm == 0) {
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
    gyrokinetic_cuts_check(app, gk->parallelism.comm, gk->parallelism.cuts, stdout);

    // Create decomp.
    app->decomp = gkyl_rect_decomp_new_from_cuts(app->cdim, gk->parallelism.cuts, &app->global);

    // Create a new communicator with the decomposition in it.
    app->comm = gkyl_comm_split_comm(gk->parallelism.comm, 0, app->decomp);

    // Create local and local_ext.
    int rank;
    gkyl_comm_get_rank(app->comm, &rank);
    gkyl_create_ranges(&app->decomp->ranges[rank], ghost, &app->local_ext, &app->local);
  }

  // Create plane communicators.
  if (app->cdim == 1) {
    app->decomp_plane[0] = gkyl_rect_decomp_acquire(app->decomp);
    app->comm_plane[0] = gkyl_comm_acquire(app->comm);
  }
  else {
    for (int dir=0; dir<app->cdim; ++dir) {
      // Identify ranks on the same plane as this one.
      int num_ranks_plane = 0;
      int ranks_plane[app->decomp->ndecomp]; 
      for (int i=0; i<app->decomp->ndecomp; i++) {
        if (app->decomp->ranges[i].lower[dir] == app->local.lower[dir]) {
          ranks_plane[num_ranks_plane] = i;
          num_ranks_plane++;
        }
      }
      // Create a range tangentially global, and local in perp direction.
      int lower_plane[app->cdim], upper_plane[app->cdim];
      for (int d=0; d<app->cdim; ++d) {
        lower_plane[d] = app->global.lower[d];
        upper_plane[d] = app->global.upper[d];
      }
      lower_plane[dir] = app->local.lower[dir];
      upper_plane[dir] = app->local.upper[dir];
      struct gkyl_range range_plane;
      gkyl_range_init(&range_plane, app->cdim, lower_plane, upper_plane);
  
      // Create decomp.
      int cuts_plane[GKYL_MAX_CDIM];
      for (int d=0; d<app->cdim; ++d)
        cuts_plane[d] = gk->parallelism.cuts[d];
      cuts_plane[dir] = 1;
      app->decomp_plane[dir] = gkyl_rect_decomp_new_from_cuts(app->cdim, cuts_plane, &range_plane);
  
      // Create a new communicator with ranks on plane.
      bool is_comm_valid;
      app->comm_plane[dir] = gkyl_comm_create_comm_from_ranks(app->comm, num_ranks_plane,
        ranks_plane, app->decomp_plane[dir], &is_comm_valid);
      assert(is_comm_valid);
    }
  }

  // Skin and ghost ranges for configuration space fields.
  for (int dir=0; dir<cdim; ++dir) {
    gkyl_skin_ghost_ranges(&app->local_lower_skin[dir], &app->local_lower_ghost[dir],
      dir, GKYL_LOWER_EDGE, &app->local_ext, ghost); 
    gkyl_skin_ghost_ranges(&app->local_upper_skin[dir], &app->local_upper_ghost[dir],
      dir, GKYL_UPPER_EDGE, &app->local_ext, ghost);
    gkyl_skin_ghost_ranges(&app->global_lower_skin[dir], &app->global_lower_ghost[dir],
      dir, GKYL_LOWER_EDGE, &app->global_ext, ghost); 
    gkyl_skin_ghost_ranges(&app->global_upper_skin[dir], &app->global_upper_ghost[dir],
      dir, GKYL_UPPER_EDGE, &app->global_ext, ghost);
  }

  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);

  // Configuration space geometry initialization
  app->position_map = gkyl_position_map_new(gk->geometry.position_map_info, app->grid, app->local, 
      app->local_ext, app->global, app->global_ext, app->basis);

  // Initialize the input struct from user side input struct
  struct gkyl_gk_geometry_inp geometry_inp = {
    .geometry_id  = gk->geometry.geometry_id,
    .c2p_ctx = gk->geometry.c2p_ctx,
    .mapc2p = gk->geometry.mapc2p,
    .bfield_ctx = gk->geometry.bfield_ctx,
    .bfield_func = gk->geometry.bfield_func,
    .efit_info = gk->geometry.efit_info,
    .tok_grid_info = gk->geometry.tok_grid_info,
    .mirror_grid_info = gk->geometry.mirror_grid_info,
    .position_map = app->position_map,
    .grid = app->grid,
    .local = app->local,
    .local_ext = app->local_ext,
    .global = app->global,
    .global_ext = app->global_ext,
    .basis = app->basis,
    .comm = app->comm,
    .has_LCFS = gk->geometry.has_LCFS,
    .x_LCFS = gk->geometry.x_LCFS,
  };
  strcpy(geometry_inp.geometry_path, gk->geometry.geometry_path);
  for(int i = 0; i<3; i++)
    geometry_inp.world[i] = gk->geometry.world[i];

  if (app->cdim < 3){
    geometry_inp.geo_grid = gkyl_gk_geometry_augment_grid(app->grid, geometry_inp);
    switch (gk->basis_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        gkyl_cart_modal_serendip(&geometry_inp.geo_basis, 3, poly_order);
        break;
      default:
        assert(false);
        break;
    }

    int ghost[] = { 1, 1, 1 };
    gkyl_create_grid_ranges(&geometry_inp.geo_grid, ghost, &geometry_inp.geo_global_ext, &geometry_inp.geo_global);
    if (comm_sz > 1) {
      // Create local and local_ext from user-supplied local range.
      gkyl_gk_geometry_augment_local(&app->local, ghost, &geometry_inp.geo_local_ext, &geometry_inp.geo_local);
    }
    else {
      // Global and local ranges are same, and so just copy.
      memcpy(&geometry_inp.geo_local, &geometry_inp.geo_global, sizeof(struct gkyl_range));
      memcpy(&geometry_inp.geo_local_ext, &geometry_inp.geo_global_ext, sizeof(struct gkyl_range));
    }
  }
  else {
    geometry_inp.geo_grid = app->grid;
    geometry_inp.geo_local = app->local;
    geometry_inp.geo_local_ext = app->local_ext;
    geometry_inp.geo_global = app->global;
    geometry_inp.geo_global_ext = app->global_ext;
    geometry_inp.geo_basis = app->basis;
  }

  struct gk_geometry* gk_geom_3d;
  switch (geometry_inp.geometry_id) {
    case GKYL_GEOMETRY_FROMFILE:
      gk_geom_3d = gkyl_gk_geometry_new(app->gk_geom, &geometry_inp, false);
      break;
    case GKYL_GEOMETRY_TOKAMAK:
      gk_geom_3d = gkyl_gk_geometry_tok_new(&geometry_inp);
      break;
    case GKYL_GEOMETRY_MIRROR:
      gk_geom_3d = gkyl_gk_geometry_mirror_new(&geometry_inp);
      break;
    case GKYL_GEOMETRY_MAPC2P:
    case GKYL_GEOMETRY_NONE:
      gk_geom_3d = gkyl_gk_geometry_mapc2p_new(&geometry_inp);
      break;
  }

  if (geometry_inp.geometry_id != GKYL_GEOMETRY_FROMFILE) {
    // Deflate geometry if necessary.
    if (app->cdim < 3)
      app->gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, &geometry_inp);
    else
      app->gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);
  }
  else {
    app->gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);
    gkyl_gyrokinetic_app_read_geometry(app, &geometry_inp);
    gkyl_gk_geometry_populate_nodal(app->gk_geom);
  }

  gkyl_gk_geometry_release(gk_geom_3d); // Release temporary 3d geometry.

  double bmag_min_local, bmag_min_global;
  bmag_min_local = gkyl_gk_geometry_reduce_bmag(app->gk_geom, GKYL_MIN);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &bmag_min_local, &bmag_min_global);

  double bmag_max_local, bmag_max_global;
  bmag_max_local = gkyl_gk_geometry_reduce_bmag(app->gk_geom, GKYL_MAX);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 1, &bmag_max_local, &bmag_max_global);

  app->bmag_ref = (bmag_max_global + bmag_min_global)/2.0;

  gkyl_position_map_set_mc2nu(app->position_map, app->gk_geom->geo_corn.mc2nu_pos);

  const struct gkyl_dg_geom_inp dg_geom_inp = {
    .grid = &app->grid,
    .range = &app->local_ext,
    .nquad = 2,
  };

  const struct gkyl_gk_dg_geom_inp gk_dg_geom_inp = {
    .grid = &app->grid,
    .range = &app->local_ext,
    .nquad = 2,
  };

  // Populate dg_geom and gk_dg_geom objects with geometric quantities
  app->dg_geom = gkyl_dg_geom_new(&dg_geom_inp);
  app->gk_dg_geom = gkyl_gk_dg_geom_new(&gk_dg_geom_inp);

  gkyl_gk_dg_geom_populate_vol(app->dg_geom, app->gk_dg_geom, app->gk_geom);
  gkyl_gk_dg_geom_populate_surf(app->dg_geom, app->gk_dg_geom, app->gk_geom);

  // If we are on the gpu, copy from host.
  if (app->use_gpu) {
    struct gkyl_dg_geom* dg_geom_dev = gkyl_dg_geom_new_from_host(&dg_geom_inp, app->dg_geom, true);
    struct gkyl_gk_dg_geom* gk_dg_geom_dev = gkyl_gk_dg_geom_new_from_host(&gk_dg_geom_inp, app->gk_dg_geom, true);
    struct gk_geometry* gk_geom_dev = gkyl_gk_geometry_new(app->gk_geom, &geometry_inp, app->use_gpu);

    gkyl_gk_geometry_release(app->gk_geom);
    app->gk_geom = gkyl_gk_geometry_acquire(gk_geom_dev);
    gkyl_gk_geometry_release(gk_geom_dev);

    gkyl_dg_geom_release(app->dg_geom);
    app->dg_geom = gkyl_dg_geom_acquire(dg_geom_dev);
    gkyl_dg_geom_release(dg_geom_dev);

    gkyl_gk_dg_geom_release(app->gk_dg_geom);
    app->gk_dg_geom = gkyl_gk_dg_geom_acquire(gk_dg_geom_dev);
    gkyl_gk_dg_geom_release(gk_dg_geom_dev);
  }

  // Basic metadata for I/O.
  const char* build_id = GIT_COMMIT_ID;
  const char* build_date = GKYL_BUILD_DATE;
  struct gkyl_msgpack_map_elem io_meta_basic[] = {
   { .key = "changeset", .elem_type = GKYL_MP_STRING, .cval = (char *)build_id },
   { .key = "builddate", .elem_type = GKYL_MP_STRING, .cval = (char *)build_date },
   { .key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = 0.0 },
   { .key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
  };
  app->io_meta_basic_len = sizeof(io_meta_basic)/sizeof(io_meta_basic[0]);
  app->io_meta_basic = gkyl_msgpack_map_elem_clone(app->io_meta_basic_len, io_meta_basic);
  // Metadata for GK app.
  struct gkyl_msgpack_map_elem io_meta[] = {
    { .key = "poly_order", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = app->basis.poly_order },
    { .key = "basis_type", .elem_type = GKYL_MP_STRING, .cval = app->basis.id }
  };
  app->io_meta_len = sizeof(io_meta)/sizeof(io_meta[0]);
  app->io_meta = gkyl_msgpack_map_elem_clone(app->io_meta_len, io_meta);

  gkyl_gyrokinetic_app_write_geometry(app, &geometry_inp);

  // Allocate 1/(J.B) using weak mul/div.
  struct gkyl_array *tmp = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  app->jacobtot_inv_weak = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  gkyl_dg_mul_op_range(app->basis, 0, tmp, 0, app->gk_geom->geo_int.bmag, 0, app->gk_geom->geo_int.jacobgeo, &app->local); 
  gkyl_dg_inv_op_range(app->basis, 0, app->jacobtot_inv_weak, 0, tmp, &app->local); 
  gkyl_array_release(tmp);

  if (gk->geometry.has_LCFS) {
    // Simulation spans the last-closed flux surface (LCFS). Create core and SOL global ranges.
    int idx_LCFS_lo = app->gk_geom->idx_LCFS_lo;
    // Length of lower and upper x ranges (one is core, the other SOL).
    int len_lo = idx_LCFS_lo;
    int len_up = app->global.upper[0]-len_lo;
    // Lower and upper x ranges.
    struct gkyl_range *global_lo_r, *global_up_r, *global_ext_lo_r, *global_ext_up_r;
    struct gkyl_range *global_lower_skin_par_lo_r , *global_upper_skin_par_lo_r ,
                      *global_lower_ghost_par_lo_r, *global_upper_ghost_par_lo_r;
    struct gkyl_range *global_lower_skin_par_up_r , *global_upper_skin_par_up_r ,
                      *global_lower_ghost_par_up_r, *global_upper_ghost_par_up_r;
    struct gkyl_range *local_lo_r, *local_up_r, *local_ext_lo_r, *local_ext_up_r;
    struct gkyl_range *local_lower_skin_par_lo_r , *local_upper_skin_par_lo_r ,
                      *local_lower_ghost_par_lo_r, *local_upper_ghost_par_lo_r;
    struct gkyl_range *local_lower_skin_par_up_r , *local_upper_skin_par_up_r ,
                      *local_lower_ghost_par_up_r, *local_upper_ghost_par_up_r;
    if (app->gk_geom->geqdsk_sign_convention == 0) {
      // x increases towards SOL.
      global_lo_r                 = &app->global_core;
      global_up_r                 = &app->global_sol;
      global_ext_lo_r             = &app->global_ext_core;
      global_ext_up_r             = &app->global_ext_sol;
      global_lower_skin_par_lo_r  = &app->global_lower_skin_par_core;
      global_lower_ghost_par_lo_r = &app->global_lower_ghost_par_core;
      global_lower_skin_par_up_r  = &app->global_lower_skin_par_sol;
      global_lower_ghost_par_up_r = &app->global_lower_ghost_par_sol;
      global_upper_skin_par_lo_r  = &app->global_upper_skin_par_core;
      global_upper_ghost_par_lo_r = &app->global_upper_ghost_par_core;
      global_upper_skin_par_up_r  = &app->global_upper_skin_par_sol;
      global_upper_ghost_par_up_r = &app->global_upper_ghost_par_sol;
      local_lo_r                  = &app->local_core;
      local_up_r                  = &app->local_sol;
      local_ext_lo_r              = &app->local_ext_core;
      local_ext_up_r              = &app->local_ext_sol;
      local_lower_skin_par_lo_r   = &app->local_lower_skin_par_core;
      local_lower_ghost_par_lo_r  = &app->local_lower_ghost_par_core;
      local_lower_skin_par_up_r   = &app->local_lower_skin_par_sol;
      local_lower_ghost_par_up_r  = &app->local_lower_ghost_par_sol;
      local_upper_skin_par_lo_r   = &app->local_upper_skin_par_core;
      local_upper_ghost_par_lo_r  = &app->local_upper_ghost_par_core;
      local_upper_skin_par_up_r   = &app->local_upper_skin_par_sol;
      local_upper_ghost_par_up_r  = &app->local_upper_ghost_par_sol;
    }
    else {
      // x increases towards core.
      global_lo_r                 = &app->global_sol;
      global_up_r                 = &app->global_core;
      global_ext_lo_r             = &app->global_ext_sol;
      global_ext_up_r             = &app->global_ext_core;
      global_lower_skin_par_lo_r  = &app->global_lower_skin_par_sol;
      global_lower_ghost_par_lo_r = &app->global_lower_ghost_par_sol;
      global_lower_skin_par_up_r  = &app->global_lower_skin_par_core;
      global_lower_ghost_par_up_r = &app->global_lower_ghost_par_core;
      global_upper_skin_par_lo_r  = &app->global_upper_skin_par_sol;
      global_upper_ghost_par_lo_r = &app->global_upper_ghost_par_sol;
      global_upper_skin_par_up_r  = &app->global_upper_skin_par_core;
      global_upper_ghost_par_up_r = &app->global_upper_ghost_par_core;
      local_lo_r                  = &app->local_sol;
      local_up_r                  = &app->local_core;
      local_ext_lo_r              = &app->local_ext_sol;
      local_ext_up_r              = &app->local_ext_core;
      local_lower_skin_par_lo_r   = &app->local_lower_skin_par_sol;
      local_lower_ghost_par_lo_r  = &app->local_lower_ghost_par_sol;
      local_lower_skin_par_up_r   = &app->local_lower_skin_par_core;
      local_lower_ghost_par_up_r  = &app->local_lower_ghost_par_core;
      local_upper_skin_par_lo_r   = &app->local_upper_skin_par_sol;
      local_upper_ghost_par_lo_r  = &app->local_upper_ghost_par_sol;
      local_upper_skin_par_up_r   = &app->local_upper_skin_par_core;
      local_upper_ghost_par_up_r  = &app->local_upper_ghost_par_core;
    }

    // Lower and upper x ranges.
    gkyl_range_shorten_from_above(global_lo_r, &app->global, 0, len_lo);
    gkyl_range_shorten_from_below(global_up_r, &app->global, 0, len_up);
    gkyl_range_shorten_from_above(local_lo_r, &app->local, 0, len_lo);
    gkyl_range_shorten_from_below(local_up_r, &app->local, 0, len_up);

    // Extended lower and upper x ranges.
    int len_lo_ext = idx_LCFS_lo+1;
    int len_up_ext = app->global_ext.upper[0]-len_lo;
    gkyl_range_shorten_from_above(global_ext_lo_r, &app->global_ext, 0, len_lo_ext);
    gkyl_range_shorten_from_below(global_ext_up_r, &app->global_ext, 0, len_up_ext);
    gkyl_range_shorten_from_above(local_ext_lo_r, &app->local_ext, 0, len_lo_ext);
    gkyl_range_shorten_from_below(local_ext_up_r, &app->local_ext, 0, len_up_ext);

    // Parallel skin and ghost ranges, limited to the lower and upper x range.
    int par_dir = app->cdim-1;
    for (int e=0; e<2; e++) {
      gkyl_range_shorten_from_above(e==0? global_lower_skin_par_lo_r        : global_upper_skin_par_lo_r,
                                    e==0? &app->global_lower_skin[par_dir]  : &app->global_upper_skin[par_dir], 0, len_lo);
      gkyl_range_shorten_from_above(e==0? global_lower_ghost_par_lo_r       : global_upper_ghost_par_lo_r,
                                    e==0? &app->global_lower_ghost[par_dir] : &app->global_upper_ghost[par_dir], 0, len_lo);
      gkyl_range_shorten_from_below(e==0? global_lower_skin_par_up_r        : global_upper_skin_par_up_r,
                                    e==0? &app->global_lower_skin[par_dir]  : &app->global_upper_skin[par_dir], 0, len_up);
      gkyl_range_shorten_from_below(e==0? global_lower_ghost_par_up_r       : global_upper_ghost_par_up_r,
                                    e==0? &app->global_lower_ghost[par_dir] : &app->global_upper_ghost[par_dir], 0, len_up);
      gkyl_range_shorten_from_above(e==0? local_lower_skin_par_lo_r         : local_upper_skin_par_lo_r,
                                    e==0? &app->local_lower_skin[par_dir]   : &app->local_upper_skin[par_dir], 0, len_lo);
      gkyl_range_shorten_from_above(e==0? local_lower_ghost_par_lo_r        : local_upper_ghost_par_lo_r,
                                    e==0? &app->local_lower_ghost[par_dir]  : &app->local_upper_ghost[par_dir], 0, len_lo);
      gkyl_range_shorten_from_below(e==0? local_lower_skin_par_up_r         : local_upper_skin_par_up_r,
                                    e==0? &app->local_lower_skin[par_dir]   : &app->local_upper_skin[par_dir], 0, len_up);
      gkyl_range_shorten_from_below(e==0? local_lower_ghost_par_up_r        : local_upper_ghost_par_up_r,
                                    e==0? &app->local_lower_ghost[par_dir]  : &app->local_upper_ghost[par_dir], 0, len_up);
    }

    // Core range extended in the parallel direction.
    int ndim = app->cdim;
    int lower_bcdir_ext[ndim], upper_bcdir_ext[ndim];
    for (int i=0; i<ndim; i++) {
      lower_bcdir_ext[i] = app->global_core.lower[i];
      upper_bcdir_ext[i] = app->global_core.upper[i];
    }
    lower_bcdir_ext[par_dir] = app->global_ext_core.lower[par_dir];
    upper_bcdir_ext[par_dir] = app->global_ext_core.upper[par_dir];
    gkyl_sub_range_init(&app->global_par_ext_core, &app->global_ext_core, lower_bcdir_ext, upper_bcdir_ext);

    for (int i=0; i<ndim; i++) {
      lower_bcdir_ext[i] = app->local_core.lower[i];
      upper_bcdir_ext[i] = app->local_core.upper[i];
    }
    lower_bcdir_ext[par_dir] = app->local_ext_core.lower[par_dir];
    upper_bcdir_ext[par_dir] = app->local_ext_core.upper[par_dir];
    gkyl_sub_range_init(&app->local_par_ext_core, &app->local_ext_core, lower_bcdir_ext, upper_bcdir_ext);
  }

  return app;
}

static void
gyrokinetic_calc_field_enabled(gkyl_gyrokinetic_app* app, double tcurr,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  struct timespec wtm = gkyl_wall_clock();
  // Compute electrostatic potential from gyrokinetic Poisson's equation.
  gk_field_accumulate_rho_c(app, app->field, fin, bflux);

  // Compute biased wall potential if present and time-dependent.
  // Note: biased wall potential use eval_on_nodes. 
  // so does copy to GPU every call if app->use_gpu = true.
  if (app->field->phi_wall_lo_evolve || app->field->phi_wall_up_evolve)
    gk_field_calc_phi_wall(app, app->field, tcurr);

  // Solve the field equation.
  gk_field_rhs(app, app->field);
  app->stat.field_tm += gkyl_time_diff_now_sec(wtm);
}

static void
gyrokinetic_calc_field_disabled(gkyl_gyrokinetic_app* app, double tcurr,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  // Do nothing.
}

static void
gkyl_gyrokinetic_app_omegaH_init(gkyl_gyrokinetic_app *app)
{
  // Compute the geometric and field-model dependent part of omega_H.
  // Each species computes its own omega_H as:
  //   omega_H = q_e*sqrt(n_{s0}/m_s) * omegaH_gf
  // where
  //   - n_{s0} is either a reference, average or max density.
  //   - omegaH_gf = (cmag/(jacobgeo*B^_\parallel))*kpar_max / 
  //                 min(sqrt(k_x^2*eps_xx+k_x*k_y*eps_xy+k_y^2*eps_yy+)).
  // and k_x,k_y,k_par are wavenumbers in computational space, and eps_ij is
  // the polarization weight in our field equation.

  app->omegaH_gf = 1.0/DBL_MAX;

  if (!(app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN || app->field->gkfield_id == GKYL_GK_FIELD_ADIABATIC)) {
    // Compute parfac = (cmag/(jacobgeo*B^_\parallel))*kpar_max.
    struct gkyl_array *parfac = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_dg_mul_op_range(app->basis, 0, parfac, 0, app->gk_geom->geo_int.cmag, 0, app->gk_geom->geo_int.jacobtot_inv, &app->local); 
    double kpar_max = M_PI*(app->poly_order+1)/app->grid.dx[app->cdim-1];
    gkyl_array_scale_range(parfac, kpar_max, &app->local);

    // Compute perpfac_inv = 1/sqrt(k_x^2*eps_xx+k_x*k_y*eps_xy+k_y^2*eps_yy+)).
    struct gkyl_array *perpfac = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    struct gkyl_array *perpfac_inv = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    double kx_min = M_PI/(app->grid.upper[0]-app->grid.lower[0]);
    double kx_sq = app->cdim == 1? 1.0 : pow(kx_min,2); // kperp_sq included in epsilon for cdim=1.
    gkyl_array_accumulate_offset_range(perpfac, kx_sq, app->field->epsilon, 0*app->basis.num_basis, &app->local);
    if (app->cdim > 2) {
      double ky_min = M_PI/(app->grid.upper[1]-app->grid.lower[1]);
      gkyl_array_accumulate_offset_range(perpfac, kx_min*ky_min, app->field->epsilon, 1*app->basis.num_basis, &app->local);
      gkyl_array_accumulate_offset_range(perpfac, pow(ky_min,2), app->field->epsilon, 2*app->basis.num_basis, &app->local);
    }
    gkyl_proj_powsqrt_on_basis* proj_sqrt = gkyl_proj_powsqrt_on_basis_new(&app->basis, app->poly_order+1, app->use_gpu);
    gkyl_proj_powsqrt_on_basis_advance(proj_sqrt, &app->local, -1.0, perpfac, perpfac_inv);
    gkyl_proj_powsqrt_on_basis_release(proj_sqrt);

    // Compute max(parfac*perpfac_inv) (using cell centers).
    struct gkyl_array *omegaH_gf_grid = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    gkyl_dg_mul_op_range(app->basis, 0, omegaH_gf_grid, 0, parfac, 0, perpfac_inv, &app->local); 
    double *omegaH_gf_red;
    if (app->use_gpu)
      omegaH_gf_red = gkyl_cu_malloc(app->basis.num_basis*sizeof(double));
    else 
      omegaH_gf_red = gkyl_malloc(app->basis.num_basis*sizeof(double));

    gkyl_array_reduce_range(omegaH_gf_red, omegaH_gf_grid, GKYL_MAX, &app->local);

    if (app->use_gpu)
      gkyl_cu_memcpy(&app->omegaH_gf, omegaH_gf_red, sizeof(double), GKYL_CU_MEMCPY_D2H);
    else
      app->omegaH_gf = omegaH_gf_red[0];
    app->omegaH_gf *= 1.0/pow(sqrt(2.0),app->cdim);

    if (app->use_gpu)
      gkyl_cu_free(omegaH_gf_red);
    else 
      gkyl_free(omegaH_gf_red);
    gkyl_array_release(omegaH_gf_grid);
    gkyl_array_release(perpfac_inv);
    gkyl_array_release(perpfac);
    gkyl_array_release(parfac);
  }
}

static void
gyrokinetic_post_positivity_quasineut_disabled(gkyl_gyrokinetic_app *app, struct gkyl_array *fout[])
{
  // Empty.
}

static void
gyrokinetic_post_positivity_quasineut_enabled(gkyl_gyrokinetic_app *app, struct gkyl_array *fout[])
{
  struct timespec wst = gkyl_wall_clock();
  // Enforce quasineutrality after applying positivity shift to charged species.
  gkyl_array_clear(app->ps_delta_m0_ions, 0.0);
  gkyl_array_clear(app->ps_delta_m0_elcs, 0.0);
  for (int i=0; i<app->num_species; ++i) {
    // Accumulate the shift density of all like-species:
    struct gk_species *gks = &app->species[i];
    struct gk_positivity *pos = &gks->positivity;
    if (pos->quasineut_rescale)
      gkyl_array_accumulate(pos->delta_m0s_tot, 1.0, pos->delta_m0);
  }
  // Rescale each species to enforce quasineutrality.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    struct gk_positivity *pos = &gks->positivity;
    if (pos->quasineut_rescale) {
      gkyl_positivity_shift_gyrokinetic_quasineutrality_scale(pos->shift_op_gk, &app->local, &gks->local,
        pos->delta_m0, pos->delta_m0s_tot, pos->delta_m0r_tot, gks->m0.marr, fout[i]);

      gkyl_array_accumulate(pos->fbuffer_ptr, 1.0, fout[i]);
    }
  }
  app->stat.pos_shift_quasineut_tm += gkyl_time_diff_now_sec(wst);
}

void
gyrokinetic_post_positivity_quasineut_init(gkyl_gyrokinetic_app* app)
{
  int num_quasineut_rescale = gk_species_positivity_num_species_in_quasineut(app);
  if (num_quasineut_rescale > 1) {
    app->post_positivity_quasineut = true;

    // Number density of the positivity shift added over all the ions.
    // Needed before species_init because species store pointers to these.
    app->ps_delta_m0_ions = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
    app->ps_delta_m0_elcs = mkarr(app->use_gpu, app->basis.num_basis, app->local_ext.volume);
  
    app->post_pos_quasineut_func = gyrokinetic_post_positivity_quasineut_enabled;
  }
}

void
gyrokinetic_post_positivity_quasineut(gkyl_gyrokinetic_app *app, struct gkyl_array *fout[])
{
  app->post_pos_quasineut_func(app, fout);
}

void
gyrokinetic_post_positivity_quasineut_release(gkyl_gyrokinetic_app* app)
{
  if (app->post_positivity_quasineut) {
    gkyl_array_release(app->ps_delta_m0_ions);
    gkyl_array_release(app->ps_delta_m0_elcs);
  }
  app->post_positivity_quasineut = false;
  app->post_pos_quasineut_func = gyrokinetic_post_positivity_quasineut_disabled;
}

void
gkyl_gyrokinetic_app_new_solver(struct gkyl_gk *gk, gkyl_gyrokinetic_app *app)
{
  int ns = app->num_species = gk->num_species;
  int neuts = app->num_neut_species = gk->num_neut_species;

  // Allocate space to store species and neutral species objects
  app->species = ns>0 ? gkyl_malloc(sizeof(struct gk_species[ns])) : 0;
  app->neut_species = neuts>0 ? gkyl_malloc(sizeof(struct gk_neut_species[neuts])) : 0;

  // Copy input parameters for each species
  for (int i=0; i<ns; ++i)
    app->species[i].info = gk->species[i];

  for (int i=0; i<neuts; ++i)
    app->neut_species[i].info = gk->neut_species[i];

  app->field = gk_field_new(gk, app); // Initialize field, even if we are  skipping field updates.

  // Choose the function that updates the fields in time.
  if (app->field->update_field)
    app->calc_field_func = gyrokinetic_calc_field_enabled;
  else
    app->calc_field_func = gyrokinetic_calc_field_disabled;

  // Initialize the post-positivity quasineutrality enforcement.
  app->post_positivity_quasineut = false;
  app->post_pos_quasineut_func = gyrokinetic_post_positivity_quasineut_disabled;
  gyrokinetic_post_positivity_quasineut_init(app);

  // Initialize each species.
  for (int i=0; i<ns; ++i)
    gk_species_init(gk, app, &app->species[i]);

  for (int i=0; i<neuts; ++i)
    gk_neut_species_init(gk, app, &app->neut_species[i]);

  // Initialize each species cross-collisions terms.
  for (int i=0; i<ns; ++i) {
    struct gk_species *gk_s = &app->species[i];

    // Initialize cross-species elastic collisions.
    gk_species_lbo_cross_init(app, &app->species[i], &gk_s->lbo);
    gk_species_bgk_cross_init(app, &app->species[i], &gk_s->bgk);

    // Initialize cross-species reactions (ionization, recombination, charge exchange) with charged species.
    gk_species_react_cross_init(app, &app->species[i], &gk_s->react);

    // Initialize cross-species reactions (ionization, recombination, charge exchange) with neutral species.
    gk_species_react_cross_init(app, &app->species[i], &gk_s->react_neut);

    // Initialize line radiation.
    gk_species_radiation_init(app, &app->species[i], &gk_s->rad);

    // Initialize cross-species part of the object that scales the species.
    gk_species_scaling_cross_init(app, gk_s, &gk_s->sca);
  }

  for (int i=0; i<neuts; ++i) {
    // Initialize neutral species cross-species reactions with charged species.
    struct gk_neut_species *gkns = &app->neut_species[i]; 

    gk_neut_species_react_cross_init(app, gkns, &gkns->react_neut);

    // Initialize cross-species part of the object that scales the species.
    gk_neut_species_scaling_cross_init(app, gkns, &gkns->sca);
    
    // Initialize wall emission terms.
    for (int d=0; d<app->cdim; ++d) {
      if (gkns->bc_is_np[d]) {
        if (gkns->lower_bc[d].type == GKYL_BC_GK_SPECIES_RECYCLE)
          gk_neut_species_recycle_cross_init(app, gkns, &gkns->bc_recycle_lo);
        if (gkns->upper_bc[d].type == GKYL_BC_GK_SPECIES_RECYCLE)
          gk_neut_species_recycle_cross_init(app, gkns, &gkns->bc_recycle_up);
      }
    }
  }

  // Initialize source terms. Done here as sources may initialize
  // a boundary flux updater for their source species.
  for (int i=0; i<ns; ++i) {
    gk_species_source_init(app, &app->species[i], &app->species[i].src);
  }
  for (int i=0; i<neuts; ++i) {
    gk_neut_species_source_init(app, &app->neut_species[i], &app->neut_species[i].src);
  }

  // Use implicit BGK collisions if specified
  bool has_implicit_coll_scheme = false;
  for (int i=0; i<ns; ++i){
    if (gk->species[i].collisions.is_implicit){
      has_implicit_coll_scheme = true;
    }
  }
  for (int i=0; i<neuts; ++i){
    if (gk->neut_species[i].collisions.is_implicit){
      has_implicit_coll_scheme = true;
    }
  }

  // Initialize EIRENE
  app->eirene = gk_eirene_init(app, gk);

  // Set the appropriate update function for taking a single time step
  // If we have implicit BGK collisions for either the gyrokinetic or neutral species, 
  // we perform a first-order operator split and treat those terms implicitly.
  // Otherwise, we default to an SSP-RK3 method. 
  if (has_implicit_coll_scheme) {
    app->update_func = gyrokinetic_update_op_split;
  }
  else {
    app->update_func = gyrokinetic_update_ssp_rk3;
  }

  // Pre-compute time-independent factors in omega_H.
  gkyl_gyrokinetic_app_omegaH_init(app); 

  // initialize stat object
  app->stat = (struct gkyl_gyrokinetic_stat) {
    .use_gpu = app->use_gpu,
    .stage_2_dt_diff = { DBL_MAX, 0.0 },
    .stage_3_dt_diff = { DBL_MAX, 0.0 },
  };

  app->dts = gkyl_dynvec_new(GKYL_DOUBLE, 1); // Dynvector to store time steps.
  app->is_first_dt_write_call = true;
}


gkyl_gyrokinetic_app*
gkyl_gyrokinetic_app_new(struct gkyl_gk *gk)
{
  gkyl_gyrokinetic_app* app = gkyl_gyrokinetic_app_new_geom(gk);
  gkyl_gyrokinetic_app_new_solver(gk, app);

  return app;
}

void
gyrokinetic_calc_field(gkyl_gyrokinetic_app* app, double tcurr,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  app->calc_field_func(app, tcurr, fin, bflux);
}

void
gyrokinetic_calc_field_and_apply_bc(gkyl_gyrokinetic_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array **bflux[], struct gkyl_array *distf_neut[])
{
  // Compute fields and apply BCs.

  // Compute the field.
  // MF 2024/09/27/: Need the cast here for consistency. Fixing
  // this may require removing 'const' from a lot of places.
  gyrokinetic_calc_field(app, tcurr, (const struct gkyl_array **) distf, bflux);

  // Apply boundary conditions.
  struct timespec wst = gkyl_wall_clock();
  for (int i=0; i<app->num_species; ++i) {
    gk_species_apply_bc(app, &app->species[i], distf[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    gk_neut_species_apply_bc(app, &app->neut_species[i], distf_neut[i]);
  }
  app->stat.bc_tm += gkyl_time_diff_now_sec(wst);
}

struct gk_species *
gk_find_species(const gkyl_gyrokinetic_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return &app->species[i];
  return 0;
}

int
gk_find_species_idx(const gkyl_gyrokinetic_app *app, const char *nm)
{
  for (int i=0; i<app->num_species; ++i)
    if (strcmp(nm, app->species[i].info.name) == 0)
      return i;
  return -1;
}

struct gk_neut_species *
gk_find_neut_species(const gkyl_gyrokinetic_app *app, const char *nm)
{
  for (int i=0; i<app->num_neut_species; ++i)
    if (strcmp(nm, app->neut_species[i].info.name) == 0)
      return &app->neut_species[i];
  return 0;
}

int
gk_find_neut_species_idx(const gkyl_gyrokinetic_app *app, const char *nm)
{
  for (int i=0; i<app->num_neut_species; ++i)
    if (strcmp(nm, app->neut_species[i].info.name) == 0)
      return i;
  return -1;
}

void
gkyl_gyrokinetic_app_apply_ic(gkyl_gyrokinetic_app* app, double t0)
{
  app->tcurr = t0;
  // Apply ICs that don't depend on other species.
  for (int i=0; i<app->num_species; ++i)
    gkyl_gyrokinetic_app_apply_ic_species(app, i, t0);

  for (int i=0; i<app->num_neut_species; ++i)
    gkyl_gyrokinetic_app_apply_ic_neut_species(app, i, t0);

  // Apply ICs that depend on other species.
  for (int i=0; i<app->num_species; ++i)
    gkyl_gyrokinetic_app_apply_ic_cross_species(app, i, t0);

  for (int i=0; i<app->num_neut_species; ++i)
    gkyl_gyrokinetic_app_apply_ic_cross_neut_species(app, i, t0);

  // Compute the fields and apply BCs.
  struct gkyl_array *distf[app->num_species];
  struct gkyl_array **bflux[app->num_species];
  struct gkyl_array *distf_neut[app->num_neut_species];
  for (int i=0; i<app->num_species; ++i) {
    distf[i] = app->species[i].f;
    bflux[i] = app->species[i].bflux.f;
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    distf_neut[i] = app->neut_species[i].f;
  }
  if (app->field->calc_init_field) {
    if (app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
      for (int i=0; i<app->num_species; ++i) {
        struct gk_species *s = &app->species[i];

        // Compute the collisionless flux so we can compute the initial boundary flux.
        gk_species_collisionless_flux(app, s, &s->collisionless, distf[i]);

        // Compute and store (in the ghost cell of out) the boundary fluxes.
        gk_species_bflux_rhs(app, &s->bflux, distf[i], distf[i]);
        // Compute moments of the boundary fluxes.
        gk_species_bflux_calc_moms(app, &s->bflux, distf[i], bflux[i]);

      }
    }

    if (app->field->info.init_from_file.type == 0 && app->field->info.init_field_profile == 0)
      // Compute the field.
      // MF 2024/09/27/: Need the cast here for consistency. Fixing
      // this may require removing 'const' from a lot of places.
      gyrokinetic_calc_field_enabled(app, t0, (const struct gkyl_array **) distf, bflux);
    else {
      if (app->field->info.init_field_profile == 0)
        // Read the field.
        gk_field_file_import_init(app, app->field->info.init_from_file);
      else
        // Project the field.
        gk_field_project_init(app);
    }

  }

  // Compute the phase-space advection speeds and boundary fluxes as t=0
  // diagnostics and emission BCs may need them.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_collisionless_flux(app, gks, &gks->collisionless, distf[i]);
    gk_species_bflux_rhs(app, &gks->bflux, gks->f, gks->f);
    gk_species_bflux_calc_moms(app, &gks->bflux, gks->f, bflux[i]);
  }

  // Apply boundary conditions.
  for (int i=0; i<app->num_species; ++i) {
    gk_species_apply_bc(app, &app->species[i], distf[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    if (!app->neut_species[i].info.is_static) {
      gk_neut_species_apply_bc(app, &app->neut_species[i], distf_neut[i]);
    }
  }
}

void
gkyl_gyrokinetic_app_apply_ic_species(gkyl_gyrokinetic_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  struct gk_species *gk_s = &app->species[sidx];

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  gk_species_apply_ic(app, gk_s, t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_gyrokinetic_app_apply_ic_neut_species(gkyl_gyrokinetic_app* app, int sidx, double t0)
{
  assert(sidx < app->num_neut_species);

  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  gk_neut_species_apply_ic(app, gk_ns, t0);
  app->stat.init_neut_species_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_gyrokinetic_app_apply_ic_cross_species(gkyl_gyrokinetic_app* app, int sidx, double t0)
{
  assert(sidx < app->num_species);

  struct gk_species *gk_s = &app->species[sidx];

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  gk_species_apply_ic_cross(app, gk_s, t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);
}

void
gkyl_gyrokinetic_app_apply_ic_cross_neut_species(gkyl_gyrokinetic_app* app, int sidx, double t0)
{
  assert(sidx < app->num_neut_species);

  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  app->tcurr = t0;
  struct timespec wtm = gkyl_wall_clock();
  gk_neut_species_apply_ic_cross(app, gk_ns, t0);
  app->stat.init_species_tm += gkyl_time_diff_now_sec(wtm);
}

//
// ............. Geometry outputs ............... //
// 
static void
gyrokinetic_app_geometry_copy_and_write(gkyl_gyrokinetic_app* app, struct gkyl_array *arr,
  struct gkyl_array *arr_host, char *varNm, struct gkyl_msgpack_data *mt)
{
  gkyl_array_copy(arr_host, arr);

  const char *fmt = "%s-%s.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, varNm);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, varNm);

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, arr_host, fileNm);
}

static void
gyrokinetic_app_geometry_copy_and_write_surf(gkyl_gyrokinetic_app* app, struct gkyl_array *arr,
  struct gkyl_array *arr_host, struct gkyl_array *arr_host_doubled, char *varNm, int dir, struct gkyl_msgpack_data *mt)
{
  gkyl_array_copy(arr_host, arr);

  gkyl_array_set_offset(arr_host_doubled, 1.0, arr_host, 0);
  gkyl_array_copy_range_to_range(arr_host, arr_host, &app->local_upper_skin[dir], &app->local_upper_ghost[dir]);
  gkyl_array_set_offset(arr_host_doubled, 1.0, arr_host, arr_host->ncomp);

  const char *fmt = "%s-%s_dir%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, app->name, varNm, dir);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, app->name, varNm, dir);

  gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, arr_host_doubled, fileNm);
}

void
gkyl_gyrokinetic_app_write_geometry(gkyl_gyrokinetic_app* app, struct gkyl_gk_geometry_inp *geometry_inp)
{
  // Package metadata.
  int io_meta_len[] = {app->io_meta_basic_len, app->io_meta_len, app->gk_geom->io_meta_len};
  const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic, app->io_meta, app->gk_geom->io_meta};
  struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);


  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (rank == 0 && geometry_inp->geometry_id == GKYL_GEOMETRY_TOKAMAK && gyrokinetic_str_ends_in_b0(app->name))
    gkyl_gk_geometry_write_efit(geometry_inp, app->io_meta_basic, app->io_meta_basic_len);

  // Gather geo into a global array
  struct gkyl_array* arr_ho1 = mkarr(false,   app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_hocdim = mkarr(false, app->cdim*app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho3 = mkarr(false, 3*app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho6 = mkarr(false, 6*app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho9 = mkarr(false, 9*app->basis.num_basis, app->local_ext.volume);

  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_corn.mc2p       , arr_ho3, "mapc2p", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_corn.mc2nu_pos  , arr_ho3, "mc2nu_pos", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_corn.bmag       , arr_ho1, "bmag_corn", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_corn.bmag_inv   , arr_ho1, "bmag_inv_corn", mt);
  if (app->cdim < 3) {
    if (geometry_inp->geometry_id == GKYL_GEOMETRY_MIRROR || geometry_inp->geometry_id == GKYL_GEOMETRY_TOKAMAK)
      gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_corn.mc2p_deflated, arr_hocdim, "mapc2p_deflated", mt);  
    gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_corn.mc2nu_pos_deflated, arr_hocdim, "mc2nu_pos_deflated", mt);  
  }

  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.mc2p       , arr_ho3, "mapc2pint", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.bmag        , arr_ho1, "bmag", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.g_ij        , arr_ho6, "g_ij", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.g_ij_neut   , arr_ho6, "g_ij_neut", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.dxdz        , arr_ho9, "dxdz", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.dzdx        , arr_ho9, "dzdx", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.normals     , arr_ho9, "normals", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.jacobgeo    , arr_ho1, "jacobgeo", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.jacobgeo_inv, arr_ho1, "jacobgeo_inv", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.gij         , arr_ho6, "gij", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.gij_neut    , arr_ho6, "gij_neut", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.b_i         , arr_ho3, "b_i", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.bcart       , arr_ho3, "bcart", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.cmag        , arr_ho1, "cmag", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.jacobtot    , arr_ho1, "jacobtot", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.jacobtot_inv, arr_ho1, "jacobtot_inv", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.gxxj        , arr_ho1, "gxxj", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.gxyj        , arr_ho1, "gxyj", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.gyyj        , arr_ho1, "gyyj", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.gxzj        , arr_ho1, "gxzj", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.eps2        , arr_ho1, "eps2", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.qprofile    , arr_ho1, "qprofile", mt);

  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.rtg33inv, arr_ho1, "rtg33inv", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.dualcurlbhatoverB, arr_ho3, "dualcurlbhatoverB", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.bioverJB, arr_ho3, "bioverJB", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.B3, arr_ho1, "B3", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.dualcurlbhat, arr_ho3, "dualcurlbhat", mt);
  gyrokinetic_app_geometry_copy_and_write(app, app->gk_geom->geo_int.bimpactangle, arr_ho1, "bimpactangle", mt);


  // Write surface quantities
  struct gkyl_array* arr_surf_ho1 = mkarr(false,   app->gk_geom->num_surf_basis, app->local_ext.volume);
  struct gkyl_array* arr_surf_ho2 = mkarr(false, 2*app->gk_geom->num_surf_basis, app->local_ext.volume);
  struct gkyl_array* arr_surf_ho3 = mkarr(false, 3*app->gk_geom->num_surf_basis, app->local_ext.volume);
  struct gkyl_array* arr_surf_ho6 = mkarr(false, 6*app->gk_geom->num_surf_basis, app->local_ext.volume);
  struct gkyl_array* arr_surf_ho9 = mkarr(false, 9*app->gk_geom->num_surf_basis, app->local_ext.volume);
  struct gkyl_array* arr_surf_ho18 = mkarr(false, 18*app->gk_geom->num_surf_basis, app->local_ext.volume);
  for (int dir = 0; dir<app->cdim; dir++ ) {
    gyrokinetic_app_geometry_copy_and_write_surf(app, app->gk_geom->geo_surf[dir].jacobgeo       , arr_surf_ho1, arr_surf_ho2, "jacobgeo", dir, mt);
    gyrokinetic_app_geometry_copy_and_write_surf(app, app->gk_geom->geo_surf[dir].jacobtot_inv   , arr_surf_ho1, arr_surf_ho2, "jacobtot_inv", dir, mt);
    gyrokinetic_app_geometry_copy_and_write_surf(app, app->gk_geom->geo_surf[dir].b_i            , arr_surf_ho3, arr_surf_ho6, "b_i", dir, mt);
    gyrokinetic_app_geometry_copy_and_write_surf(app, app->gk_geom->geo_surf[dir].cmag           , arr_surf_ho1, arr_surf_ho2, "cmag", dir, mt);
    gyrokinetic_app_geometry_copy_and_write_surf(app, app->gk_geom->geo_surf[dir].bmag           , arr_surf_ho1, arr_surf_ho2, "bmag", dir, mt);
    gyrokinetic_app_geometry_copy_and_write_surf(app, app->gk_geom->geo_surf[dir].B3           , arr_surf_ho1, arr_surf_ho2, "B3", dir, mt);
    gyrokinetic_app_geometry_copy_and_write_surf(app, app->gk_geom->geo_surf[dir].normcurlbhat           , arr_surf_ho1, arr_surf_ho2, "normcurlbhat", dir, mt);
    gyrokinetic_app_geometry_copy_and_write_surf(app, app->gk_geom->geo_surf[dir].normals           , arr_surf_ho9, arr_surf_ho18, "normals", dir, mt);
    gyrokinetic_app_geometry_copy_and_write_surf(app, app->gk_geom->geo_surf[dir].lenr           , arr_surf_ho1, arr_surf_ho2, "lenr", dir, mt);
  }

  // Write out nodes. This has to be done from rank 0 so we need to gather mc2p.
  struct gkyl_array *mc2p_global = mkarr(app->use_gpu, app->gk_geom->geo_corn.mc2p->ncomp, app->global_ext.volume);
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, app->gk_geom->geo_corn.mc2p, mc2p_global);
  struct gkyl_array *mc2p_global_ho = mkarr(false, mc2p_global->ncomp, mc2p_global->size);
  gkyl_array_copy(mc2p_global_ho, mc2p_global);

  struct gkyl_array *mc2pint_global = mkarr(app->use_gpu, app->gk_geom->geo_int.mc2p->ncomp, app->global_ext.volume);
  gkyl_comm_array_allgather(app->comm, &app->local, &app->global, app->gk_geom->geo_int.mc2p, mc2pint_global);
  struct gkyl_array *mc2pint_global_ho = mkarr(false, mc2pint_global->ncomp, mc2pint_global->size);
  gkyl_array_copy(mc2pint_global_ho, mc2pint_global);

  if (rank == 0) {
    // Create Nodal Range and Grid and Write Nodal Coordinates
    struct gkyl_range nrange;
    gkyl_gk_geometry_init_nodal_range(&nrange, &app->global, app->poly_order);
    struct gkyl_array* mc2p_nodal = mkarr(false, 3, nrange.volume);
    struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&app->basis, &app->grid, false);
    gkyl_nodal_ops_m2n(n2m, &app->basis, &app->grid, &nrange, &app->global, 3, mc2p_nodal, mc2p_global_ho, false);
    struct gkyl_rect_grid ngrid;
    gkyl_gk_geometry_init_nodal_grid(&ngrid, &app->grid, &nrange);

    const char *fmt = "%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, "nodes");
    char fileNm[sz+1]; // ensures no buffer overflow
    sprintf(fileNm, fmt, app->name, "nodes");

    // Package metadata for node file.
    int io_meta_nodes_len[] = {app->io_meta_basic_len, app->gk_geom->io_meta_len};
    const struct gkyl_msgpack_map_elem* io_meta_nodes[] = {app->io_meta_basic, app->gk_geom->io_meta};
    struct gkyl_msgpack_data *mt_nodes = gkyl_msgpack_create_union(sizeof(io_meta_nodes_len)/sizeof(int), io_meta_nodes_len, io_meta_nodes);

    gkyl_grid_sub_array_write(&ngrid, &nrange, mt_nodes,  mc2p_nodal, fileNm);

    gkyl_msgpack_data_release(mt_nodes);
    gkyl_nodal_ops_release(n2m);
    gkyl_array_release(mc2p_nodal);
  }

  if (rank == 0) {
    // Create Interior Nodal Range and Grid and Write Nodal Coordinates
    // nrange_int must be sized from app->global (not the per-rank local nrange_int)
    int num_quad_points_int = app->poly_order + 1;
    int num_nodes_int[GKYL_MAX_CDIM];
    for (int d=0; d<app->grid.ndim; ++d)
      num_nodes_int[d] = gkyl_range_shape(&app->global, d) * num_quad_points_int;
    struct gkyl_range nrange_int_global;
    gkyl_range_init_from_shape(&nrange_int_global, app->grid.ndim, num_nodes_int);
    struct gkyl_array* mc2pint_nodal = mkarr(false, 3, nrange_int_global.volume);
    struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&app->basis, &app->grid, false);
    gkyl_nodal_ops_m2n(n2m, &app->basis, &app->grid, &nrange_int_global, &app->global, 3, mc2pint_nodal, mc2pint_global_ho, true);

    struct gkyl_rect_grid ngrid_quad;
    gkyl_gk_geometry_init_nodal_grid(&ngrid_quad, &app->grid, &nrange_int_global);

    const char *fmt = "%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, "nodesint");
    char fileNm[sz+1]; // ensures no buffer overflow
    sprintf(fileNm, fmt, app->name, "nodesint");

    gkyl_grid_sub_array_write(&ngrid_quad, &nrange_int_global, 0, mc2pint_nodal, fileNm);

    gkyl_nodal_ops_release(n2m);
    gkyl_array_release(mc2pint_nodal);
  }

  gkyl_array_release(mc2p_global);
  gkyl_array_release(mc2p_global_ho);
  gkyl_array_release(mc2pint_global);
  gkyl_array_release(mc2pint_global_ho);
  gkyl_array_release(arr_ho1);
  gkyl_array_release(arr_hocdim);
  gkyl_array_release(arr_ho3);
  gkyl_array_release(arr_ho6);
  gkyl_array_release(arr_ho9);
  gkyl_array_release(arr_surf_ho1);
  gkyl_array_release(arr_surf_ho2);
  gkyl_array_release(arr_surf_ho3);
  gkyl_array_release(arr_surf_ho6);
  gkyl_array_release(arr_surf_ho9);
  gkyl_array_release(arr_surf_ho18);

  gkyl_msgpack_data_release(mt);
}

//
// ............. Field outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_field(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  if (app->field->update_field || frame == 0) {
    struct timespec wtm = gkyl_wall_clock();
    // Copy data from device to host before writing it out.
    if (app->use_gpu) {
      gkyl_array_copy(app->field->phi_host, app->field->phi_smooth);
    }

    // Package metadata.
    gkyl_msgpack_map_elem_set_double(app->io_meta_basic_len, app->io_meta_basic, "time", tm);
    gkyl_msgpack_map_elem_set_uint(app->io_meta_basic_len, app->io_meta_basic, "frame", frame);
    int io_meta_len[] = {app->io_meta_basic_len, app->io_meta_len, app->gk_geom->io_meta_len};
    const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic, app->io_meta, app->gk_geom->io_meta};
    struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

    const char *fmt = "%s-field_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, frame);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);

    gkyl_comm_array_write(app->comm, &app->grid, &app->local, mt, app->field->phi_host, fileNm);

    gkyl_msgpack_data_release(mt);

    app->stat.field_io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.n_field_io += 1;
  }
}

void
gkyl_gyrokinetic_app_calc_field_energy(gkyl_gyrokinetic_app* app, double tm)
{
  if (app->field->update_field) {
    struct timespec wst = gkyl_wall_clock();

    gk_field_calc_energy(app, tm, app->field);

    app->stat.field_diag_calc_tm += gkyl_time_diff_now_sec(wst);
    app->stat.n_field_diag += 1;
  }
}

void
gkyl_gyrokinetic_app_write_field_energy(gkyl_gyrokinetic_app* app)
{
  if (app->field->update_field) {
    struct timespec wtm = gkyl_wall_clock();
    // Write out the field energy.
    const char *fmt0 = "%s-field_energy.gkyl";
    int sz0 = gkyl_calc_strlen(fmt0, app->name);
    char fileNm0[sz0+1]; // ensures no buffer overflow
    snprintf(fileNm0, sizeof fileNm0, fmt0, app->name);

    int rank;
    gkyl_comm_get_rank(app->comm, &rank);

    if (rank == 0) {
      if (app->field->is_first_energy_write_call) {
        // Write to a new file (this ensure previous output is removed).
        gkyl_dynvec_write(app->field->integ_energy, fileNm0);
        app->field->is_first_energy_write_call = false;
      }
      else {
        // Append to existing file.
        gkyl_dynvec_awrite(app->field->integ_energy, fileNm0);
      }
      app->stat.n_field_diag_io += 1;
    }
    gkyl_dynvec_clear(app->field->integ_energy);

    if (app->field->info.time_rate_diagnostics) {
      // Write out the time rate of change of the field energy.
      const char *fmt1 = "%s-field_energy_dot.gkyl";
      int sz1 = gkyl_calc_strlen(fmt1, app->name);
      char fileNm1[sz1+1]; // ensures no buffer overflow
      snprintf(fileNm1, sizeof fileNm1, fmt1, app->name);

      if (rank == 0) {
        if (app->field->is_first_energy_dot_write_call) {
          // Write to a new file (this ensure previous output is removed).
          gkyl_dynvec_write(app->field->integ_energy_dot, fileNm1);
          app->field->is_first_energy_dot_write_call = false;
        }
        else {
          // Append to existing file.
          gkyl_dynvec_awrite(app->field->integ_energy_dot, fileNm1);
        }
        app->stat.n_field_diag_io += 1;
      }
      gkyl_dynvec_clear(app->field->integ_energy_dot);
    }

    app->stat.field_diag_io_tm += gkyl_time_diff_now_sec(wtm);
  }
}

//
// ............. Species outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_write(app, gks, tm, frame);
}

void
gkyl_gyrokinetic_app_write_neut_species(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_write(app, gkns, tm, frame);
}

void
gkyl_gyrokinetic_app_write_species_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_write_mom(app, gks, tm, frame);
}

void
gkyl_gyrokinetic_app_write_neut_species_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_write_mom(app, gkns, tm, frame);
}

void
gkyl_gyrokinetic_app_calc_species_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_calc_integrated_mom(app, gks, tm);
}

void
gkyl_gyrokinetic_app_calc_neut_species_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_calc_integrated_mom(app, gkns, tm);
}

void
gkyl_gyrokinetic_app_calc_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_bflux_calc_integrated_mom(app, gks, &gks->bflux, tm);
}

void
gkyl_gyrokinetic_app_calc_neut_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_bflux_calc_integrated_mom(app, gkns, &gkns->bflux, tm);
}

void
gkyl_gyrokinetic_app_write_species_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_write_integrated_mom(app, gks);
}

void
gkyl_gyrokinetic_app_write_neut_species_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_write_integrated_mom(app, gkns);
}

void
gkyl_gyrokinetic_app_calc_species_L2norm(gkyl_gyrokinetic_app *app, int sidx, double tm)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_calc_L2norm(app, gks, tm);
}

void
gkyl_gyrokinetic_app_write_species_L2norm(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_write_L2norm(app, gks);
}

void
gkyl_gyrokinetic_app_write_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_bflux_write_integrated_mom(app, gks, &gks->bflux);
}

void
gkyl_gyrokinetic_app_write_species_boundary_flux_mom(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_bflux_write_mom(app, gks, &gks->bflux, tm, frame);
}

void
gkyl_gyrokinetic_app_write_neut_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_bflux_write_integrated_mom(app, gkns, &gkns->bflux);
}

void
gkyl_gyrokinetic_app_write_neut_species_boundary_flux_mom(gkyl_gyrokinetic_app *app, int sidx, double tm, int frame)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_bflux_write_mom(app, gkns, &gkns->bflux, tm, frame);
}

//
// ............. Source outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_source(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_source_write(app, gks, tm, frame);
}

void
gkyl_gyrokinetic_app_write_neut_species_source(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_source_write(app, gkns, tm, frame);
}

void
gkyl_gyrokinetic_app_write_species_source_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_source_write_mom(app, gks, tm, frame);
}

void
gkyl_gyrokinetic_app_write_neut_species_source_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_source_write_mom(app, gkns, tm, frame);
}

void
gkyl_gyrokinetic_app_calc_species_source_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_source_calc_integrated_mom(app, gks, tm);
}

void
gkyl_gyrokinetic_app_calc_neut_species_source_integrated_mom(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_source_calc_integrated_mom(app, gkns, tm);
}

void
gkyl_gyrokinetic_app_write_species_source_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_source_write_integrated_mom(app, gks);
}

void
gkyl_gyrokinetic_app_write_neut_species_source_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_source_write_integrated_mom(app, gkns);
}

//
// ............. Damping outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_damping(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_damping_write(app, gks, tm, frame);
}

//
// ............. df/dt multiplier outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_fdot_multiplier(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_fdot_multiplier_write(app, gks, tm, frame);
}

//
// ............. BGK Source outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_source_bgk_diagnostics(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_source_bgk_write_diags(app, gks, &gks->bgk_src, tm, frame);
}

void
gkyl_gyrokinetic_app_calc_species_source_bgk_integrated_diagnostics(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_source_bgk_calc_integrated_diags(app, gks, &gks->bgk_src, tm);
}

void
gkyl_gyrokinetic_app_write_species_source_bgk_integrated_diagnostics(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_source_bgk_write_integrated_diags(app, gks, &gks->bgk_src);
}

//
// ............. EIRENE outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_eirene_diagnostics(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  gk_eirene_write(app, tm, frame);
}

void
gkyl_gyrokinetic_app_calc_eirene_integrated_diagnostics(gkyl_gyrokinetic_app* app, double tm)
{
  gk_eirene_calc_integrated_diagnostics(app, tm);
}

void
gkyl_gyrokinetic_app_write_eirene_integrated_diagnostics(gkyl_gyrokinetic_app* app)
{
  gk_eirene_write_integrated_diagnostics(app);
}

// ............. Collisionless outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_collisionless_diagnostics(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_collisionless_write_diags(app, gks, &gks->collisionless, tm, frame);
}

// ............. Positivity outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_positivity_diagnostics(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_positivity_write_diags(app, gks, &gks->positivity, tm, frame);
}

void
gkyl_gyrokinetic_app_calc_species_positivity_integrated_diagnostics(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_positivity_calc_integrated_diags(app, gks, &gks->positivity, tm);
}

void
gkyl_gyrokinetic_app_write_species_positivity_integrated_diagnostics(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_positivity_write_integrated_diags(app, gks, &gks->positivity);
}

void
gkyl_gyrokinetic_app_write_neut_species_positivity_diagnostics(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_positivity_write_diags(app, gkns, &gkns->positivity, tm, frame);
}

void
gkyl_gyrokinetic_app_calc_neut_species_positivity_integrated_diagnostics(gkyl_gyrokinetic_app* app, int sidx, double tm)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_positivity_calc_integrated_diags(app, gkns, &gkns->positivity, tm);
}

void
gkyl_gyrokinetic_app_write_neut_species_positivity_integrated_diagnostics(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_positivity_write_integrated_diags(app, gkns, &gkns->positivity);
}

//
// ............. LTE outputs ............... //
// 

void
gkyl_gyrokinetic_app_write_species_lte_max_corr_status(gkyl_gyrokinetic_app* app, int sidx)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_lte_write_max_corr_status(app, gks);  
}

void
gkyl_gyrokinetic_app_write_neut_species_lte_max_corr_status(gkyl_gyrokinetic_app* app, int sidx)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  gk_neut_species_lte_write_max_corr_status(app, gkns);  
}

//
// ............. Collision outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_lbo_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_lbo_write_mom(app, gks, tm, frame);
}

void
gkyl_gyrokinetic_app_write_species_bgk_mom(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_bgk_write_mom(app, gks, tm, frame);
}

//
// ............. Radiation outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_rad_drag(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_radiation_write_drag(app, gks, tm, frame);
}

void
gkyl_gyrokinetic_app_write_species_rad_emissivity(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_radiation_write_emissivity(app, gks, tm, frame);
}

void
gkyl_gyrokinetic_app_calc_species_rad_integrated_mom(gkyl_gyrokinetic_app *app, int sidx, double tm)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_radiation_calc_integrated_mom(app, gks, tm);
}

void
gkyl_gyrokinetic_app_write_species_rad_integrated_mom(gkyl_gyrokinetic_app *app, int sidx)
{
  struct gk_species *gks = &app->species[sidx];
  gk_species_radiation_write_integrated_mom(app, gks);
}

//
// ............. Neutral reaction outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_species_react(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  struct gk_react *gkr = &gks->react;
  gk_species_react_write(app, gks, gkr, ridx, tm, frame);
}

void
gkyl_gyrokinetic_app_write_species_react_neut(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame)
{
  struct gk_species *gks = &app->species[sidx];
  struct gk_react *gkr = &gks->react_neut;
  gk_species_react_write(app, gks, gkr, ridx, tm, frame);
}

void
gkyl_gyrokinetic_app_write_neut_species_react_neut(gkyl_gyrokinetic_app* app, int sidx, int ridx, double tm, int frame)
{
  struct gk_neut_species *gkns = &app->neut_species[sidx];
  struct gk_react *gkr = &gkns->react_neut;
  gk_neut_species_react_write(app, gkns, gkr, ridx, tm, frame);
}

//
// ............. Functions that group several outputs for a single species ............... //
// 
void
gkyl_gyrokinetic_app_write_species_phase(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_source(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_damping(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_fdot_multiplier(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_rad_drag(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_collisionless_diagnostics(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_app_write_neut_species_phase(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_neut_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_neut_species_source(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_app_write_species_conf(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_species_mom(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_source_mom(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_source_bgk_diagnostics(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_positivity_diagnostics(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_lbo_mom(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_bgk_mom(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_species_rad_emissivity(app, sidx, tm, frame);

  struct gk_species *gks = &app->species[sidx];
  for (int j=0; j<gks->react.num_react; ++j) {
    gkyl_gyrokinetic_app_write_species_react(app, sidx, j, tm, frame);
  }
  for (int j=0; j<gks->react_neut.num_react; ++j) {
    gkyl_gyrokinetic_app_write_species_react_neut(app, sidx, j, tm, frame);
  }

  gkyl_gyrokinetic_app_write_species_boundary_flux_mom(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_app_write_neut_species_conf(gkyl_gyrokinetic_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_neut_species_mom(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_neut_species_source_mom(app, sidx, tm, frame);

  gkyl_gyrokinetic_app_write_neut_species_positivity_diagnostics(app, sidx, tm, frame);

  struct gk_neut_species *gkns = &app->neut_species[sidx];

  for (int j=0; j<gkns->react_neut.num_react; ++j) {
    gkyl_gyrokinetic_app_write_neut_species_react_neut(app, sidx, j, tm, frame);
  }

  if (gkns->lower_bc[app->cdim-1].type == GKYL_BC_GK_SPECIES_RECYCLE)
    gk_neut_species_recycle_write_flux(app, gkns, &gkns->bc_recycle_lo, tm, frame);
  if (gkns->upper_bc[app->cdim-1].type == GKYL_BC_GK_SPECIES_RECYCLE)
    gk_neut_species_recycle_write_flux(app, gkns, &gkns->bc_recycle_up, tm, frame);

  gkyl_gyrokinetic_app_write_neut_species_boundary_flux_mom(app, sidx, tm, frame);
}

//
// ............. Functions that group several species outputs ............... //
// 
void
gkyl_gyrokinetic_app_write_mom(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_write_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_species_source_mom(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_species_source_bgk_diagnostics(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_species_positivity_diagnostics(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_species_lbo_mom(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_species_bgk_mom(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_species_rad_emissivity(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_species_boundary_flux_mom(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_app_write_neut_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_neut_species_source_mom(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_neut_species_positivity_diagnostics(app, i, tm, frame);
    gkyl_gyrokinetic_app_write_neut_species_boundary_flux_mom(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_app_calc_integrated_mom(gkyl_gyrokinetic_app* app, double tm)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_calc_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_app_calc_species_rad_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_app_calc_species_boundary_flux_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_app_calc_species_source_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_app_calc_species_source_bgk_integrated_diagnostics(app, i, tm);
    gkyl_gyrokinetic_app_calc_species_positivity_integrated_diagnostics(app, i, tm);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_app_calc_neut_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_app_calc_neut_species_source_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_app_calc_neut_species_positivity_integrated_diagnostics(app, i, tm);
    gkyl_gyrokinetic_app_calc_neut_species_boundary_flux_integrated_mom(app, i, tm);
  }

  gkyl_gyrokinetic_app_calc_eirene_integrated_diagnostics(app, tm);
}

void
gkyl_gyrokinetic_app_calc_L2norm(gkyl_gyrokinetic_app* app, double tm)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_calc_species_L2norm(app, i, tm);
  }
}

void
gkyl_gyrokinetic_app_write_integrated_mom(gkyl_gyrokinetic_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_write_species_integrated_mom(app, i);
    gkyl_gyrokinetic_app_write_species_source_integrated_mom(app, i);
    gkyl_gyrokinetic_app_write_species_source_bgk_integrated_diagnostics(app, i);
    gkyl_gyrokinetic_app_write_species_positivity_integrated_diagnostics(app, i);
    gkyl_gyrokinetic_app_write_species_lte_max_corr_status(app, i);
    gkyl_gyrokinetic_app_write_species_rad_integrated_mom(app, i);
    gkyl_gyrokinetic_app_write_species_boundary_flux_integrated_mom(app, i);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_app_write_neut_species_integrated_mom(app, i);
    gkyl_gyrokinetic_app_write_neut_species_source_integrated_mom(app, i);
    gkyl_gyrokinetic_app_write_neut_species_positivity_integrated_diagnostics(app, i);
    gkyl_gyrokinetic_app_write_neut_species_lte_max_corr_status(app, i);
    gkyl_gyrokinetic_app_write_neut_species_boundary_flux_integrated_mom(app, i);
  }

  gkyl_gyrokinetic_app_write_eirene_integrated_diagnostics(app);
}

void
gkyl_gyrokinetic_app_write_L2norm(gkyl_gyrokinetic_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_write_species_L2norm(app, i);
  }
}

void
gkyl_gyrokinetic_app_write_conf(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_field(app, tm, frame);

  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_write_species_conf(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_app_write_neut_species_conf(app, i, tm, frame);
  }

  gkyl_gyrokinetic_app_write_eirene_diagnostics(app, tm, frame);
}

void
gkyl_gyrokinetic_app_write_phase(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_app_write_species_phase(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_app_write_neut_species_phase(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_app_write(gkyl_gyrokinetic_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_app_write_phase(app, tm, frame);

  gkyl_gyrokinetic_app_write_conf(app, tm, frame);
}
//
// ............. End of write functions ............... //
// 

void
gyrokinetic_rhs(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  const struct gkyl_array *fin[], const struct gkyl_array *fbar_in[],
  struct gkyl_array *fout[], struct gkyl_array **bflux_out[], 
  const struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], struct gkyl_array **bflux_out_neut[], 
  struct gkyl_update_status *st)
{
  double dtmin = DBL_MAX;

  // Compute moments needed by various modules.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gk_s = &app->species[i];
    gk_species_lbo_moms(app, gk_s, &gk_s->lbo, fin[i]);
    gk_species_bgk_moms(app, gk_s, &gk_s->bgk, fin[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gk_ns = &app->neut_species[i];
    gk_neut_species_bgk_moms(app, gk_ns, &gk_ns->bgk, fin_neut[i]);
  }

  // Compute cross-species moments needed by various modules.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gk_s = &app->species[i];

    // Elastic collisions.
    gk_species_lbo_cross_moms(app, gk_s, &gk_s->lbo, fin[i]);        
    gk_species_bgk_cross_moms(app, gk_s, &gk_s->bgk, fin[i]);        

    // Reactions (e.g. ionization, recombination charge exchange).
    gk_species_react_cross_moms(app, gk_s, &gk_s->react, fin, fin_neut);
    gk_species_react_cross_moms(app, gk_s, &gk_s->react_neut, fin, fin_neut);

    // Line radiation.
    gk_species_radiation_moms(app, gk_s, &gk_s->rad, fin, fin_neut);

    // Scaling.
    gk_species_scaling_cross_moms(app, gk_s, &gk_s->sca, fin, fin_neut);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gk_ns = &app->neut_species[i];

    // Reactions (e.g. ionization, recombination charge exchange).
    gk_neut_species_react_cross_moms(app, gk_ns, &gk_ns->react_neut, fin, fin_neut);

    // Scaling.
    gk_neut_species_scaling_cross_moms(app, gk_ns, &gk_ns->sca, fin, fin_neut);
  }

  // Compute df/dt (not including sources).
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gk_s = &app->species[i];
    double dt1 = gk_species_rhs(app, gk_s, fin[i], fbar_in[i], fout[i], bflux_out[i]);
    dtmin = fmin(dtmin, dt1);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gk_ns = &app->neut_species[i];
    double dt1 = gk_neut_species_rhs(app, gk_ns, fin_neut[i], fout_neut[i], bflux_out_neut[i]);
    dtmin = fmin(dtmin, dt1);
  }

  gk_eirene_rhs(app, fin, fout);

  // Sources. Done after df/dt in case boundary fluxes are needed.
  for (int i=0; i<app->num_species; ++i) {
    gk_species_source_rhs(app, &app->species[i], 
      &app->species[i].src, fin[i], fout[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    gk_neut_species_source_rhs(app, &app->neut_species[i], 
      &app->neut_species[i].src, fin_neut[i], fout_neut[i]);
  }

  // Multiply dfdt by a factor.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gks = &app->species[i];
    gk_species_fdot_multiplier_advance_times_rate(app, gks, &gks->fdot_mult, app->field->phi_smooth, fout[i]);
  }

  struct timespec wtm = gkyl_wall_clock();
  double dt_max_rel_diff = 0.01;
  // Check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // Compute minimum time-step across all processors.
  double dtmin_local = dtmin, dtmin_global;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dtmin_local, &dtmin_global);
  dtmin = dtmin_global;
  
  // Don't take a time-step larger that input dt.
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;
  app->stat.dfdt_dt_reduce_tm += gkyl_time_diff_now_sec(wtm);
}

void
gyrokinetic_rhs_implicit(gkyl_gyrokinetic_app* app, double tcurr, double dt,
  struct gkyl_array *fin[], struct gkyl_array *fout[], struct gkyl_array **bflux_out[], 
  struct gkyl_array *fin_neut[], struct gkyl_array *fout_neut[], struct gkyl_array **bflux_out_neut[], 
  struct gkyl_update_status *st)
{
  // Compute moments needed by various modules.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gk_s = &app->species[i];
    gk_species_bgk_moms_implicit(app, gk_s, &gk_s->bgk, fin[i]);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gk_ns = &app->neut_species[i];
    gk_neut_species_bgk_moms_implicit(app, gk_ns, &gk_ns->bgk, fin_neut[i]);
  }

  // Compute cross-species moments needed by various modules.
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gk_s = &app->species[i];
    // Elastic collisions.
    gk_species_bgk_cross_moms_implicit(app, gk_s, &gk_s->bgk, fin[i]);        
  }

  // Compute df/dt from implicit terms (not including sources).
  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gk_s = &app->species[i];
    gk_species_rhs_implicit(app, gk_s, fin[i], fout[i], bflux_out[i], dt);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gk_ns = &app->neut_species[i];
    gk_neut_species_rhs_implicit(app, gk_ns, fin_neut[i], fout_neut[i], bflux_out_neut[i], dt);
  }
}

struct gkyl_update_status
gkyl_gyrokinetic_update(gkyl_gyrokinetic_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status = app->update_func(app, dt);
  app->tcurr += status.dt_actual;

  app->stat.time_loop_tm += gkyl_time_diff_now_sec(wst);
  // Check for any CUDA errors during time step
  if (app->use_gpu)
    checkCuda(cudaGetLastError());
  return status;
}

struct gkyl_gyrokinetic_stat
gkyl_gyrokinetic_app_stat(gkyl_gyrokinetic_app* app)
{
  struct gkyl_gyrokinetic_stat *stat = &app->stat;

  // Timers not yet computed in app directly.
  stat->time_rate_diags_tm = stat->fdot_tm + stat->phidot_tm;
  stat->pos_shift_tm = stat->species_pos_shift_tm + stat->neut_species_pos_shift_tm + stat->pos_shift_quasineut_tm;
  stat->io_tm = stat->species_io_tm + stat->species_diag_calc_tm + stat->species_diag_io_tm + stat->neut_species_io_tm
    + stat->neut_species_diag_calc_tm + stat->neut_species_diag_io_tm + stat->field_io_tm + stat->field_diag_calc_tm
    + stat->field_diag_io_tm + stat->app_io_tm;

  // Additions of several timers.
  stat->fwd_euler_sum_tm = stat->species_coll_mom_tm + stat->species_react_mom_tm + stat->neut_species_coll_mom_tm
    + stat->neut_species_react_mom_tm + stat->species_rad_mom_tm + stat->species_gyroavg_tm + stat->species_collisionless_tm
    + stat->species_coll_tm + stat->species_damp_tm + stat->species_fdot_mult_tm + stat->species_diffusion_tm
    + stat->species_rad_tm + stat->species_react_tm + stat->species_bflux_calc_tm + stat->species_bflux_moms_tm
    + stat->species_omega_cfl_tm + stat->species_src_tm + stat->species_source_bgk_tm
    + stat->neut_species_collisionless_tm + stat->neut_species_coll_tm + stat->neut_species_react_tm
    + stat->neut_species_omega_cfl_tm + stat->neut_species_src_tm + stat->dfdt_dt_reduce_tm + stat->fwd_euler_step_f_tm;
  stat->field_sum_tm = stat->field_phi_rhs_tm  + stat->field_phi_solve_tm;
  stat->bc_sum_tm = stat->species_bc_tm + stat->neut_species_bc_tm;
  stat->time_rate_diags_sum_tm = stat->fdot_tm + stat->phidot_tm;
  stat->pos_shift_sum_tm = stat->species_pos_shift_tm + stat->neut_species_pos_shift_tm + stat->pos_shift_quasineut_tm;
  stat->time_stepper_sum_tm = stat->fwd_euler_tm  + stat->field_tm + stat->bc_tm + stat->time_rate_diags_tm
    + stat->pos_shift_tm + stat->time_stepper_arithmetic_tm;
  stat->io_sum_tm = stat->species_io_tm + stat->species_diag_calc_tm + stat->species_diag_io_tm + stat->neut_species_io_tm
    + stat->neut_species_diag_calc_tm + stat->neut_species_diag_io_tm + stat->field_io_tm + stat->field_diag_calc_tm
    + stat->field_diag_io_tm + stat->app_io_tm;

  return *stat;
}

static inline
double
ratio_to_percent(double num, double den, double alt)
{
  return den > 1e-12 ? 100.*num/den : alt;
}

void
gkyl_gyrokinetic_app_print_timings(gkyl_gyrokinetic_app* app, FILE *iostream)
{
  struct gkyl_gyrokinetic_stat *stat = &app->stat;

  double bflux_tm = stat->species_bflux_calc_tm+stat->species_bflux_moms_tm;

  gkyl_gyrokinetic_app_cout(app, iostream, "Timing:\n");
  gkyl_gyrokinetic_app_cout(app, iostream, "  - Time loop:                         %.4e sec.\n", stat->time_loop_tm);
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Forward Euler:                   %.4e sec. / %4.2f %%.\n", stat->fwd_euler_tm, ratio_to_percent(stat->fwd_euler_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Collision moments (charged):   %.4e sec. / %4.2f %%.\n", stat->species_coll_mom_tm          , ratio_to_percent(stat->species_coll_mom_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Reaction moments (charged):    %.4e sec. / %4.2f %%.\n", stat->species_react_mom_tm         , ratio_to_percent(stat->species_react_mom_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Collision moments (neutral):   %.4e sec. / %4.2f %%.\n", stat->neut_species_coll_mom_tm     , ratio_to_percent(stat->neut_species_coll_mom_tm     ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Reaction moments (neutral):    %.4e sec. / %4.2f %%.\n", stat->neut_species_react_mom_tm    , ratio_to_percent(stat->neut_species_react_mom_tm    ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Radiation moments:             %.4e sec. / %4.2f %%.\n", stat->species_rad_mom_tm           , ratio_to_percent(stat->species_rad_mom_tm           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Species gyroaverage:           %.4e sec. / %4.2f %%.\n", stat->species_gyroavg_tm           , ratio_to_percent(stat->species_gyroavg_tm           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Species LTE:                   %.4e sec. / %4.2f %%.\n", stat->species_lte_tm               , ratio_to_percent(stat->species_lte_tm               ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Collisionless terms (charged): %.4e sec. / %4.2f %%.\n", stat->species_collisionless_tm     , ratio_to_percent(stat->species_collisionless_tm     ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Collision terms (charged):     %.4e sec. / %4.2f %%.\n", stat->species_coll_tm              , ratio_to_percent(stat->species_coll_tm              ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Damping (charged):             %.4e sec. / %4.2f %%.\n", stat->species_damp_tm              , ratio_to_percent(stat->species_damp_tm              ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ df/dt multiplier (charged):    %.4e sec. / %4.2f %%.\n", stat->species_fdot_mult_tm         , ratio_to_percent(stat->species_fdot_mult_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Diffusion (charged):           %.4e sec. / %4.2f %%.\n", stat->species_diffusion_tm         , ratio_to_percent(stat->species_diffusion_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Radiation terms:               %.4e sec. / %4.2f %%.\n", stat->species_rad_tm               , ratio_to_percent(stat->species_rad_tm               ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Reaction terms (charged):      %.4e sec. / %4.2f %%.\n", stat->species_react_tm             , ratio_to_percent(stat->species_react_tm             ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Boundary fluxes (charged):     %.4e sec. / %4.2f %%.\n", bflux_tm                           , ratio_to_percent(bflux_tm                           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ omega_cfl (charged):           %.4e sec. / %4.2f %%.\n", stat->species_omega_cfl_tm         , ratio_to_percent(stat->species_omega_cfl_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Sources (charged):             %.4e sec. / %4.2f %%.\n", stat->species_src_tm               , ratio_to_percent(stat->species_src_tm               ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ BGK Sources (charged):             %.4e sec. / %4.2f %%.\n", stat->species_source_bgk_tm               , ratio_to_percent(stat->species_source_bgk_tm               ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Collisionless terms (neutral): %.4e sec. / %4.2f %%.\n", stat->neut_species_collisionless_tm, ratio_to_percent(stat->neut_species_collisionless_tm,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Species LTE (neutral):         %.4e sec. / %4.2f %%.\n", stat->neut_species_lte_tm          , ratio_to_percent(stat->neut_species_lte_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Boundary fluxes (neutral):     %.4e sec. / %4.2f %%.\n", stat->neut_species_bflux_calc_tm +stat->neut_species_bflux_moms_tm, ratio_to_percent(stat->neut_species_bflux_calc_tm +stat->neut_species_bflux_moms_tm,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Collision terms (neutral):     %.4e sec. / %4.2f %%.\n", stat->neut_species_coll_tm         , ratio_to_percent(stat->neut_species_coll_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Reaction terms (neutral):      %.4e sec. / %4.2f %%.\n", stat->neut_species_react_tm        , ratio_to_percent(stat->neut_species_react_tm        ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ omega_cfl (neutral):           %.4e sec. / %4.2f %%.\n", stat->neut_species_omega_cfl_tm    , ratio_to_percent(stat->neut_species_omega_cfl_tm    ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Sources (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_src_tm          , ratio_to_percent(stat->neut_species_src_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Time step reduction:           %.4e sec. / %4.2f %%.\n", stat->dfdt_dt_reduce_tm            , ratio_to_percent(stat->dfdt_dt_reduce_tm            ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Step f:                        %.4e sec. / %4.2f %%.\n", stat->fwd_euler_step_f_tm          , ratio_to_percent(stat->fwd_euler_step_f_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->fwd_euler_sum_tm, stat->fwd_euler_tm, 100.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Field solves:                    %.4e sec. / %4.2f %%.\n", stat->field_tm, ratio_to_percent(stat->field_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Phi eqn RHS:                   %.4e sec. / %4.2f %%.\n", stat->field_phi_rhs_tm  , ratio_to_percent(stat->field_phi_rhs_tm  ,stat->field_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Phi eqn solve:                 %.4e sec. / %4.2f %%.\n", stat->field_phi_solve_tm, ratio_to_percent(stat->field_phi_solve_tm,stat->field_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->field_sum_tm, stat->field_tm, 100.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Boundary conditions::            %.4e sec. / %4.2f %%.\n", stat->bc_tm, ratio_to_percent(stat->bc_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Species (charged):             %.4e sec. / %4.2f %%.\n", stat->species_bc_tm     , ratio_to_percent(stat->species_bc_tm     ,stat->bc_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Species (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_bc_tm, ratio_to_percent(stat->neut_species_bc_tm,stat->bc_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->bc_sum_tm, stat->bc_tm, 100.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Time rate diagnostics:           %.4e sec. / %4.2f %%.\n", stat->time_rate_diags_tm, ratio_to_percent(stat->time_rate_diags_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Charged species:               %.4e sec. / %4.2f %%.\n", stat->fdot_tm  , ratio_to_percent(stat->fdot_tm  ,stat->time_rate_diags_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Phi:                           %.4e sec. / %4.2f %%.\n", stat->phidot_tm, ratio_to_percent(stat->phidot_tm,stat->time_rate_diags_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->time_rate_diags_sum_tm, stat->time_rate_diags_tm, 100.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Positivity:                      %.4e sec. / %4.2f %%.\n", stat->pos_shift_tm, ratio_to_percent(stat->pos_shift_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Species (charged):             %.4e sec. / %4.2f %%.\n", stat->species_pos_shift_tm     , ratio_to_percent(stat->species_pos_shift_tm     ,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Species (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_pos_shift_tm, ratio_to_percent(stat->neut_species_pos_shift_tm,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Quasineutrality:               %.4e sec. / %4.2f %%.\n", stat->pos_shift_quasineut_tm   , ratio_to_percent(stat->pos_shift_quasineut_tm   ,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->pos_shift_sum_tm, stat->pos_shift_tm, 100.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Time stepper arithmetic:         %.4e sec. / %4.2f %%.\n", stat->time_stepper_arithmetic_tm, ratio_to_percent(stat->time_stepper_arithmetic_tm,stat->time_loop_tm, 0.0));

  gkyl_gyrokinetic_app_cout(app, iostream, "    * Accounted for:                   %4.2f %%.\n", ratio_to_percent(stat->time_stepper_sum_tm, stat->time_loop_tm, 100.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "  - I/O:                               %.4e sec.\n", stat->io_tm);
  gkyl_gyrokinetic_app_cout(app, iostream, "    * f write (charged):               %.4e sec. / %4.2f %%.\n", stat->species_io_tm            , ratio_to_percent(stat->species_io_tm            , stat->io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Species diag calc (charged):     %.4e sec. / %4.2f %%.\n", stat->species_diag_calc_tm     , ratio_to_percent(stat->species_diag_calc_tm     , stat->io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Species diag write (charged):    %.4e sec. / %4.2f %%.\n", stat->species_diag_io_tm       , ratio_to_percent(stat->species_diag_io_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * f write (neutral):               %.4e sec. / %4.2f %%.\n", stat->neut_species_io_tm       , ratio_to_percent(stat->neut_species_io_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Species diag calc (neutral):     %.4e sec. / %4.2f %%.\n", stat->neut_species_diag_calc_tm, ratio_to_percent(stat->neut_species_diag_calc_tm, stat->io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Species diag write (neutral):    %.4e sec. / %4.2f %%.\n", stat->neut_species_diag_io_tm  , ratio_to_percent(stat->neut_species_diag_io_tm  , stat->io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Field write:                     %.4e sec. / %4.2f %%.\n", stat->field_io_tm              , ratio_to_percent(stat->field_io_tm              , stat->io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Field diag calc:                 %.4e sec. / %4.2f %%.\n", stat->field_diag_calc_tm       , ratio_to_percent(stat->field_diag_calc_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Field diag write:                %.4e sec. / %4.2f %%.\n", stat->field_diag_io_tm         , ratio_to_percent(stat->field_diag_io_tm         , stat->io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Common write:                    %.4e sec. / %4.2f %%.\n", stat->app_io_tm                , ratio_to_percent(stat->app_io_tm                , stat->io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, iostream, "    * Accounted for:                   %4.2f %%.\n", ratio_to_percent(stat->io_sum_tm, stat->io_tm, 100.0));
}

static void
range_stat_write(gkyl_gyrokinetic_app* app, const char *nm, const struct gkyl_range *r, FILE *fp)
{
  gkyl_gyrokinetic_app_cout(app, fp, " %s_cells : [ ", nm);
  for (int i=0; i<r->ndim; ++i)
    gkyl_gyrokinetic_app_cout(app, fp, " %d, ", gkyl_range_shape(r, i));
  gkyl_gyrokinetic_app_cout(app, fp, " ],\n");
}

// ensure stats across processors are made consistent
static void
comm_reduce_app_stat(const gkyl_gyrokinetic_app* app,
  const struct gkyl_gyrokinetic_stat *local, struct gkyl_gyrokinetic_stat *global)
{
  int comm_sz;
  gkyl_comm_get_size(app->comm, &comm_sz);
  if (comm_sz == 1) {
    memcpy(global, local, sizeof(struct gkyl_gyrokinetic_stat));
    return;
  }

  global->use_gpu = local->use_gpu;

  enum { NUP, NFEULER, NSTAGE_2_FAIL, NSTAGE_3_FAIL, L_END };
  int64_t l_red[] = {
    [NUP] = local->nup,
    [NFEULER] = local->nfeuler,
    [NSTAGE_2_FAIL] = local->nstage_2_fail,
    [NSTAGE_3_FAIL] = local->nstage_3_fail, 
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

  if (app->num_neut_species > 0) {
    int64_t l_red_neut_n_iter_corr[app->num_neut_species];
    int64_t l_red_neut_num_corr[app->num_neut_species];
    for (int s=0; s<app->num_neut_species; ++s) {
      l_red_neut_n_iter_corr[s] = local->neut_n_iter_corr[s];
      l_red_neut_num_corr[s] = local->neut_num_corr[s];
    }

    int64_t l_red_global_neut_n_iter_corr[app->num_neut_species];
    int64_t l_red_global_neut_num_corr[app->num_neut_species];
    gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, app->num_neut_species, 
      l_red_neut_n_iter_corr, l_red_global_neut_n_iter_corr);
    gkyl_comm_allreduce_host(app->comm, GKYL_INT_64, GKYL_MAX, app->num_neut_species, 
      l_red_neut_num_corr, l_red_global_neut_num_corr);

    for (int s=0; s<app->num_neut_species; ++s) {
      global->neut_n_iter_corr[s] = l_red_global_neut_n_iter_corr[s];
      global->neut_num_corr[s] = l_red_global_neut_num_corr[s];
    }
  }

  enum {
    INIT_SPECIES_TM, INIT_NEUT_SPECIES_TM,
    TIME_LOOP_TM, FWD_EULER_TM, FWD_EULER_STEP_F_TM, DFDT_DT_REDUCE_TM,
    SPECIES_COLLISIONLESS_TM, SPECIES_LTE_TM, SPECIES_GYROAVG_TM,
    SPECIES_BFLUX_CALC_TM, SPECIES_BFLUX_MOMS_TM, SPECIES_DAMP_TM, SPECIES_FDOT_MULT_TM, SPECIES_DIFFUSION_TM,
    SPECIES_COLL_MOM_TM, SPECIES_COLL_TM, 
    SPECIES_RAD_MOM_TM, SPECIES_RAD_TM, SPECIES_REACT_MOM_TM, SPECIES_REACT_TM, SPECIES_SRC_TM, SPECIES_SOURCE_BGK_TM, SPECIES_OMEGA_CFL_TM,
    NEUT_SPECIES_COLLISIONLESS_TM, NEUT_SPECIES_LTE_TM, NEUT_SPECIES_BFLUX_CALC_TM, NEUT_SPECIES_BFLUX_MOMS_TM,
    NEUT_SPECIES_COLL_MOM_TM, NEUT_SPECIES_COLL_TM, 
    NEUT_SPECIES_REACT_MOM_TM,  NEUT_SPECIES_REACT_TM, NEUT_SPECIES_SRC_TM, NEUT_SPECIES_OMEGA_CFL_TM,
    FDOT_TM, PHIDOT_TM, FIELD_TM, FIELD_PHI_RHS_TM, FIELD_PHI_SOLVE_TM,
    BC_TM, SPECIES_BC_TM, NEUT_SPECIES_BC_TM,
    TIME_STEPPER_ARITHMETIC_TM,
    SPECIES_POS_SHIFT_TM, NEUT_SPECIES_POS_SHIFT_TM, POS_SHIFT_QUASINEUT_TM,
    SPECIES_IO_TM, SPECIES_DIAG_CALC_TM, SPECIES_DIAG_IO_TM, 
    NEUT_SPECIES_IO_TM, NEUT_SPECIES_DIAG_CALC_TM, NEUT_SPECIES_DIAG_IO_TM, 
    FIELD_IO_TM, FIELD_DIAG_CALC_TM, FIELD_DIAG_IO_TM, 
    APP_IO_TM, IO_TM, 
    D_END
  };

  double d_red[D_END] = {
    [INIT_SPECIES_TM] = local->init_species_tm,
    [INIT_NEUT_SPECIES_TM] = local->init_neut_species_tm,
    [TIME_LOOP_TM] = local->time_loop_tm,
    [FWD_EULER_STEP_F_TM] = local->fwd_euler_step_f_tm,
    [DFDT_DT_REDUCE_TM] = local->dfdt_dt_reduce_tm,
    [SPECIES_COLLISIONLESS_TM] = local->species_collisionless_tm,
    [SPECIES_LTE_TM] = local->species_lte_tm,
    [SPECIES_GYROAVG_TM] = local->species_gyroavg_tm,
    [SPECIES_BFLUX_CALC_TM] = local->species_bflux_calc_tm,
    [SPECIES_BFLUX_MOMS_TM] = local->species_bflux_moms_tm,
    [SPECIES_DAMP_TM] = local->species_damp_tm,
    [SPECIES_FDOT_MULT_TM] = local->species_fdot_mult_tm,
    [SPECIES_DIFFUSION_TM] = local->species_diffusion_tm,
    [SPECIES_COLL_MOM_TM] = local->species_coll_mom_tm,
    [SPECIES_COLL_TM] = local->species_coll_tm,
    [SPECIES_RAD_MOM_TM] = local->species_rad_mom_tm,
    [SPECIES_RAD_TM] = local->species_rad_tm,
    [SPECIES_REACT_MOM_TM] = local->species_react_mom_tm,
    [SPECIES_REACT_TM] = local->species_react_tm,
    [SPECIES_SRC_TM] = local->species_src_tm,
    [SPECIES_SOURCE_BGK_TM] = local->species_source_bgk_tm,
    [SPECIES_OMEGA_CFL_TM] = local->species_omega_cfl_tm,
    [NEUT_SPECIES_COLLISIONLESS_TM] = local->neut_species_collisionless_tm,
    [NEUT_SPECIES_LTE_TM] = local->neut_species_lte_tm,
    [NEUT_SPECIES_BFLUX_CALC_TM] = local->neut_species_bflux_calc_tm,
    [NEUT_SPECIES_BFLUX_MOMS_TM] = local->neut_species_bflux_moms_tm,
    [NEUT_SPECIES_COLL_MOM_TM] = local->neut_species_coll_mom_tm,
    [NEUT_SPECIES_COLL_TM] = local->neut_species_coll_tm,
    [NEUT_SPECIES_REACT_MOM_TM] = local->neut_species_react_mom_tm,
    [NEUT_SPECIES_REACT_TM] = local->neut_species_react_tm,
    [NEUT_SPECIES_SRC_TM] = local->neut_species_src_tm,
    [NEUT_SPECIES_OMEGA_CFL_TM] = local->neut_species_omega_cfl_tm,
    [FDOT_TM] = local->fdot_tm,
    [PHIDOT_TM] = local->phidot_tm,
    [FIELD_TM] = local->field_tm,
    [FIELD_PHI_RHS_TM] = local->field_phi_rhs_tm,
    [FIELD_PHI_SOLVE_TM] = local->field_phi_solve_tm,
    [BC_TM] = local->bc_tm,
    [SPECIES_BC_TM] = local->species_bc_tm,
    [NEUT_SPECIES_BC_TM] = local->neut_species_bc_tm,
    [TIME_STEPPER_ARITHMETIC_TM] = local->time_stepper_arithmetic_tm,
    [SPECIES_POS_SHIFT_TM] = local->species_pos_shift_tm,
    [NEUT_SPECIES_POS_SHIFT_TM] = local->neut_species_pos_shift_tm,
    [POS_SHIFT_QUASINEUT_TM] = local->pos_shift_quasineut_tm,
    [SPECIES_IO_TM] = local->species_io_tm,
    [SPECIES_DIAG_CALC_TM] = local->species_diag_calc_tm,
    [SPECIES_DIAG_IO_TM] = local->species_diag_io_tm,
    [NEUT_SPECIES_IO_TM] = local->neut_species_io_tm,
    [NEUT_SPECIES_DIAG_CALC_TM] = local->neut_species_diag_calc_tm,
    [NEUT_SPECIES_DIAG_IO_TM] = local->neut_species_diag_io_tm,
    [FIELD_IO_TM] = local->field_io_tm,
    [FIELD_DIAG_CALC_TM] = local->field_diag_calc_tm,
    [FIELD_DIAG_IO_TM] = local->field_diag_io_tm,
    [APP_IO_TM] = local->app_io_tm,
    [IO_TM] = local->io_tm,
  };

  double d_red_global[D_END];
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, D_END, d_red, d_red_global);
  
  global->init_species_tm = d_red_global[INIT_SPECIES_TM];
  global->init_neut_species_tm = d_red_global[INIT_NEUT_SPECIES_TM];

  global->time_loop_tm = d_red_global[TIME_LOOP_TM];
  global->fwd_euler_step_f_tm = d_red_global[FWD_EULER_STEP_F_TM];
  global->dfdt_dt_reduce_tm = d_red_global[DFDT_DT_REDUCE_TM];

  global->species_collisionless_tm = d_red_global[SPECIES_COLLISIONLESS_TM];
  global->species_lte_tm = d_red_global[SPECIES_LTE_TM];
  global->species_gyroavg_tm = d_red_global[SPECIES_GYROAVG_TM];
  global->species_bflux_calc_tm = d_red_global[SPECIES_BFLUX_CALC_TM];
  global->species_bflux_moms_tm = d_red_global[SPECIES_BFLUX_MOMS_TM];
  global->species_damp_tm = d_red_global[SPECIES_DAMP_TM];
  global->species_fdot_mult_tm = d_red_global[SPECIES_FDOT_MULT_TM];
  global->species_diffusion_tm = d_red_global[SPECIES_DIFFUSION_TM];
  global->species_coll_mom_tm = d_red_global[SPECIES_COLL_MOM_TM];
  global->species_coll_tm = d_red_global[SPECIES_COLL_TM];
  global->species_rad_mom_tm = d_red_global[SPECIES_RAD_MOM_TM];
  global->species_rad_tm = d_red_global[SPECIES_RAD_TM];
  global->species_react_mom_tm = d_red_global[SPECIES_REACT_MOM_TM];
  global->species_react_tm = d_red_global[SPECIES_REACT_TM];
  global->species_src_tm = d_red_global[SPECIES_SRC_TM];
  global->species_source_bgk_tm = d_red_global[SPECIES_SOURCE_BGK_TM];
  global->species_omega_cfl_tm = d_red_global[SPECIES_OMEGA_CFL_TM];

  global->neut_species_collisionless_tm = d_red_global[NEUT_SPECIES_COLLISIONLESS_TM];
  global->neut_species_lte_tm = d_red_global[NEUT_SPECIES_LTE_TM];
  global->neut_species_bflux_calc_tm = d_red_global[NEUT_SPECIES_BFLUX_CALC_TM];
  global->neut_species_bflux_moms_tm = d_red_global[NEUT_SPECIES_BFLUX_MOMS_TM];
  global->neut_species_coll_mom_tm = d_red_global[NEUT_SPECIES_COLL_MOM_TM];
  global->neut_species_coll_tm = d_red_global[NEUT_SPECIES_COLL_TM];
  global->neut_species_react_mom_tm = d_red_global[NEUT_SPECIES_REACT_MOM_TM];
  global->neut_species_react_tm = d_red_global[NEUT_SPECIES_REACT_TM];
  global->neut_species_src_tm = d_red_global[NEUT_SPECIES_SRC_TM];
  global->neut_species_omega_cfl_tm = d_red_global[NEUT_SPECIES_OMEGA_CFL_TM];

  global->fdot_tm = d_red_global[FDOT_TM];
  global->phidot_tm = d_red_global[PHIDOT_TM];
  global->field_tm = d_red_global[FIELD_TM];
  global->field_phi_rhs_tm = d_red_global[FIELD_PHI_RHS_TM];
  global->field_phi_solve_tm = d_red_global[FIELD_PHI_SOLVE_TM];

  global->bc_tm = d_red_global[BC_TM];
  global->species_bc_tm = d_red_global[SPECIES_BC_TM];
  global->neut_species_bc_tm = d_red_global[NEUT_SPECIES_BC_TM];

  global->fdot_tm = d_red_global[FDOT_TM];
  global->phidot_tm = d_red_global[PHIDOT_TM];

  global->species_pos_shift_tm = d_red_global[SPECIES_POS_SHIFT_TM];
  global->neut_species_pos_shift_tm = d_red_global[NEUT_SPECIES_POS_SHIFT_TM];
  global->pos_shift_quasineut_tm = d_red_global[POS_SHIFT_QUASINEUT_TM];

  global->time_stepper_arithmetic_tm = d_red_global[TIME_STEPPER_ARITHMETIC_TM];

  global->species_io_tm = d_red_global[SPECIES_IO_TM];
  global->species_diag_calc_tm = d_red_global[SPECIES_DIAG_CALC_TM];
  global->species_diag_io_tm = d_red_global[SPECIES_DIAG_IO_TM];

  global->neut_species_io_tm = d_red_global[NEUT_SPECIES_IO_TM];
  global->neut_species_diag_calc_tm = d_red_global[NEUT_SPECIES_DIAG_CALC_TM];
  global->neut_species_diag_io_tm = d_red_global[NEUT_SPECIES_DIAG_IO_TM];

  global->field_io_tm = d_red_global[FIELD_IO_TM];
  global->field_diag_calc_tm = d_red_global[FIELD_DIAG_CALC_TM];
  global->field_diag_io_tm = d_red_global[FIELD_DIAG_IO_TM];

  global->app_io_tm = d_red_global[APP_IO_TM];
  global->io_tm = d_red_global[IO_TM];

  // misc data needing reduction

  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_2_dt_diff,
    global->stage_2_dt_diff);
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MAX, 2, local->stage_3_dt_diff,
    global->stage_3_dt_diff);
}

void
gkyl_gyrokinetic_app_stat_write(gkyl_gyrokinetic_app* app)
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

  for (int i=0; i<app->num_species; ++i) {
    struct gk_species *gk_s = &app->species[i];
    gk_species_n_iter_corr(app, gk_s, i); 
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    struct gk_neut_species *gk_ns = &app->neut_species[i];
    gk_neut_species_n_iter_corr(app, gk_ns, i);
  }

  struct gkyl_gyrokinetic_stat stat = { };
  comm_reduce_app_stat(app, &app->stat, &stat);
  
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  // append to existing file so we have a history of different runs
  FILE *fp = 0;
  if (rank == 0) fp = fopen(fileNm, "a");

  gkyl_gyrokinetic_app_cout(app, fp, "{\n");

  if (strftime(buff, sizeof buff, "%c", &curr_tm))
    gkyl_gyrokinetic_app_cout(app, fp, " date : %s,\n", buff);

  gkyl_gyrokinetic_app_cout(app, fp, " use_gpu : %d,\n", stat.use_gpu);
  gkyl_gyrokinetic_app_cout(app, fp, " num_ranks : %d,\n", num_ranks); 
  
  for (int s=0; s<app->num_species; ++s)
    range_stat_write(app, app->species[s].info.name, &app->species[s].global, fp);
  
  gkyl_gyrokinetic_app_cout(app, fp, " nup : %ld,\n", stat.nup);
  gkyl_gyrokinetic_app_cout(app, fp, " nfeuler : %ld,\n", stat.nfeuler);
  gkyl_gyrokinetic_app_cout(app, fp, " nstage_2_fail : %ld,\n", stat.nstage_2_fail);
  gkyl_gyrokinetic_app_cout(app, fp, " nstage_3_fail : %ld,\n", stat.nstage_3_fail);

  gkyl_gyrokinetic_app_cout(app, fp, " stage_2_dt_diff : [ %lg, %lg ],\n",
    stat.stage_2_dt_diff[0], stat.stage_2_dt_diff[1]);
  gkyl_gyrokinetic_app_cout(app, fp, " stage_3_dt_diff : [ %lg, %lg ],\n",
    stat.stage_3_dt_diff[0], stat.stage_3_dt_diff[1]);

  gkyl_gyrokinetic_app_cout(app, fp, " init_species_tm : %lg,\n", stat.init_species_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " init_neut_species_tm : %lg,\n", stat.init_neut_species_tm);

  for (int s=0; s<app->num_species; ++s) {
    gkyl_gyrokinetic_app_cout(app, fp, " n_iter_corr[%d] : %ld,\n", s, 
      stat.n_iter_corr[s]);
    gkyl_gyrokinetic_app_cout(app, fp, " num_corr[%d] : %ld,\n", s, 
      stat.num_corr[s]);          
  }
  for (int s=0; s<app->num_neut_species; ++s) {
    gkyl_gyrokinetic_app_cout(app, fp, " neut_n_iter_corr[%d] : %ld,\n", s, 
      stat.neut_n_iter_corr[s]);
    gkyl_gyrokinetic_app_cout(app, fp, " neut_num_corr[%d] : %ld,\n", s, 
      stat.neut_num_corr[s]);    
  }

  double bflux_tm = stat.species_bflux_calc_tm + stat.species_bflux_moms_tm;

  // Timers for the time-stepping loop
  gkyl_gyrokinetic_app_cout(app, fp, " time_loop_tm : %.4e,\n", stat.time_loop_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " fwd_euler_tm : %.4e,\n", stat.fwd_euler_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " fwd_euler_pct : %.2f,\n", ratio_to_percent(stat.fwd_euler_tm, stat.time_loop_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_coll_mom_tm : %.4e,\n", stat.species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_coll_mom_pct : %.2f,\n", ratio_to_percent(stat.species_coll_mom_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_react_mom_tm : %.4e,\n", stat.species_react_mom_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_react_mom_pct : %.2f,\n", ratio_to_percent(stat.species_react_mom_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_coll_mom_tm : %.4e,\n", stat.neut_species_coll_mom_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_coll_mom_pct : %.2f,\n", ratio_to_percent(stat.neut_species_coll_mom_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_react_mom_tm : %.4e,\n", stat.neut_species_react_mom_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_react_mom_pct : %.2f,\n", ratio_to_percent(stat.neut_species_react_mom_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_rad_mom_tm : %.4e,\n", stat.species_rad_mom_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_rad_mom_pct : %.2f,\n", ratio_to_percent(stat.species_rad_mom_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_gyroavg_tm : %.4e,\n", stat.species_gyroavg_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_gyroavg_pct : %.2f,\n", ratio_to_percent(stat.species_gyroavg_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_lte_tm : %.4e,\n", stat.species_lte_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_lte_pct : %.2f,\n", ratio_to_percent(stat.species_lte_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_collisionless_tm : %.4e,\n", stat.species_collisionless_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_collisionless_pct : %.2f,\n", ratio_to_percent(stat.species_collisionless_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_coll_tm : %.4e,\n", stat.species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_coll_pct : %.2f,\n", ratio_to_percent(stat.species_coll_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_diffusion_tm : %.4e,\n", stat.species_diffusion_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_diffusion_pct : %.2f,\n", ratio_to_percent(stat.species_diffusion_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_rad_tm : %.4e,\n", stat.species_rad_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_rad_pct : %.2f,\n", ratio_to_percent(stat.species_rad_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_react_tm : %.4e,\n", stat.species_react_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_react_pct : %.2f,\n", ratio_to_percent(stat.species_react_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_damp_tm : %.4e,\n", stat.species_damp_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_damp_pct : %.2f,\n", ratio_to_percent(stat.species_damp_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_fdot_mult_tm : %.4e,\n", stat.species_fdot_mult_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_fdot_mult_pct : %.2f,\n", ratio_to_percent(stat.species_fdot_mult_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_bflux_calc_tm : %.4e,\n", stat.species_bflux_calc_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_bflux_moms_tm : %.4e,\n", stat.species_bflux_moms_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_bflux_tm : %.4e,\n", bflux_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_bflux_pct : %.2f,\n", ratio_to_percent(bflux_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_omega_cfl_tm : %.4e,\n", stat.species_omega_cfl_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_omega_cfl_pct : %.2f,\n", ratio_to_percent(stat.species_omega_cfl_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_src_tm : %.4e,\n", stat.species_src_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_src_pct : %.2f,\n", ratio_to_percent(stat.species_src_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_source_bgk_tm : %.4e,\n", stat.species_source_bgk_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_source_bgk_pct : %.2f,\n", ratio_to_percent(stat.species_source_bgk_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_collisionless_tm : %.4e,\n", stat.neut_species_collisionless_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_collisionless_pct : %.2f,\n", ratio_to_percent(stat.neut_species_collisionless_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_lte_tm : %.4e,\n", stat.neut_species_lte_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_lte_pct : %.2f,\n", ratio_to_percent(stat.neut_species_lte_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_bflux_calc_tm : %.4e,\n", stat.neut_species_bflux_calc_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_bflux_moms_pct : %.2f,\n", ratio_to_percent(stat.neut_species_bflux_calc_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_bflux_moms_tm : %.4e,\n", stat.neut_species_bflux_moms_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_bflux_moms_pct : %.2f,\n", ratio_to_percent(stat.neut_species_bflux_moms_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_coll_tm : %.4e,\n", stat.neut_species_coll_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_coll_pct : %.2f,\n", ratio_to_percent(stat.neut_species_coll_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_react_tm : %.4e,\n", stat.neut_species_react_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_react_pct : %.2f,\n", ratio_to_percent(stat.neut_species_react_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_omega_cfl_tm : %.4e,\n", stat.neut_species_omega_cfl_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_omega_cfl_pct : %.2f,\n", ratio_to_percent(stat.neut_species_omega_cfl_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_src_tm : %.4e,\n", stat.neut_species_src_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_src_pct : %.2f,\n", ratio_to_percent(stat.neut_species_src_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " dfdt_dt_reduce_tm : %.4e,\n", stat.dfdt_dt_reduce_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " dfdt_dt_reduce_pct : %.2f,\n", ratio_to_percent(stat.dfdt_dt_reduce_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " fwd_euler_step_f_tm : %.4e,\n", stat.fwd_euler_step_f_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " fwd_euler_step_f_pct : %.2f,\n", ratio_to_percent(stat.fwd_euler_step_f_tm, stat.fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " fwd_euler_sum_tm : %.4e,\n", stat.fwd_euler_sum_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " fwd_euler_accounted_pct : %.2f,\n", ratio_to_percent(stat.fwd_euler_sum_tm, stat.fwd_euler_tm, 100.0));

  // Field solve timers
  gkyl_gyrokinetic_app_cout(app, fp, " field_tm : %.4e,\n", stat.field_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " field_pct : %.2f,\n", ratio_to_percent(stat.field_tm, stat.time_loop_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " field_phi_rhs_tm : %.4e,\n", stat.field_phi_rhs_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " field_phi_rhs_pct : %.2f,\n", ratio_to_percent(stat.field_phi_rhs_tm, stat.field_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " field_phi_solve_tm : %.4e,\n", stat.field_phi_solve_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " field_phi_solve_pct : %.2f,\n", ratio_to_percent(stat.field_phi_solve_tm, stat.field_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " field_sum_tm : %.4e,\n", stat.field_sum_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " field_accounted_pct : %.2f,\n", ratio_to_percent(stat.field_sum_tm, stat.field_tm, 100.0));

  // Boundary condition timers
  gkyl_gyrokinetic_app_cout(app, fp, " bc_tm : %.4e,\n", stat.bc_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " bc_pct : %.2f,\n", ratio_to_percent(stat.bc_tm, stat.time_loop_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_bc_tm : %.4e,\n", stat.species_bc_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_bc_pct : %.2f,\n", ratio_to_percent(stat.species_bc_tm, stat.bc_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_bc_tm : %.4e,\n", stat.neut_species_bc_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_bc_pct : %.2f,\n", ratio_to_percent(stat.neut_species_bc_tm, stat.bc_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " bc_sum_tm : %.4e,\n", stat.bc_sum_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " bc_accounted_pct : %.2f,\n", ratio_to_percent(stat.bc_sum_tm, stat.bc_tm, 100.0));

  // Time rate diagnostic timers
  gkyl_gyrokinetic_app_cout(app, fp, " time_rate_diags_tm : %.4e,\n", stat.time_rate_diags_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " time_rate_diags_pct : %.2f,\n", ratio_to_percent(stat.time_rate_diags_tm, stat.time_loop_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " fdot_tm : %.4e,\n", stat.fdot_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " fdot_pct : %.2f,\n", ratio_to_percent(stat.fdot_tm, stat.time_rate_diags_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " phidot_tm : %.4e,\n", stat.phidot_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " phidot_pct : %.2f,\n", ratio_to_percent(stat.phidot_tm, stat.time_rate_diags_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " time_rate_diags_sum_tm : %.4e,\n", stat.time_rate_diags_sum_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " time_rate_diags_accounted_pct : %.2f,\n", ratio_to_percent(stat.time_rate_diags_sum_tm, stat.time_rate_diags_tm, 100.0));

  // Position shift timers
  gkyl_gyrokinetic_app_cout(app, fp, " pos_shift_tm : %.4e,\n", stat.pos_shift_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " pos_shift_pct : %.2f,\n", ratio_to_percent(stat.pos_shift_tm, stat.time_loop_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_pos_shift_tm : %.4e,\n", stat.species_pos_shift_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_pos_shift_pct : %.2f,\n", ratio_to_percent(stat.species_pos_shift_tm, stat.pos_shift_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_pos_shift_tm : %.4e,\n", stat.neut_species_pos_shift_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_pos_shift_pct : %.2f,\n", ratio_to_percent(stat.neut_species_pos_shift_tm, stat.pos_shift_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " pos_shift_quasineut_tm : %.4e,\n", stat.pos_shift_quasineut_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " pos_shift_quasineut_pct : %.2f,\n", ratio_to_percent(stat.pos_shift_quasineut_tm, stat.pos_shift_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " pos_shift_sum_tm : %.4e,\n", stat.pos_shift_sum_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " pos_shift_accounted_pct : %.2f,\n", ratio_to_percent(stat.pos_shift_sum_tm, stat.pos_shift_tm, 100.0));

  // Time stepper arithmetic timers
  gkyl_gyrokinetic_app_cout(app, fp, " time_stepper_arithmetic_tm : %.4e,\n", stat.time_stepper_arithmetic_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " time_stepper_arithmetic_pct : %.2f,\n", ratio_to_percent(stat.time_stepper_arithmetic_tm, stat.time_loop_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " time_stepper_sum_tm : %.4e,\n", stat.time_stepper_sum_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " time_stepper_accounted_pct : %.2f,\n", ratio_to_percent(stat.time_stepper_sum_tm, stat.time_loop_tm, 100.0));

  // I/O timers
  gkyl_gyrokinetic_app_cout(app, fp, " io_tm : %.4e,\n", stat.io_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_io_tm : %.4e,\n", stat.species_io_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_io_pct : %.2f,\n", ratio_to_percent(stat.species_io_tm, stat.io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_diag_calc_tm : %.4e,\n", stat.species_diag_calc_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_diag_calc_pct : %.2f,\n", ratio_to_percent(stat.species_diag_calc_tm, stat.io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " species_diag_io_tm : %.4e,\n", stat.species_diag_io_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " species_diag_io_pct : %.2f,\n", ratio_to_percent(stat.species_diag_io_tm, stat.io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_io_tm : %.4e,\n", stat.neut_species_io_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_io_pct : %.2f,\n", ratio_to_percent(stat.neut_species_io_tm, stat.io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_diag_calc_tm : %.4e,\n", stat.neut_species_diag_calc_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_diag_calc_pct : %.2f,\n", ratio_to_percent(stat.neut_species_diag_calc_tm, stat.io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_diag_io_tm : %.4e,\n", stat.neut_species_diag_io_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " neut_species_diag_io_pct : %.2f,\n", ratio_to_percent(stat.neut_species_diag_io_tm, stat.io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " field_io_tm : %.4e,\n", stat.field_io_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " field_io_pct : %.2f,\n", ratio_to_percent(stat.field_io_tm, stat.io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " field_diag_calc_tm : %.4e,\n", stat.field_diag_calc_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " field_diag_calc_pct : %.2f,\n", ratio_to_percent(stat.field_diag_calc_tm, stat.io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " field_diag_io_tm : %.4e,\n", stat.field_diag_io_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " field_diag_io_pct : %.2f,\n", ratio_to_percent(stat.field_diag_io_tm, stat.io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " app_io_tm : %.4e,\n", stat.app_io_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " app_io_pct : %.2f,\n", ratio_to_percent(stat.app_io_tm, stat.io_tm, 0.0));
  gkyl_gyrokinetic_app_cout(app, fp, " io_sum_tm : %.4e,\n", stat.io_sum_tm);
  gkyl_gyrokinetic_app_cout(app, fp, " io_accounted_pct : %.2f,\n", ratio_to_percent(stat.io_sum_tm, stat.io_tm, 100.0));

  gkyl_gyrokinetic_app_cout(app, fp, " n_species_omega_cfl : %ld,\n", stat.n_species_omega_cfl);
  gkyl_gyrokinetic_app_cout(app, fp, " n_mom : %ld,\n", stat.n_mom);
  gkyl_gyrokinetic_app_cout(app, fp, " n_diag : %ld,\n", stat.n_diag);
  gkyl_gyrokinetic_app_cout(app, fp, " n_io : %ld,\n", stat.n_io);
  gkyl_gyrokinetic_app_cout(app, fp, " n_diag_io : %ld,\n", stat.n_diag_io);

  gkyl_gyrokinetic_app_cout(app, fp, " n_neut_species_omega_cfl : %ld,\n", stat.n_neut_species_omega_cfl);
  gkyl_gyrokinetic_app_cout(app, fp, " n_neut_mom : %ld,\n", stat.n_neut_mom);
  gkyl_gyrokinetic_app_cout(app, fp, " n_neut_diag : %ld,\n", stat.n_neut_diag);
  gkyl_gyrokinetic_app_cout(app, fp, " n_neut_io : %ld,\n", stat.n_neut_io);
  gkyl_gyrokinetic_app_cout(app, fp, " n_neut_diag_io : %ld,\n", stat.n_neut_diag_io);

  gkyl_gyrokinetic_app_cout(app, fp, " n_field_diag : %ld,\n", stat.n_field_diag);
  gkyl_gyrokinetic_app_cout(app, fp, " n_field_io : %ld,\n", stat.n_field_io);
  gkyl_gyrokinetic_app_cout(app, fp, " n_field_diag_io : %ld,\n", stat.n_field_diag_io);

  gkyl_gyrokinetic_app_cout(app, fp, "}\n");

  if (rank == 0)
    fclose(fp);  

}

void
gkyl_gyrokinetic_app_save_dt(gkyl_gyrokinetic_app* app, double tm, double dt)
{
  gkyl_dynvec_append(app->dts, tm, &dt);
}

void
gkyl_gyrokinetic_app_write_dt(gkyl_gyrokinetic_app* app)
{
  struct timespec wtm = gkyl_wall_clock();
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, "dt");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, "dt");

    if (app->is_first_dt_write_call) {
      gkyl_dynvec_write(app->dts, fileNm);
      app->is_first_dt_write_call = false;
    }
    else {
      gkyl_dynvec_awrite(app->dts, fileNm);
    }
  }
  gkyl_dynvec_clear(app->dts);

  app->stat.app_io_tm += gkyl_time_diff_now_sec(wtm);
  app->stat.n_diag_io += 1;
}

static struct gkyl_app_restart_status
header_from_file(gkyl_gyrokinetic_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = { .io_status = GKYL_ARRAY_RIO_FOPEN_FAILED };
  
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    struct gkyl_rect_grid grid;
    struct gkyl_array_header_info hdr;
    rstat.io_status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);

    if (GKYL_ARRAY_RIO_SUCCESS == rstat.io_status) {
      if (hdr.etype != GKYL_DOUBLE)
        rstat.io_status = GKYL_ARRAY_RIO_DATA_MISMATCH;
    }

    struct gkyl_msgpack_map_elem elem_list[] = {
      { .key = "frame", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
      { .key = "time", .elem_type = GKYL_MP_DOUBLE, .cval = 0 },
    };
    int elem_list_len = sizeof(elem_list)/sizeof(elem_list[0]);
    gkyl_msgpack_to_map_elem_list(&(struct gkyl_msgpack_data) {
        .meta = hdr.meta,
        .meta_sz = hdr.meta_size
      }, elem_list_len, elem_list);

    rstat.frame = gkyl_msgpack_map_elem_get_uint(elem_list_len, elem_list, "frame");
    rstat.stime = gkyl_msgpack_map_elem_get_double(elem_list_len, elem_list, "time");

    gkyl_grid_sub_array_header_release(&hdr);
  }
  
  return rstat;
}

//
// ............. Reading functions ............... //
// 
static void
gyrokinetic_app_geometry_read_and_copy(gkyl_gyrokinetic_app* app, struct gkyl_array *arr,
  struct gkyl_array *arr_host, char *path, char *varNm)
{
  cstr fileNm = cstr_from_fmt("%s%s-%s.gkyl", path, app->name, varNm);

  struct gkyl_app_restart_status rstat = header_from_file(app, fileNm.str);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, arr_host, fileNm.str);
    gkyl_array_copy(arr, arr_host);
  }
  else {
    gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read geometry file! (%s)\n",
        gkyl_array_rio_status_msg(rstat.io_status));
    assert(false);
  }
  cstr_drop(&fileNm);
}

static void
gyrokinetic_app_geometry_read_and_copy_surf(gkyl_gyrokinetic_app* app, struct gkyl_array *arr,
  struct gkyl_array *arr_host, struct gkyl_array *arr_host_doubled, char *path, char *varNm, int dir)
{
  cstr fileNm = cstr_from_fmt("%s%s-%s_dir%d.gkyl", path, app->name, varNm, dir);

  struct gkyl_app_restart_status rstat = header_from_file(app, fileNm.str);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, arr_host_doubled, fileNm.str);
    gkyl_array_set_offset_range(arr_host, 1.0, arr_host_doubled, 0, &app->local);
    gkyl_array_copy_range_to_range(arr_host_doubled, arr_host_doubled,
      &app->local_upper_ghost[dir], &app->local_upper_skin[dir]);
    gkyl_array_set_offset_range(arr_host, 1.0, arr_host_doubled, arr_host->ncomp, &app->local_upper_ghost[dir]);
    gkyl_array_copy(arr, arr_host);

  }
  else {
    gkyl_gyrokinetic_app_cout(app, stderr, "*** Failed to read geometry file! (%s)\n",
        gkyl_array_rio_status_msg(rstat.io_status));
    assert(false);
  }
  cstr_drop(&fileNm);
}

void
gkyl_gyrokinetic_app_read_geometry(gkyl_gyrokinetic_app* app, struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gkyl_array* arr_ho1 = mkarr(false,   app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho3 = mkarr(false, 3*app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho6 = mkarr(false, 6*app->basis.num_basis, app->local_ext.volume);
  struct gkyl_array* arr_ho9 = mkarr(false, 9*app->basis.num_basis, app->local_ext.volume);

  cstr fileNm = cstr_from_fmt("%s%s-%s.gkyl", geometry_inp->geometry_path, app->name, "jacobgeo");
  struct gkyl_array_header_info hdr;

  FILE *fp;
  with_file(fp, fileNm.str, "r") {

    int status = gkyl_grid_sub_array_header_read_fp(&app->grid, &hdr, fp);

    // Read geometry ID from header.
    struct gkyl_msgpack_map_elem elem_list[] = {
      { .key = "geometry_type", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
    };
    int elem_list_len = sizeof(elem_list)/sizeof(elem_list[0]);
    gkyl_msgpack_to_map_elem_list(&(struct gkyl_msgpack_data) {
        .meta = hdr.meta,
        .meta_sz = hdr.meta_size
      }, elem_list_len, elem_list);
    app->gk_geom->geometry_id = gkyl_msgpack_map_elem_get_uint(elem_list_len, elem_list, "geometry_type");

    if ((app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK) || (app->gk_geom->geometry_id == GKYL_GEOMETRY_MIRROR)) {
      // Read other metadata for numerical equilibrium from header.
      struct gkyl_msgpack_map_elem elem_list_numeq[] = {
        { .key = "geqdsk_sign_convention", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
      };
      int elem_list_numeq_len = sizeof(elem_list)/sizeof(elem_list[0]);
      gkyl_msgpack_to_map_elem_list(&(struct gkyl_msgpack_data) {
          .meta = hdr.meta,
          .meta_sz = hdr.meta_size
        }, elem_list_numeq_len, elem_list_numeq);
      app->gk_geom->geqdsk_sign_convention = gkyl_msgpack_map_elem_get_uint(elem_list_numeq_len,
        elem_list_numeq, "geqdsk_sign_convention");
    }

    if (app->gk_geom->geometry_id == GKYL_GEOMETRY_TOKAMAK) {
      // Read other metadata for numerical equilibrium from header.
      struct gkyl_msgpack_map_elem elem_list_numeq[] = {
        { .key = "half_domain", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 0 },
      };
      int elem_list_numeq_len = sizeof(elem_list)/sizeof(elem_list[0]);
      gkyl_msgpack_to_map_elem_list(&(struct gkyl_msgpack_data) {
          .meta = hdr.meta,
          .meta_sz = hdr.meta_size
        }, elem_list_numeq_len, elem_list_numeq);
      app->gk_geom->half_domain = gkyl_msgpack_map_elem_get_uint(elem_list_numeq_len,
        elem_list_numeq, "half_domain");
    }

    gkyl_grid_sub_array_header_release(&hdr);
  }
  cstr_drop(&fileNm);

  gkyl_gk_geometry_reset_io_meta(app->gk_geom); // Update metadata inside geo object.

  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_corn.mc2p        , arr_ho3, geometry_inp->geometry_path, "mapc2p");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_corn.mc2nu_pos   , arr_ho3, geometry_inp->geometry_path, "mc2nu_pos");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_corn.bmag        , arr_ho1, geometry_inp->geometry_path, "bmag_corn");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_corn.bmag_inv   , arr_ho1,  geometry_inp->geometry_path, "bmag_inv_corn");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.mc2p         , arr_ho3, geometry_inp->geometry_path, "mapc2pint");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.bmag        , arr_ho1,  geometry_inp->geometry_path, "bmag");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.g_ij        , arr_ho6,  geometry_inp->geometry_path, "g_ij");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.dxdz        , arr_ho9,  geometry_inp->geometry_path, "dxdz");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.dzdx        , arr_ho9,  geometry_inp->geometry_path, "dzdx");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.normals     , arr_ho9,  geometry_inp->geometry_path, "normals");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.jacobgeo    , arr_ho1,  geometry_inp->geometry_path, "jacobgeo");
  gkyl_array_copy(app->gk_geom->geo_int.jacobgeo_ghost, app->gk_geom->geo_int.jacobgeo);
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.jacobgeo_inv, arr_ho1,  geometry_inp->geometry_path, "jacobgeo_inv");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.gij         , arr_ho6,  geometry_inp->geometry_path, "gij");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.b_i         , arr_ho3,  geometry_inp->geometry_path, "b_i");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.bcart       , arr_ho3,  geometry_inp->geometry_path, "bcart");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.cmag        , arr_ho1,  geometry_inp->geometry_path, "cmag");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.jacobtot    , arr_ho1,  geometry_inp->geometry_path, "jacobtot");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.jacobtot_inv, arr_ho1,  geometry_inp->geometry_path, "jacobtot_inv");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.gxxj        , arr_ho1,  geometry_inp->geometry_path, "gxxj");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.gxyj        , arr_ho1,  geometry_inp->geometry_path, "gxyj");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.gyyj        , arr_ho1,  geometry_inp->geometry_path, "gyyj");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.gxzj        , arr_ho1,  geometry_inp->geometry_path, "gxzj");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.eps2        , arr_ho1,  geometry_inp->geometry_path, "eps2");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.qprofile    , arr_ho1,  geometry_inp->geometry_path, "qprofile");

  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.rtg33inv         , arr_ho1, geometry_inp->geometry_path, "rtg33inv");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.dualcurlbhatoverB, arr_ho3, geometry_inp->geometry_path, "dualcurlbhatoverB");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.bioverJB         , arr_ho3, geometry_inp->geometry_path, "bioverJB");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.B3               , arr_ho1, geometry_inp->geometry_path, "B3");
  gyrokinetic_app_geometry_read_and_copy(app, app->gk_geom->geo_int.dualcurlbhat     , arr_ho3, geometry_inp->geometry_path, "dualcurlbhat");

  struct gkyl_array* arr_surf_ho1 = mkarr(false,   app->gk_geom->num_surf_basis, app->local_ext.volume);
  struct gkyl_array* arr_surf_ho2 = mkarr(false, 2*app->gk_geom->num_surf_basis, app->local_ext.volume);
  struct gkyl_array* arr_surf_ho3 = mkarr(false, 3*app->gk_geom->num_surf_basis, app->local_ext.volume);
  struct gkyl_array* arr_surf_ho6 = mkarr(false, 6*app->gk_geom->num_surf_basis, app->local_ext.volume);
  struct gkyl_array* arr_surf_ho9 = mkarr(false, 9*app->gk_geom->num_surf_basis, app->local_ext.volume);
  struct gkyl_array* arr_surf_ho18 = mkarr(false, 18*app->gk_geom->num_surf_basis, app->local_ext.volume);
  for (int dir = 0; dir<app->cdim; dir++ ) {
    gyrokinetic_app_geometry_read_and_copy_surf(app, app->gk_geom->geo_surf[dir].jacobgeo     , arr_surf_ho1, arr_surf_ho2, geometry_inp->geometry_path, "jacobgeo", dir);
    gyrokinetic_app_geometry_read_and_copy_surf(app, app->gk_geom->geo_surf[dir].jacobtot_inv , arr_surf_ho1, arr_surf_ho2, geometry_inp->geometry_path, "jacobtot_inv", dir);
    gyrokinetic_app_geometry_read_and_copy_surf(app, app->gk_geom->geo_surf[dir].b_i          , arr_surf_ho3, arr_surf_ho6, geometry_inp->geometry_path, "b_i", dir);
    gyrokinetic_app_geometry_read_and_copy_surf(app, app->gk_geom->geo_surf[dir].cmag         , arr_surf_ho1, arr_surf_ho2, geometry_inp->geometry_path, "cmag", dir);
    gyrokinetic_app_geometry_read_and_copy_surf(app, app->gk_geom->geo_surf[dir].bmag         , arr_surf_ho1, arr_surf_ho2, geometry_inp->geometry_path, "bmag", dir);
    gyrokinetic_app_geometry_read_and_copy_surf(app, app->gk_geom->geo_surf[dir].normcurlbhat , arr_surf_ho1, arr_surf_ho2, geometry_inp->geometry_path, "normcurlbhat", dir);
    gyrokinetic_app_geometry_read_and_copy_surf(app, app->gk_geom->geo_surf[dir].B3           , arr_surf_ho1, arr_surf_ho2, geometry_inp->geometry_path, "B3", dir);
    gyrokinetic_app_geometry_read_and_copy_surf(app, app->gk_geom->geo_surf[dir].normals      , arr_surf_ho9, arr_surf_ho18,geometry_inp->geometry_path,  "normals", dir);
    gyrokinetic_app_geometry_read_and_copy_surf(app, app->gk_geom->geo_surf[dir].lenr         , arr_surf_ho1, arr_surf_ho2, geometry_inp->geometry_path, "lenr", dir);
    // jacobgeo_ratio is not used in single block.
    gkyl_array_clear(app->gk_geom->geo_surf[dir].jacobgeo_ratio, 0.0);
    gkyl_array_shiftc(app->gk_geom->geo_surf[dir].jacobgeo_ratio, pow(sqrt(2.0),app->cdim), 0);
  }

  gkyl_array_release(arr_ho1);
  gkyl_array_release(arr_ho3);
  gkyl_array_release(arr_ho6);
  gkyl_array_release(arr_ho9);
  gkyl_array_release(arr_surf_ho1);
  gkyl_array_release(arr_surf_ho2);
  gkyl_array_release(arr_surf_ho3);
  gkyl_array_release(arr_surf_ho6);
  gkyl_array_release(arr_surf_ho9);
  gkyl_array_release(arr_surf_ho18);
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_file_field(gkyl_gyrokinetic_app *app, const char *fname)
{
  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(app->comm, &app->grid, &app->local, app->field->phi_host, fname);
    if (app->use_gpu)
      gkyl_array_copy(app->field->phi_smooth, app->field->phi_host);
  }
  
  return rstat;
}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_app_from_file_species(gkyl_gyrokinetic_app *app, int sidx,
  const char *fname)
{
  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  struct gk_species *gk_s = &app->species[sidx];
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status = gkyl_comm_array_read(gk_s->comm, &gk_s->grid, &gk_s->local, gk_s->f_host, fname);
    if (app->use_gpu)
      gkyl_array_copy(gk_s->f, gk_s->f_host);

    if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
      gk_species_source_calc(app, gk_s, &gk_s->src, gk_s->lte.f_lte, 0.0);
      // Read volume and time integrated boundary flux diagnostics.
      gk_species_bflux_read_voltime_integrated_mom(app, gk_s, &gk_s->bflux);
    }
  }

  return rstat;
}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_app_from_file_neut_species(gkyl_gyrokinetic_app *app, int sidx,
  const char *fname)
{
  struct gkyl_app_restart_status rstat = header_from_file(app, fname);

  struct gk_neut_species *gk_ns = &app->neut_species[sidx];
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    rstat.io_status =
      gkyl_comm_array_read(gk_ns->comm, &gk_ns->grid, &gk_ns->local, gk_ns->f_host, fname);
    if (app->use_gpu)
      gkyl_array_copy(gk_ns->f, gk_ns->f_host);
    if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
      gk_neut_species_source_calc(app, gk_ns, &gk_ns->src, gk_ns->lte.f_lte, 0.0);
      // Read volume and time integrated boundary flux diagnostics.
      gk_neut_species_bflux_read_voltime_integrated_mom(app, gk_ns, &gk_ns->bflux);
    }
  }

  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_frame_field(gkyl_gyrokinetic_app *app, int frame)
{
  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, "field", frame);
  struct gkyl_app_restart_status rstat = gkyl_gyrokinetic_app_from_file_field(app, fileNm.str);
  app->field->is_first_energy_write_call = false; // Append to existing diagnostic.
  app->field->is_first_energy_dot_write_call = false; // Append to existing diagnostic.
  cstr_drop(&fileNm);
  
  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_frame_species(gkyl_gyrokinetic_app *app, int sidx, int frame)
{
  struct gk_species *gk_s = &app->species[sidx];

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, gk_s->info.name, frame);
  struct gkyl_app_restart_status rstat = gkyl_gyrokinetic_app_from_file_species(app, sidx, fileNm.str);
  cstr_drop(&fileNm);

  // Append to existing integrated diagnostics.
  app->is_first_dt_write_call = false;
  gk_s->is_first_integ_write_call = false;
  gk_s->is_first_L2norm_write_call = false;
  for (int b=0; b<gk_s->bflux.num_boundaries; ++b)
    gk_s->bflux.is_first_intmom_write_call[b] = false;
  if (gk_s->info.time_rate_diagnostics)
    gk_s->is_first_fdot_integ_write_call = false;
  if (gk_s->positivity.type)
    gk_s->positivity.is_first_integ_write_call = false;
  if (gk_s->rad.radiation_id == GKYL_GK_RADIATION)
    gk_s->rad.is_first_integ_write_call = false;
  if (gk_s->src.source_id)
    gk_s->src.is_first_integ_write_call = false;
  if (gk_s->lte.correct_all_moms)
    gk_s->lte.is_first_corr_status_write_call = false;

  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_app_from_frame_neut_species(gkyl_gyrokinetic_app *app, int sidx, int frame)
{
  struct gk_neut_species *gk_ns = &app->neut_species[sidx];

  cstr fileNm = cstr_from_fmt("%s-%s_%d.gkyl", app->name, gk_ns->info.name, frame);
  struct gkyl_app_restart_status rstat = gkyl_gyrokinetic_app_from_file_neut_species(app, sidx, fileNm.str);
  gk_ns->is_first_integ_write_call = false; // append to existing diagnostic
  cstr_drop(&fileNm);
  
  // Append to existing integrated diagnostics.
  gk_ns->is_first_integ_write_call = false;
  if (gk_ns->src.source_id) {
    gk_ns->src.is_first_integ_write_call = false;
  }
  if (gk_ns->positivity.type)
    gk_ns->positivity.is_first_integ_write_call = false;
  if (gk_ns->lte.correct_all_moms) {
    gk_ns->lte.is_first_corr_status_write_call = false;
  }

  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_app_read_from_frame(gkyl_gyrokinetic_app *app, int frame)
{
  struct gkyl_app_restart_status rstat;
  for (int i=0; i<app->num_neut_species; i++) {
    if (app->neut_species[i].info.is_static) {
      gk_neut_species_apply_ic(app, &app->neut_species[i], 0.0);
    }
    else {
      rstat = gkyl_gyrokinetic_app_from_frame_neut_species(app, i, frame);
    }
  }
  for (int i=0; i<app->num_species; i++) {
    if (app->species[i].info.is_static) {
      gk_species_apply_ic(app, &app->species[i], 0.0);
    }
    else {
      rstat = gkyl_gyrokinetic_app_from_frame_species(app, i, frame);
    }
  }

  // Apply ICs that depend on other species.
  for (int i=0; i<app->num_neut_species; ++i)
    gkyl_gyrokinetic_app_apply_ic_cross_neut_species(app, i, 0.0);
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    // Compute the fields and apply BCs.
    struct gkyl_array *distf[app->num_species];
    struct gkyl_array **bflux[app->num_species];
    struct gkyl_array *distf_neut[app->num_neut_species];
    for (int i=0; i<app->num_species; ++i) {
      distf[i] = app->species[i].f;
      bflux[i] = app->species[i].bflux.f;
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      distf_neut[i] = app->neut_species[i].f;
    }
    if (app->field->update_field) {
      if (app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
        for (int i=0; i<app->num_species; ++i) {
          struct gk_species *s = &app->species[i];

          // Compute the collisionless flux so we can compute the initial boundary flux.
          gk_species_collisionless_flux(app, s, &s->collisionless, distf[i]);

          // Compute and store (in the ghost cell of of out) the boundary fluxes.
          gk_species_bflux_rhs(app, &s->bflux, distf[i], distf[i]);
          // Compute moments of the boundary fluxes.
          gk_species_bflux_calc_moms(app, &s->bflux, distf[i], bflux[i]);
        }
      }

      // Compute the field.
      // MF 2024/09/27/: Need the cast here for consistency. Fixing
      // this may require removing 'const' from a lot of places.
      gyrokinetic_calc_field(app, rstat.stime, (const struct gkyl_array **) distf, bflux);
    }
    else {
      // Read the t=0 field.
      gkyl_gyrokinetic_app_from_frame_field(app, 0);
    }

    // Compute boundary fluxes, for recycling and diagnostics and adapt the source.
    for (int i=0; i<app->num_species; ++i) {
      struct gk_species *s = &app->species[i];

      // Compute the collisionless flux so we can compute the initial boundary flux.
      gk_species_collisionless_flux(app, s, &s->collisionless, distf[i]);

      // Compute and store (in the ghost cell of of out) the boundary fluxes.
      gk_species_bflux_rhs(app, &s->bflux, distf[i], distf[i]);
      // Compute moments of the boundary fluxes.
      gk_species_bflux_calc_moms(app, &s->bflux, distf[i], bflux[i]);

    }

    // Apply boundary conditions.
    for (int i=0; i<app->num_species; ++i) {
      gk_species_apply_bc(app, &app->species[i], distf[i]);
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      gk_neut_species_apply_bc(app, &app->neut_species[i], distf_neut[i]);
    }
  }
  app->field->is_first_energy_write_call = false; // Append to existing diagnostic.
  app->field->is_first_energy_dot_write_call = false; // Append to existing diagnostic.
  return rstat;
}

// private function to handle variable argument list for printing
static void
v_gk_app_cout(const gkyl_gyrokinetic_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank, r = 0;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp) {
    vfprintf(fp, fmt, argp);
    fflush(fp);
  }
}

void
gkyl_gyrokinetic_app_cout(const gkyl_gyrokinetic_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_gk_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_gyrokinetic_app_release_geom(gkyl_gyrokinetic_app* app)
{
  gkyl_gk_geometry_release(app->gk_geom);
  gkyl_position_map_release(app->position_map);
  for (int dir=0; dir<app->cdim; ++dir) {
    gkyl_rect_decomp_release(app->decomp_plane[dir]);
    gkyl_comm_release(app->comm_plane[dir]);
  }
  gkyl_comm_release(app->comm);
  gkyl_rect_decomp_release(app->decomp);

  if (app->use_gpu) {
    gkyl_cu_free(app->basis_on_dev);
  }
}

void
gkyl_gyrokinetic_app_release(gkyl_gyrokinetic_app* app)
{
  gyrokinetic_post_positivity_quasineut_release(app);

  gkyl_array_release(app->jacobtot_inv_weak);
  gkyl_gk_geometry_release(app->gk_geom);
  gkyl_dg_geom_release(app->dg_geom);
  gkyl_gk_dg_geom_release(app->gk_dg_geom);

  gk_field_release(app, app->field);
  gk_eirene_release(app, app->eirene);

  gkyl_position_map_release(app->position_map);

  for (int i=0; i<app->num_species; ++i)
    gk_species_release(app, &app->species[i]);
  if (app->num_species > 0)
    gkyl_free(app->species);

  for (int i=0; i<app->num_neut_species; ++i)
    gk_neut_species_release(app, &app->neut_species[i]);
  if (app->num_neut_species > 0)
    gkyl_free(app->neut_species);

  for (int dir=0; dir<app->cdim; ++dir) {
    gkyl_rect_decomp_release(app->decomp_plane[dir]);
    gkyl_comm_release(app->comm_plane[dir]);
  }
  gkyl_comm_release(app->comm);
  gkyl_rect_decomp_release(app->decomp);

  if (app->use_gpu) {
    gkyl_cu_free(app->basis_on_dev);
  }

  gkyl_dynvec_release(app->dts);

  gkyl_msgpack_map_elem_release(app->io_meta_basic_len, app->io_meta_basic);
  gkyl_msgpack_map_elem_release(app->io_meta_len, app->io_meta);

  gkyl_free(app);
}

void
gkyl_gyrokinetic_app_reset_cfl_frac_omegaH(gkyl_gyrokinetic_app* app, double tm,
  double cfl_frac_omegaH)
{
  double new_cfl_frac_omegaH = fabs(cfl_frac_omegaH) < 1e-16 ? 1.7 : cfl_frac_omegaH;
  app->cfl_omegaH = new_cfl_frac_omegaH;
}

void
gkyl_gyrokinetic_app_reset_species_fdot_multiplier(gkyl_gyrokinetic_app* app, double tm,
  const char* species_name, struct gkyl_gyrokinetic_fdot_multiplier fdot_mult_inp)
{
  struct gk_species *gks = gk_find_species(app, species_name);
  gk_species_fdot_multiplier_reset(app, tm, gks, &gks->fdot_mult, fdot_mult_inp);
}

void
gkyl_gyrokinetic_app_reset_species_collisionless(gkyl_gyrokinetic_app* app, double tm,
  const char *species_name, struct gkyl_gyrokinetic_collisionless gkcls_inp)
{
  struct gk_species *gks = gk_find_species(app, species_name);
  gk_species_collisionless_reset(app, tm, gks, &gks->collisionless, gkcls_inp);
}

void
gkyl_gyrokinetic_app_reset_species_positivity(gkyl_gyrokinetic_app* app, double tm,
  const char *species_name, struct gkyl_gyrokinetic_positivity pos_inp)
{
  struct gk_species *gks = gk_find_species(app, species_name);
  gk_species_positivity_reset(app, tm, gks, &gks->positivity, pos_inp);
}

void
gkyl_gyrokinetic_app_reset_species_damping(gkyl_gyrokinetic_app* app, double tm,
  const char *species_name, struct gkyl_gyrokinetic_damping damping_inp)
{
  struct gk_species *gks = gk_find_species(app, species_name);
  gk_species_damping_reset(app, tm, gks, &gks->damping, damping_inp);
}

void
gkyl_gyrokinetic_app_reset_field(gkyl_gyrokinetic_app* app, double tm,
  struct gkyl_gyrokinetic_field field_inp)
{
  app->field->info.is_static = field_inp.is_static;
  app->field->update_field = !field_inp.is_static;
  if (app->field->update_field)
    app->calc_field_func = gyrokinetic_calc_field_enabled;
  else
    app->calc_field_func = gyrokinetic_calc_field_disabled;
}
