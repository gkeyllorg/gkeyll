#include <gkyl_flash_block_geom.h>
#include <gkyl_alloc.h>
#include <gkyl_flash_multib.h>
#include <gkyl_flash_multib_priv.h>
#include <gkyl_flash.h>
#include <gkyl_flash_priv.h>
#include <gkyl_util.h>

#include <assert.h>
#include <stdarg.h>

// Compute total number of ranges specified by cuts.
static inline int
calc_cuts(int ndim, const int *cuts)
{
  int tc = 1;
  for (int d=0; d<ndim; ++d) tc *= cuts[d];
  return tc;
}

// Simple linear search to check if val occurs in lst.
static bool
has_int(int n, int val, const int *lst)
{
  for (int i=0; i<n; ++i)
    if (val == lst[i])
      return true;
  return false;
}

/**
 * Initialize a new gkyl_flash_app for a single block. This routine is 
 * assumed to be called by a parent app that has already initialized MPI
 * and performed the grid blockwise decomposition.
 * 
 * @param inp 
 * @return gkyl_flash_app* 
 */
static struct gkyl_flash_app *
gkyl_flash_singleb_app_new(const struct gkyl_flash_multib_inp *mbinp, int block_id,
  const struct gkyl_flash_multib_app *mbapp) 
{
  int cdim = mbapp->cdim, vdim = mbapp->vdim;

  const struct gkyl_flash_block_geom_info *bgi = gkyl_flash_block_geom_get_block(mbapp->block_geom, block_id);

  struct gkyl_flash_inp inp = {
    .cdim = cdim, .vdim = vdim
   };
   snprintf(inp.name, 128*sizeof(char), "%s_b%d", mbinp->name, block_id);

  for (int d=0; d<cdim; ++d) {
    inp.cells_fluid[d] = bgi->cells[d];
    inp.lower_conf[d] = bgi->lower[d];
    inp.upper_conf[d] = bgi->upper[d];
    inp.nghost_fluid[d] = mbinp->nghost_flu[d];
    inp.nghost_kinetic[d] = mbinp->nghost_kin[d];
  }

  inp.cfl_frac = mbinp->cfl_frac == 0.0 ? 1.0 : mbinp->cfl_frac;

  // Skip initializing a fluid species?
  inp.skip_fluid_species = mbinp->skip_fluid_species;

  // Set up species input parameters
  int num_species = inp.num_species = mbinp->num_species;
  for (int n=0; n<num_species; ++n) {
    const struct gkyl_flash_multib_species *mb_spec = &mbinp->species[n];

    struct gkyl_flash_species spec_inp = { };
    strcpy(spec_inp.name_kinetic, mb_spec->name_kin);
    strcpy(spec_inp.name_fluid, mb_spec->name_flu);

    spec_inp.mass = mb_spec->mass;
    spec_inp.gas_gamma = mb_spec->gas_gamma;
    spec_inp.limiter = mb_spec->limiter;
    // spec_inp.collisions = mb_spec->collisions;
    memcpy(&spec_inp.collisions, &mb_spec->collisions, sizeof(mb_spec->collisions));

    // Velocity space grid
    for (int v=0; v<vdim; ++v) {
      spec_inp.lower[v] = mb_spec->lower[v];
      spec_inp.upper[v] = mb_spec->upper[v];
      spec_inp.cells[v] = mb_spec->cells[v];
    }

    // Fetch the initial conditions input for the current block.
    const struct gkyl_flash_multib_species_inp_per_block *spec_pb = &mb_spec->blocks[0];
    if (!mb_spec->duplicate_ics) {
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        if (block_id == mb_spec->blocks[b].block_id) {
          spec_pb = &mb_spec->blocks[b];
          break;
        }
      }
    }

    // Copy initial condition parameters into species input
    spec_inp.init_dist = spec_pb->init_dist;
    spec_inp.init_fluid = spec_pb->init_fluid;
    spec_inp.ctx = spec_pb->ctx;
    spec_inp.skip_init_dist = spec_pb->skip_init_dist;
    spec_inp.init_dist_from_arrays = spec_pb->init_dist_from_arrays;
    spec_inp.init_dist_moms = spec_pb->init_dist_moms;

    // Set species block boundary conditions
    // By default, skip applying BCs in each direction
    for (int e=0; e<2; ++e) {
      spec_inp.bcx[e] = GKYL_SPECIES_SKIP;
      spec_inp.bcy[e] = GKYL_SPECIES_SKIP;
      spec_inp.bcz[e] = GKYL_SPECIES_SKIP;
    }

    // Only worry about physical (domain) BCs;
    // Inter-block communuication uses GKYL_SPECIES_SKIP
    for (int i=0; i<mb_spec->num_physical_bcs; ++i) {
      if (mb_spec->bcs[i].block_id != block_id) continue;

      int edge = mb_spec->bcs[i].edge;
      enum gkyl_species_bc_type bc_type = mb_spec->bcs[i].bc_type;
      if (mb_spec->bcs[i].dir == 0) 
        spec_inp.bcx[edge] = bc_type;
      else if (mb_spec->bcs[i].dir == 1) 
        spec_inp.bcy[edge] = bc_type;
      else if (mb_spec->bcs[i].dir == 2) 
        spec_inp.bcz[edge] = bc_type;
    }

    spec_inp.write_pkpm_diags = mb_spec->write_pkpm_diags;

    // Copy species input into app input
    memcpy(&inp.species[n], &spec_inp, sizeof(struct gkyl_flash_species));
  } 

  // Field input parameters
  const struct gkyl_flash_multib_field *mb_field = &mbinp->field;
  struct gkyl_flash_field field_inp = { };
  field_inp.is_static = mb_field->is_static;
  field_inp.epsilon0 = mb_field->epsilon0;
  field_inp.mu0 = mb_field->mu0;
  field_inp.limiter_fac = mb_field->limiter_fac;

  // Fetch field initial conditions for current block
  struct gkyl_flash_multib_field_inp_per_block *field_pb = &mb_field->blocks[0];
  if (!mb_field->duplicate_ics) {
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      if (block_id == mb_field->blocks[b].block_id) {
        field_pb = &mb_field->blocks[b];
        break;
      }
    }
  }
  field_inp.ctx = field_pb->ctx;
  field_inp.init = field_pb->init;
  field_inp.ext_em_ctx = field_pb->ext_em_ctx;
  field_inp.ext_em_init = field_pb->ext_em_init;
  for (int e=0; e<2; ++e) {
    field_inp.bcx[0] = field_pb->bcx[0];
    field_inp.bcy[0] = field_pb->bcy[0];
    field_inp.bcz[0] = field_pb->bcz[0];
  }

  // Copy field input into app input
  memcpy(&inp.field, &field_inp, sizeof(struct gkyl_flash_field));

  // Check for periodic domain boundaries
  inp.num_periodic_dir = mbinp->num_periodic_dir;
  for (int d=0; d<inp.num_periodic_dir; ++d) {
    inp.periodic_dirs[d] = mbinp->periodic_dirs[d];
  }

  inp.skip_fluid_species = mbinp->skip_fluid_species;

  return gkyl_flash_app_new(&inp);
}

static void
flash_multib_app_init_decomp(struct gkyl_flash_multib_app *mbapp, struct gkyl_flash_multib_inp *mbinp)
{
  int num_ranks;
  gkyl_comm_get_size(mbinp->comm, &num_ranks);
  int num_blocks = gkyl_flash_block_geom_num_blocks(mbapp->block_geom);

  int *branks = gkyl_malloc(num_blocks*sizeof(int));

  // Compute number of ranks in each block.
  for (int b=0; b<num_blocks; ++b) {
    const struct gkyl_flash_block_geom_info *bgi = gkyl_flash_block_geom_get_block(mbapp->block_geom, b);
    branks[b] = calc_cuts(mbapp->cdim, bgi->cuts);
  };
  
  mbapp->decomp_type = mbinp->decomp_type ? mbinp->decomp_type : GKYL_FLASH_MULTIB_DECOMP_RR;
  if (mbapp->decomp_type == GKYL_FLASH_MULTIB_DECOMP_RR) {
    mbapp->round_robin = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);
  }
  else if (mbapp->decomp_type == GKYL_FLASH_MULTIB_DECOMP_LINEAR) {
    mbapp->lin_decomp = gkyl_flash_decomp_new(num_ranks, num_blocks, mbapp->cdim, branks);
  }
  else {
    gkyl_exit("Unsupported domain decomposition type! Check decomp_type in gkyl_flash_multib_inp.");
  }

  gkyl_free(branks);
}

static void
flash_multib_app_getranks(struct gkyl_flash_multib_app *mbapp, int block_id, int *rank_list)
{

  if (mbapp->decomp_type == GKYL_FLASH_MULTIB_DECOMP_RR) {
    gkyl_rrobin_decomp_getranks(mbapp->round_robin, block_id, rank_list);
  }
  else {
    // Linear decomposition
    gkyl_flash_decomp_getranks(mbapp->lin_decomp, block_id, rank_list);
  }
}


gkyl_flash_multib_app*
gkyl_flash_multib_app_new(struct gkyl_flash_multib_inp *mbinp)
{
  struct gkyl_flash_multib_app *mbapp = gkyl_malloc(sizeof(struct gkyl_flash_multib_app));
  int my_rank, num_ranks;
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);
  gkyl_comm_get_rank(mbapp->comm, &my_rank);
  gkyl_comm_get_size(mbapp->comm, &num_ranks);

  mbapp->block_geom = gkyl_flash_block_geom_acquire(mbinp->block_geom);
  mbapp->topo = gkyl_flash_block_geom_topo(mbinp->block_geom);
  mbapp->cdim = mbinp->cdim; mbapp->vdim = mbinp->vdim;
  mbapp->cfl_frac = mbinp->cfl_frac == 0.0 ? 1.0 : mbinp->cfl_frac;
  strncpy(mbapp->name, mbinp->name, sizeof(mbinp->name));

  // Are we initializing a fluid species?
  mbapp->skip_fluid_species = mbinp->skip_fluid_species;

  int num_blocks = gkyl_flash_block_geom_num_blocks(mbapp->block_geom);
  int *branks = gkyl_malloc(num_blocks*sizeof(int));

  // Initialize global domain decomp, assigning ranks to blocks.
  for (int b=0; b<num_blocks; ++b) {
    const struct gkyl_flash_block_geom_info *bgi = gkyl_flash_block_geom_get_block(mbapp->block_geom, b);
    branks[b] = calc_cuts(mbapp->cdim, bgi->cuts);
  };
  flash_multib_app_init_decomp(mbapp, mbinp);
  mbapp->local_blocks = gkyl_malloc(num_blocks*sizeof(int));

  for (int b=0; b<num_blocks; ++b) {
    const struct gkyl_flash_block_geom_info *bgi = gkyl_flash_block_geom_get_block(mbapp->block_geom, b);
    branks[b] = calc_cuts(mbapp->cdim, bgi->cuts);
  };

  // Construct list of block communicators, one for each block.
  // Each communicator is not necessarily valid on each rank, but there are always
  // num_local_blocks valid communicators on each rank.
  mbapp->block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));

  mbapp->decomps_kin = gkyl_malloc(num_blocks*sizeof(struct gkyl_rect_decomp *));
  mbapp->decomps_flu = gkyl_malloc(num_blocks*sizeof(struct gkyl_rect_decomp *));
  int lidx = 0, num_local_blocks = 0;
  int *rank_list = gkyl_malloc(num_ranks*sizeof(int)); // This is larger than needed.
  for (int i=0; i<num_blocks; ++i) {
    flash_multib_app_getranks(mbapp, i, rank_list);

    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);

    if (is_my_rank_in_decomp) {
      mbapp->local_blocks[lidx++] = i;
      num_local_blocks += 1;
    }

    const struct gkyl_flash_block_geom_info *bgi = gkyl_flash_block_geom_get_block(mbapp->block_geom, i);
    struct gkyl_range block_range_kin, block_range_flu;

    // Kinetic range has half the number of cells as the fluid range.
    gkyl_create_global_range(mbapp->cdim, bgi->cells, &block_range_flu);

    int upper_kin[GKYL_MAX_DIM];
    for (int d=0; d<mbapp->cdim; ++d) {
      upper_kin[d] = (int)(block_range_flu.upper[d]/2);
    }
    gkyl_create_global_range(mbapp->cdim, upper_kin, &block_range_kin);

    mbapp->decomps_kin[i] = gkyl_rect_decomp_new_from_cuts(
      mbapp->cdim, bgi->cuts, &block_range_kin);

    if (!mbapp->skip_fluid_species) {
      mbapp->decomps_flu[i] = gkyl_rect_decomp_new_from_cuts(
        mbapp->cdim, bgi->cuts, &block_range_flu);
    }

    bool status;
    mbapp->block_comms[i] = gkyl_comm_create_comm_from_ranks(mbapp->comm,
      branks[i], rank_list, mbapp->decomps_kin[i], &status);
  }

  mbapp->num_local_blocks = num_local_blocks;
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    printf("Rank %d handles block %d\n", my_rank, mbapp->local_blocks[b]);
  }

  // Ghost cells in each configuration space dimension. Defaults must match
  // those of the single-block apps (1 kinetic, 2 fluid).
  int nghost_flu[3], nghost_kin[3];
  for (int d=0; d<mbapp->cdim; ++d) {
    nghost_flu[d] = mbinp->nghost_flu[d] != 0 ? mbinp->nghost_flu[d] : 2;
    nghost_kin[d] = mbinp->nghost_kin[d] != 0 ? mbinp->nghost_kin[d] : 1;
  }

  // Initialize single-block apps
  mbapp->block_apps = gkyl_malloc(mbapp->num_local_blocks*sizeof(struct gkyl_flash_app*));
  mbapp->num_species = mbinp->num_species;
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    int block_id = mbapp->local_blocks[b];
    mbapp->block_apps[b] = gkyl_flash_singleb_app_new(mbinp,  block_id, mbapp);
  }

  // Create send/recv connections for each block
  int block_id;
  mbapp->mbcc_sync_kin = gkyl_malloc(sizeof(struct gkyl_mbcc_sr));
  mbapp->mbcc_sync_kin->recv = gkyl_malloc(mbapp->num_local_blocks*sizeof(struct gkyl_multib_comm_conn*));
  mbapp->mbcc_sync_kin->send = gkyl_malloc(mbapp->num_local_blocks*sizeof(struct gkyl_multib_comm_conn*));
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    block_id = mbapp->local_blocks[b];

    flash_multib_app_getranks(mbapp, block_id, rank_list);
    int brank = -1;
    for (int i=0; i<branks[block_id]; ++i)
      if (rank_list[i] == my_rank) brank = i;
    
    struct gkyl_multib_comm_conn **recvs_kin = mbapp->mbcc_sync_kin->recv;
    struct gkyl_multib_comm_conn **sends_kin = mbapp->mbcc_sync_kin->send;
    
    sends_kin[b] = gkyl_multib_comm_conn_new_send(block_id, brank, nghost_kin, &mbapp->topo->conn[block_id], mbapp->decomps_kin);
    recvs_kin[b] = gkyl_multib_comm_conn_new_recv(block_id, brank, nghost_kin, &mbapp->topo->conn[block_id], mbapp->decomps_kin);

    // Set rank of partner block for send/recv operations
    struct gkyl_multib_comm_conn *mbcc_kin_s = sends_kin[b], *mbcc_kin_r = recvs_kin[b];
    for (int ns=0; ns<mbcc_kin_s->num_comm_conn; ++ns) {
      struct gkyl_comm_conn *ccs_kin = &mbcc_kin_s->comm_conn[ns];
      flash_multib_app_getranks(mbapp, ccs_kin->block_id, rank_list);
      ccs_kin->rank = rank_list[ccs_kin->rank];
      
      // Extend ranges to include velocity space
      struct gkyl_flash_app *sbapp = mbapp->block_apps[b];
      struct gkyl_range phase_range;
      gkyl_range_ten_prod(&phase_range, &ccs_kin->range, &sbapp->species[0].local_vel);

      // Convert data transfer range into a sub-range of the global block range (for indexing).
      gkyl_sub_range_init(&ccs_kin->range, &sbapp->species[0].local_ext, phase_range.lower, phase_range.upper);
    }

    for (int nr=0; nr<mbcc_kin_r->num_comm_conn; ++nr) {
      struct gkyl_comm_conn *ccr_kin = &mbcc_kin_r->comm_conn[nr];
      flash_multib_app_getranks(mbapp, ccr_kin->block_id, rank_list);
      ccr_kin->rank = rank_list[ccr_kin->rank];

      // Extend ranges to include velocity space
      struct gkyl_flash_app *sbapp = mbapp->block_apps[b];
      struct gkyl_range phase_range;
      gkyl_range_ten_prod(&phase_range, &ccr_kin->range, &sbapp->species[0].local_vel);

      // Convert data transfer range into a sub-range of the global block range (for indexing).
      gkyl_sub_range_init(&ccr_kin->range, &sbapp->species[0].local_ext, phase_range.lower, phase_range.upper);
    }
  }

  // Connections for configuration-space arrays on the kinetic grid
  mbapp->mbcc_sync_conf = gkyl_malloc(sizeof(struct gkyl_mbcc_sr));
  mbapp->mbcc_sync_conf->recv = gkyl_malloc(mbapp->num_local_blocks*sizeof(struct gkyl_multib_comm_conn*));
  mbapp->mbcc_sync_conf->send = gkyl_malloc(mbapp->num_local_blocks*sizeof(struct gkyl_multib_comm_conn*));
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    block_id = mbapp->local_blocks[b];
    struct gkyl_flash_app *sbapp = mbapp->block_apps[b];

    flash_multib_app_getranks(mbapp, block_id, rank_list);
    int brank = -1;
    for (int i=0; i<branks[block_id]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    struct gkyl_multib_comm_conn **sends_conf = mbapp->mbcc_sync_conf->send;
    struct gkyl_multib_comm_conn **recvs_conf = mbapp->mbcc_sync_conf->recv;

    sends_conf[b] = gkyl_multib_comm_conn_new_send(block_id, brank, nghost_kin, &mbapp->topo->conn[block_id], mbapp->decomps_kin);
    recvs_conf[b] = gkyl_multib_comm_conn_new_recv(block_id, brank, nghost_kin, &mbapp->topo->conn[block_id], mbapp->decomps_kin);

    // Set rank of partner block for send/recv operations
    struct gkyl_multib_comm_conn *mbcc_conf_s = sends_conf[b], *mbcc_conf_r = recvs_conf[b];
    for (int ns=0; ns<mbcc_conf_s->num_comm_conn; ++ns) {
      struct gkyl_comm_conn *ccs_conf = &mbcc_conf_s->comm_conn[ns];
      flash_multib_app_getranks(mbapp, ccs_conf->block_id, rank_list);
      ccs_conf->rank = rank_list[ccs_conf->rank];

      // Convert data transfer range into a sub-range of the global block range (for indexing).
      gkyl_sub_range_init(&ccs_conf->range, &sbapp->local_kinetic_ext, ccs_conf->range.lower, ccs_conf->range.upper);
    }
    for (int nr=0; nr<mbcc_conf_r->num_comm_conn; ++nr) {
      struct gkyl_comm_conn *ccr_conf = &mbcc_conf_r->comm_conn[nr];
      flash_multib_app_getranks(mbapp, ccr_conf->block_id, rank_list);
      ccr_conf->rank = rank_list[ccr_conf->rank];
      
      // Convert data transfer range into a sub-range of the global block range (for indexing).
      gkyl_sub_range_init(&ccr_conf->range, &sbapp->local_kinetic_ext, ccr_conf->range.lower, ccr_conf->range.upper);
    }
  }

  // If we have a fluid species, intialize send/recv block communications
  if (!mbapp->skip_fluid_species) {
    mbapp->mbcc_sync_flu = gkyl_malloc(sizeof(struct gkyl_mbcc_sr));
    mbapp->mbcc_sync_flu->recv = gkyl_malloc(mbapp->num_local_blocks*sizeof(struct gkyl_multib_comm_conn*));
    mbapp->mbcc_sync_flu->send = gkyl_malloc(mbapp->num_local_blocks*sizeof(struct gkyl_multib_comm_conn*));
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      block_id = mbapp->local_blocks[b];
      struct gkyl_flash_app *sbapp = mbapp->block_apps[b];

      flash_multib_app_getranks(mbapp, block_id, rank_list);
      int brank = -1;
      for (int i=0; i<branks[block_id]; ++i) {
        if (rank_list[i] == my_rank) brank = i;
      }

      struct gkyl_multib_comm_conn **recvs_flu = mbapp->mbcc_sync_flu->recv;
      struct gkyl_multib_comm_conn **sends_flu = mbapp->mbcc_sync_flu->send;

      sends_flu[b] = gkyl_multib_comm_conn_new_send(block_id, brank, nghost_flu, &mbapp->topo->conn[block_id], mbapp->decomps_flu);
      recvs_flu[b] = gkyl_multib_comm_conn_new_recv(block_id, brank, nghost_flu, &mbapp->topo->conn[block_id], mbapp->decomps_flu);

      // Set rank of partner block for send/recv operations
      struct gkyl_multib_comm_conn *mbcc_flu_s = sends_flu[b], *mbcc_flu_r = recvs_flu[b];
      for (int ns=0; ns<mbcc_flu_s->num_comm_conn; ++ns) {
        struct gkyl_comm_conn *ccs_flu = &mbcc_flu_s->comm_conn[ns];
        flash_multib_app_getranks(mbapp, ccs_flu->block_id, rank_list);
        ccs_flu->rank = rank_list[ccs_flu->rank];

        // Convert data transfer range into a sub-range of the global block range (for indexing).
        gkyl_sub_range_init(&ccs_flu->range, &sbapp->local_fluid_ext, ccs_flu->range.lower, ccs_flu->range.upper);
      }

      for (int nr=0; nr<mbcc_flu_r->num_comm_conn; ++nr) {
        struct gkyl_comm_conn *ccr_flu = &mbcc_flu_r->comm_conn[nr];
        flash_multib_app_getranks(mbapp, ccr_flu->block_id, rank_list);
        ccr_flu->rank = rank_list[ccr_flu->rank];

        // Convert data transfer range into a sub-range of the global block range (for indexing).
        gkyl_sub_range_init(&ccr_flu->range, &sbapp->local_fluid_ext, ccr_flu->range.lower, ccr_flu->range.upper);
      }
    }
  }

  mbapp->stat = (struct gkyl_flash_stat) {};

  gkyl_free(rank_list);
  gkyl_free(branks);

  return mbapp;
}

void
gkyl_flash_multib_app_apply_ic_species(gkyl_flash_multib_app *app, int sidx, double t0)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_flash_app_apply_ic_species(app->block_apps[b], sidx, t0);
  }
  gkyl_comm_barrier(app->comm);
}

void
gkyl_flash_multib_app_apply_ic_field(gkyl_flash_multib_app *app, double t0)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_flash_app_apply_ic_field(app->block_apps[b], t0);
  }
  gkyl_comm_barrier(app->comm);
}

/**
 * Apply boundary conditions to distribution function and fluid vars for each block.
 * 
 * @param app 
 * @param tcurr 
 * @param distf Pointer to array of distribution functions for each species on each block.
 * @param fluid 
 */
void
gkyl_flash_multib_app_apply_bc(gkyl_flash_multib_app *app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array *fluid[])
{
  // Apply physical boundary conditions to each block. Inter-block boundary conditions
  // are tagged with GKYL_SPECIES_SKIP and are handled after this.
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_flash_app *sbapp = app->block_apps[b];
    int lin_idx = b*app->num_species;
    for (int i=0; i<app->num_species; ++i) {
      flash_species_apply_bc(sbapp, &sbapp->species[i], distf[lin_idx+i]);

      if (!app->skip_fluid_species) 
        flash_fluid_species_apply_bc(sbapp, &sbapp->fluid_species[i], fluid[lin_idx+i]);
    }
  }

  // Communicate between blocks to fill ghost cells at inter-block boundaries
  for (int i=0; i<app->num_species; ++i) {
    // Create list of pointers to distf arrays on the local blocks of the current rank
    struct gkyl_array *distf_blocks[app->num_local_blocks];

    for (int b=0; b<app->num_local_blocks; ++b) {
      int lin_idx = b*app->num_species;
      distf_blocks[b] = distf[lin_idx+i];
    }

    // Sync distribution function
    gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
      app->mbcc_sync_kin->send, app->mbcc_sync_kin->recv, distf_blocks, distf_blocks);
  }

  // Sync inter-block boundaries for fluid variables if necessary
  if (!app->skip_fluid_species) {
    for (int i=0; i<app->num_species; ++i) {
      // Create list of pointers to distf arrays on the local blocks of the current rank
      struct gkyl_array *fluid_blocks[app->num_local_blocks];

      for (int b=0; b<app->num_local_blocks; ++b) {
        int lin_idx = b*app->num_species;
        fluid_blocks[b] = fluid[lin_idx+i];
      }

      // Sync fluid vars
      gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
        app->mbcc_sync_flu->send, app->mbcc_sync_flu->recv, fluid_blocks, fluid_blocks);
    }
  }
}

void
gkyl_flash_multib_app_apply_ic(gkyl_flash_multib_app *app, double t0)
{
  app->tcurr = t0;

  // Apply initial conditions to kinetic and fluid species (if applicable)
  for (int i=0; i<app->num_species; ++i) {
    gkyl_flash_multib_app_apply_ic_species(app, i, app->tcurr);
  }

  // Apply field initial conditions
  gkyl_flash_multib_app_apply_ic_field(app, t0);

  // Apply boundary conditions
  struct gkyl_array *distf[app->num_species * app->num_local_blocks];
  struct gkyl_array *fluid[app->num_species * app->num_local_blocks];
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_flash_app *sbapp = app->block_apps[b];
    int lin_idx = b*app->num_species;
    for (int i=0; i<app->num_species; ++i) {
      distf[lin_idx + i] = sbapp->species[i].f;
      fluid[lin_idx + i] = sbapp->fluid_species[i].fluid;
    }
  }

  gkyl_flash_multib_app_apply_bc(app, t0, distf, fluid);
}

double
gkyl_flash_multib_app_compute_dt_max(gkyl_flash_multib_app *app)
{
  // Compute maximum stable timestep as global minimum across all blocks/ranks.
  double dt_local = DBL_MAX;
  for (int b=0; b<app->num_local_blocks; ++b) {
    dt_local = fmin(dt_local, gkyl_flash_app_compute_dt_max(app->block_apps[b]));
  }

  double dt_global;
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 1, &dt_local, &dt_global);
  return dt_global;
}

struct gkyl_update_status
gkyl_flash_multib_kinetic_update(gkyl_flash_multib_app *app, double dt0)
{
  app->stat.nup += 1;

  struct timespec wst = gkyl_wall_clock();
  struct gkyl_update_status status = { .success = true };

  int ns = app->num_species;
  const struct gkyl_array *fin[ns];
  struct gkyl_array *fout[ns];

  // Compute maximum stable timestep as global minimum across all blocks/ranks.
  double dt = fmin(dt0, gkyl_flash_multib_app_compute_dt_max(app));

  // Call forward Euler update on local blocks, storing maximum stable timesteps.
  double dt_local[2] = { DBL_MAX, DBL_MAX }; // [dt_actual, dt_suggested]
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_flash_app *sbapp = app->block_apps[b];

    for (int i=0; i<ns; ++i) {
      fin[i] = sbapp->species[i].f;
      fout[i] = sbapp->species[i].f1;
    }

    struct gkyl_update_status st_b;
    flash_vlasov_forward_euler(sbapp, sbapp->tcurr, dt, fin, fout, &st_b);

    for (int i=0; i<ns; ++i) {
      gkyl_array_copy_range(sbapp->species[i].f, fout[i], &sbapp->species[i].local_ext);
    }

    dt_local[0] = fmin(dt_local[0], st_b.dt_actual);
    dt_local[1] = fmin(dt_local[1], st_b.dt_suggested);
  }

  // Global limiting timesteps across all ranks.
  double dt_global[2];
  gkyl_comm_allreduce_host(app->comm, GKYL_DOUBLE, GKYL_MIN, 2, dt_local, dt_global);
  status.dt_actual = dt_global[0];
  status.dt_suggested = dt_global[1];

  app->tcurr += status.dt_actual;

  app->stat.total_tm += gkyl_time_diff_now_sec(wst);

  // Apply boundary conditions
  struct gkyl_array *distf[app->num_species * app->num_local_blocks];
  struct gkyl_array *fluid[app->num_species * app->num_local_blocks];
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_flash_app *sbapp = app->block_apps[b];
    int lin_idx = b*app->num_species;
    for (int i=0; i<app->num_species; ++i) {
      distf[lin_idx + i] = sbapp->species[i].f;
      fluid[lin_idx + i] = sbapp->fluid_species[i].fluid;
    }
  }

  gkyl_flash_multib_app_apply_bc(app, app->tcurr, distf, fluid);
  return status;
}

struct gkyl_update_status
gkyl_flash_multib_update(gkyl_flash_multib_app *app, double dt0)
{
  app->stat.nup += 1;

  struct timespec wst = gkyl_wall_clock();
  struct gkyl_update_status status;
  if (!app->skip_fluid_species) {
    status = gkyl_flash_multib_update_ssp_rk3(app, dt0);
  } else {
    gkyl_exit("Calling gkyl_flash_multib_update without a fluid species! Try with skip_fluid_species = false in app initialization.");
  }
  app->tcurr += status.dt_actual;

  app->stat.total_tm += gkyl_time_diff_now_sec(wst);

  return status;
}

void
gkyl_flash_multib_app_write(gkyl_flash_multib_app *app, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_flash_app *sbapp = app->block_apps[b];
  
    for (int n=0; n<app->num_species; ++n) {
      gkyl_flash_app_write_species(sbapp, n, tm, frame);

      if (!app->skip_fluid_species)
        gkyl_flash_app_write_fluid_species(sbapp, n, tm, frame);

      if (sbapp->species[n].info.write_pkpm_diags) {
        gkyl_flash_app_write_pkpm_moms(sbapp, n, tm, frame);
      }
    }
  }
}

// Private function to handle variable argument list for printing.
static void
v_flash_multib_app_cout(const gkyl_flash_multib_app *app, FILE *fp, const char *fmt, va_list argp)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp) {
    vfprintf(fp, fmt, argp);
    fflush(fp);
  }
}

void
gkyl_flash_multib_app_cout(const gkyl_flash_multib_app *app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_flash_multib_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_flash_multib_app_release(gkyl_flash_multib_app *mbapp)
{
  int num_blocks = gkyl_flash_block_geom_num_blocks(mbapp->block_geom);

  // Inter-block send/recv connections.
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_kin->send[b]);
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_kin->recv[b]);
  }
  gkyl_free(mbapp->mbcc_sync_kin->recv);
  gkyl_free(mbapp->mbcc_sync_kin->send);
  gkyl_free(mbapp->mbcc_sync_kin);

  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->send[b]);
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->recv[b]);
  }
  gkyl_free(mbapp->mbcc_sync_conf->recv);
  gkyl_free(mbapp->mbcc_sync_conf->send);
  gkyl_free(mbapp->mbcc_sync_conf);

  if (!mbapp->skip_fluid_species) {
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_flu->send[b]);
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_flu->recv[b]);
    }
    gkyl_free(mbapp->mbcc_sync_flu->recv);
    gkyl_free(mbapp->mbcc_sync_flu->send);
    gkyl_free(mbapp->mbcc_sync_flu);
  }

  // Single-block apps.
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    gkyl_flash_app_release(mbapp->block_apps[b]);
  }
  gkyl_free(mbapp->block_apps);
  gkyl_free(mbapp->local_blocks);

  // Block communicators and decompositions (one per global block).
  for (int i=0; i<num_blocks; ++i) {
    gkyl_comm_release(mbapp->block_comms[i]);
    gkyl_rect_decomp_release(mbapp->decomps_kin[i]);
    if (!mbapp->skip_fluid_species)
      gkyl_rect_decomp_release(mbapp->decomps_flu[i]);
  }
  gkyl_free(mbapp->block_comms);
  gkyl_free(mbapp->decomps_kin);
  gkyl_free(mbapp->decomps_flu);

  if (mbapp->decomp_type == GKYL_FLASH_MULTIB_DECOMP_RR) {
    gkyl_rrobin_decomp_release(mbapp->round_robin);
  }
  else if (mbapp->decomp_type == GKYL_FLASH_MULTIB_DECOMP_LINEAR) {
    gkyl_flash_decomp_release(mbapp->lin_decomp);
  }

  gkyl_block_topo_release(mbapp->topo);
  gkyl_flash_block_geom_release(mbapp->block_geom);
  gkyl_comm_release(mbapp->comm);

  gkyl_free(mbapp);
}