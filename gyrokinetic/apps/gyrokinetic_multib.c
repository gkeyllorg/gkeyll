#include <gkyl_gk_block_geom.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>
#include <gkyl_multib_conn.h>

#include <mpack.h>

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

// Compute total and maximum number of cuts.
static void
calc_tot_and_max_cuts(const struct gkyl_gk_block_geom *gk_block_geom, int tot_max[2])
{
  int ndim = gkyl_gk_block_geom_ndim(gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(gk_block_geom);

  int max_cuts = 0, tot_cuts = 0;
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(gk_block_geom, i);
    int ncuts = calc_cuts(ndim, bgi->cuts);
    max_cuts = ncuts > max_cuts ? ncuts : max_cuts;
    tot_cuts += ncuts;
  }
  tot_max[0] = tot_cuts;
  tot_max[1] = max_cuts;
}

// Construct single-block App geometry for given block ID.
static struct gkyl_gyrokinetic_app *
singleb_app_new_geom(const struct gkyl_gyrokinetic_multib *mbinp, int bid,
  const struct gkyl_gyrokinetic_multib_app *mbapp)
{
  // For kinetic simulations, block dimension defined configuration-space dimensionality.
  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  const struct gkyl_gk_block_geom_info *bgi =
    gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, bid);

  // Construct top-level single-block input struct.
  struct gkyl_gk app_inp = { };

  strcpy(app_inp.name, mbinp->name);
  if (num_blocks > 1) {
    cstr app_name = cstr_from_fmt("%s_b%d", mbinp->name, bid);
    strcpy(app_inp.name, app_name.str);
    cstr_drop(&app_name);
  }

  // Set the configuration-space extents, cells, and geometry.
  app_inp.cdim = cdim;
  for (int i=0; i<cdim; ++i) {
    app_inp.lower[i] = bgi->lower[i];
    app_inp.upper[i] = bgi->upper[i];
    app_inp.cells[i] = bgi->cells[i];
  }

  // Set z dir grid extents based on tokamak global normalization
  if (bgi->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK || bgi->geometry.geometry_id == GKYL_GEOMETRY_FROMFILE) {
    gkyl_gk_geometry_tok_set_grid_extents(bgi->geometry.efit_info, bgi->geometry.tok_grid_info, &app_inp.lower[cdim-1], &app_inp.upper[cdim-1]);
    gkyl_gk_block_geom_reset_block_extents(mbapp->gk_block_geom, bid, app_inp.lower, app_inp.upper);
  }

  app_inp.geometry = bgi->geometry;
  int num_species = app_inp.num_species = mbinp->num_species;
  int num_neut_species = app_inp.num_neut_species = mbinp->num_neut_species; 

  app_inp.poly_order = mbinp->poly_order;
  app_inp.basis_type = mbinp->basis_type;
  app_inp.cfl_frac = mbinp->cfl_frac; 
  app_inp.cfl_frac_omegaH = mbinp->cfl_frac_omegaH; 

  struct gkyl_comm *comm = mbapp->block_comms[bid];

  struct gkyl_app_parallelism_inp parallel_inp = {};
  parallel_inp.use_gpu = mbinp->use_gpu;
  for (int d=0; d<cdim; ++d) parallel_inp.cuts[d] = bgi->cuts[d];
  parallel_inp.comm = comm;
  // Copy parallelism input into app input.
  memcpy(&app_inp.parallelism, &parallel_inp, sizeof(struct gkyl_app_parallelism_inp));

  app_inp.num_periodic_dir = mbinp->num_periodic_dir;
  for(int i = 0; i < mbinp->cdim; i++)
    app_inp.periodic_dirs[i] = mbinp->periodic_dirs[i];

  app_inp.metadata.num_attributes = mbapp->io_meta_basic_len;
  app_inp.metadata.attributes = mbapp->io_meta_basic;

  return gkyl_gyrokinetic_app_new_geom(&app_inp);
}

// Construct single-block App solver for given block ID.
static void
singleb_app_new_solver(const struct gkyl_gyrokinetic_multib *mbinp, int bid,
  const struct gkyl_gyrokinetic_multib_app *mbapp, struct gkyl_gyrokinetic_app *app)
{
  // For kinetic simulations, block dimension defined configuration-space dimensionality.
  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  const struct gkyl_gk_block_geom_info *bgi =
    gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, bid);

  // Construct top-level single-block input struct.
  struct gkyl_gk app_inp = { };

  // Set the configuration-space extents, cells.
  app_inp.cdim = cdim;
  for (int i=0; i<cdim; ++i) {
    app_inp.lower[i] = bgi->lower[i];
    app_inp.upper[i] = bgi->upper[i];
    app_inp.cells[i] = bgi->cells[i];
  }

  int num_species = app_inp.num_species = mbinp->num_species;
  int num_neut_species = app_inp.num_neut_species = mbinp->num_neut_species; 

  app_inp.poly_order = mbinp->poly_order;
  app_inp.basis_type = mbinp->basis_type;
  app_inp.cfl_frac = mbinp->cfl_frac; 
  app_inp.cfl_frac_omegaH = mbinp->cfl_frac_omegaH; 
  app_inp.eirene = mbinp->eirene;

  for (int i=0; i<num_species; ++i) {
    const struct gkyl_gyrokinetic_multib_species *sp = &mbinp->species[i];
    
    struct gkyl_gyrokinetic_species species_inp = { };
    strcpy(species_inp.name, sp->name);

    species_inp.charge = sp->charge;
    species_inp.mass = sp->mass;

    // Velocity-space information
    int vdim = species_inp.vdim = sp->vdim;
    for (int v=0; v<vdim; ++v) {
      species_inp.lower[v] = sp->lower[v];
      species_inp.upper[v] = sp->upper[v];
      species_inp.cells[v] = sp->cells[v];      
    }
    species_inp.mapc2p = sp->mapc2p;

    // Species physics modules.
    species_inp.collisionless = sp->collisionless;
    species_inp.collisions = sp->collisions;
    species_inp.anomalous_diffusion = sp->anomalous_diffusion;
    species_inp.radiation = sp->radiation;
    species_inp.react = sp->react;
    species_inp.react_neut = sp->react_neut; 

    // Species diagnostics
    species_inp.num_diag_moments = sp->num_diag_moments;
    for (int n=0; n<species_inp.num_diag_moments; ++n) {
      species_inp.diag_moments[n] = sp->diag_moments[n];
    }
    species_inp.num_integrated_diag_moments = sp->num_integrated_diag_moments;
    for (int n=0; n<species_inp.num_integrated_diag_moments; ++n) {
      species_inp.integrated_diag_moments[n] = sp->integrated_diag_moments[n];
    }
    species_inp.time_rate_diagnostics = sp->time_rate_diagnostics;
    species_inp.boundary_flux_diagnostics = sp->boundary_flux_diagnostics;

    // Choose proper block-specific species input.
    const struct gkyl_gyrokinetic_multib_species_pb *sp_pb = &sp->blocks[0];
    if (!sp->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i) {
        if (bid == sp->blocks[i].block_id) {
          sp_pb = &sp->blocks[i];
          break;
        }
      }
    }
    species_inp.projection = sp_pb->projection;
    species_inp.source = sp_pb->source;
    species_inp.polarization_density = sp_pb->polarization_density;  

    // By default, skip BCs altogether.
    for (int i=0; i<2*GKYL_MAX_CDIM; i++) {
      species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
      species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
    }

    int pardir = cdim-1;
    int num_below = gkyl_multib_conn_get_num_connected(mbapp->block_topo,
      bid, pardir, 0, GKYL_CONN_BELOW);
    int num_above = gkyl_multib_conn_get_num_connected(mbapp->block_topo,
      bid, pardir, 0, GKYL_CONN_ABOVE);

    int bc_count_sp[num_blocks];
    for (int i=0; i<num_blocks; i++)
      bc_count_sp[i] = 0;

    if (cdim ==3) {
      if (num_below > 0) {
        species_inp.bcs[bc_count_sp[bid]].dir = pardir;
        species_inp.bcs[bc_count_sp[bid]].edge = GKYL_LOWER_EDGE;
        species_inp.bcs[bc_count_sp[bid]].type = GKYL_BC_GK_SPECIES_TWISTSHIFT;
        bc_count_sp[bid] += 1;
      }

      if (num_above > 0) {
        species_inp.bcs[bc_count_sp[bid]].dir = pardir;
        species_inp.bcs[bc_count_sp[bid]].edge = GKYL_UPPER_EDGE;
        species_inp.bcs[bc_count_sp[bid]].type = GKYL_BC_GK_SPECIES_TWISTSHIFT;
        bc_count_sp[bid] += 1;
      }
    }

    // Set species physical BCs.
    for (int i=0; i<sp->num_physical_bcs; ++i) {
      if (bid == sp->bcs[i].bidx) {
        species_inp.bcs[bc_count_sp[bid]].dir = sp->bcs[i].dir;
        species_inp.bcs[bc_count_sp[bid]].edge = sp->bcs[i].edge;
        species_inp.bcs[bc_count_sp[bid]].type = sp->bcs[i].type;
        species_inp.bcs[bc_count_sp[bid]].aux_profile = sp->bcs[i].aux_profile;
        species_inp.bcs[bc_count_sp[bid]].aux_ctx = sp->bcs[i].aux_ctx;
        species_inp.bcs[bc_count_sp[bid]].projection = sp->bcs[i].projection;
        for (int k=0; k<3; ++k)
          species_inp.bcs[bc_count_sp[bid]].value[k] = sp->bcs[i].value[k];

        bc_count_sp[bid] += 1;
      }
    }

    // Copy species input into app input.
    memcpy(&app_inp.species[i], &species_inp, sizeof(struct gkyl_gyrokinetic_species));
  }

  for (int i=0; i<num_neut_species; ++i) {
    const struct gkyl_gyrokinetic_multib_neut_species *nsp = &mbinp->neut_species[i];
    
    struct gkyl_gyrokinetic_neut_species neut_species_inp = { };
    strcpy(neut_species_inp.name, nsp->name);

    neut_species_inp.mass = nsp->mass; 
    neut_species_inp.is_static = nsp->is_static; 

    // Velocity space information (neutrals are 3V)
    int vdim = neut_species_inp.vdim = nsp->vdim;
    for (int v=0; v<vdim; ++v) {
      neut_species_inp.lower[v] = nsp->lower[v];
      neut_species_inp.upper[v] = nsp->upper[v];
      neut_species_inp.cells[v] = nsp->cells[v];
    }
    neut_species_inp.mapc2p = nsp->mapc2p;

    // Neutral species physics modules
    neut_species_inp.react_neut = nsp->react_neut;

    // Neutral species diagnostics
    neut_species_inp.num_diag_moments = nsp->num_diag_moments;
    for (int n=0; n<neut_species_inp.num_diag_moments; ++n) {
      neut_species_inp.diag_moments[n] = nsp->diag_moments[n];
    }

    // choose proper block-specific species input
    const struct gkyl_gyrokinetic_multib_neut_species_pb *nsp_pb = &nsp->blocks[0];
    if (!nsp->duplicate_across_blocks) {
      for (int i=0; i<num_blocks; ++i) {
        if (bid == nsp->blocks[i].block_id) {
          nsp_pb = &nsp->blocks[i];
          break;
        }
      }
    }
    neut_species_inp.projection = nsp_pb->projection;
    neut_species_inp.source = nsp_pb->source;

    // By default, skip BCs altogether.
    for (int i=0; i<2*GKYL_MAX_CDIM; i++) {
      neut_species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
      neut_species_inp.bcs[i].type = GKYL_BC_GK_SKIP;
    }

    // Set species physical BCs.
    int bc_count_nsp[num_blocks];
    for (int i=0; i<num_blocks; i++)
      bc_count_nsp[i] = 0;

    for (int i=0; i<nsp->num_physical_bcs; ++i) {
      if (bid == nsp->bcs[i].bidx) {
        neut_species_inp.bcs[bc_count_nsp[bid]].dir = nsp->bcs[i].dir;
        neut_species_inp.bcs[bc_count_nsp[bid]].edge = nsp->bcs[i].edge;
        neut_species_inp.bcs[bc_count_nsp[bid]].type = nsp->bcs[i].type;
        neut_species_inp.bcs[bc_count_nsp[bid]].aux_profile = nsp->bcs[i].aux_profile;
        neut_species_inp.bcs[bc_count_nsp[bid]].aux_ctx = nsp->bcs[i].aux_ctx;
        neut_species_inp.bcs[bc_count_nsp[bid]].projection = nsp->bcs[i].projection;
        for (int k=0; k<3; ++k)
          neut_species_inp.bcs[bc_count_nsp[bid]].value[k] = nsp->bcs[i].value[k];

        bc_count_nsp[bid] += 1;
      }
    }

    // Copy neutral species input into app input.
    memcpy(&app_inp.neut_species[i], &neut_species_inp, sizeof(struct gkyl_gyrokinetic_neut_species));
  } 

  // Initialize the single-block field solver (only used for num_blocks=1).
  const struct gkyl_gyrokinetic_multib_field *fld = &mbinp->field;
  struct gkyl_gyrokinetic_field field_inp = { };
  field_inp.gkfield_id = fld->gkfield_id;
  field_inp.kperpSq = fld->kperpSq; 
  field_inp.time_rate_diagnostics = fld->time_rate_diagnostics; 
  
  // Adiabatic electron inputs.
  field_inp.electron_mass = fld->electron_mass;
  field_inp.electron_charge = fld->electron_charge;
  field_inp.electron_density = fld->electron_density; 
  field_inp.electron_temp = fld->electron_temp; 
  
  // BCs.
  for (int d=0; d<cdim; d++) {
    for (int i=0; i<2; i++) {
      // Set it to Dirichlet first, reset below. This avoids problems in
      // creating the single block field solve, which may not be used.
      field_inp.poisson_bcs[2*d+i].dir = d;
      field_inp.poisson_bcs[2*d+i].edge = i==0? GKYL_LOWER_EDGE : GKYL_UPPER_EDGE;
      field_inp.poisson_bcs[2*d+i].type = GKYL_BC_GK_FIELD_DIRICHLET;
      for (int k=0; k<3; ++k)
        field_inp.poisson_bcs[i].value[k] = -1.0e3;
    }
  }
  
  for (int i=0; i<fld->num_physical_bcs; i++) { 
    if (bid == fld->bcs[i].bidx) {
      struct gkyl_gyrokinetic_bc *bc_curr = gk_fetch_bc_with_dir_edge(field_inp.poisson_bcs, 2*cdim, fld->bcs[i].dir, fld->bcs[i].edge);
      bc_curr->type = fld->bcs[i].type;
      bc_curr->aux_profile = fld->bcs[i].aux_profile;
      bc_curr->aux_ctx = fld->bcs[i].aux_ctx;
      for (int k=0; k<3; ++k)
        bc_curr->value[k] = fld->bcs[i].value[k];
    }
  }
  
  const struct gkyl_gyrokinetic_multib_field_pb *fld_pb = &fld->blocks[0];
  // Choose proper block-specific field input.
  if (!fld->duplicate_across_blocks) {
    for (int i=0; i<num_blocks; ++i) {
      if (bid == fld->blocks[i].block_id) {
        const struct gkyl_gyrokinetic_multib_field_pb *fld_pb = &fld->blocks[i];
        break;
      }
    }
  }
  
  if (!fld->duplicate_across_blocks) {
    for (int i=0; i<num_blocks; ++i) {
      if (bid == fld->blocks[i].block_id) {
        fld_pb = &fld->blocks[i];
        break;
      }
    }
  }
  
  field_inp.polarization_bmag = fld_pb->polarization_bmag ? fld_pb->polarization_bmag : mbapp->bmag_ref;
  field_inp.kperpSq = fld_pb->kperpSq;
  field_inp.time_rate_diagnostics = fld_pb->time_rate_diagnostics; 
  
  field_inp.phi_wall_lo_ctx = fld_pb->phi_wall_lo_ctx; 
  field_inp.phi_wall_lo = fld_pb->phi_wall_lo; 
  field_inp.phi_wall_lo_evolve = fld_pb->phi_wall_lo_evolve; 
  
  field_inp.phi_wall_up_ctx = fld_pb->phi_wall_up_ctx; 
  field_inp.phi_wall_up = fld_pb->phi_wall_up; 
  field_inp.phi_wall_up_evolve = fld_pb->phi_wall_up_evolve;   
  

  // Copy field input into app input.
  memcpy(&app_inp.field, &field_inp, sizeof(struct gkyl_gyrokinetic_field));  

  gkyl_gyrokinetic_app_new_solver(&app_inp, app);
}

gkyl_gyrokinetic_multib_app*
gkyl_gyrokinetic_multib_app_new_geom(const struct gkyl_gyrokinetic_multib *mbinp)
{
  int my_rank, num_ranks;
  gkyl_comm_get_rank(mbinp->comm, &my_rank);
  gkyl_comm_get_size(mbinp->comm, &num_ranks);

  int tot_max[2];
  calc_tot_and_max_cuts(mbinp->gk_block_geom, tot_max);
  if ((num_ranks > tot_max[0]) || (num_ranks < tot_max[1])) {
    fprintf(stderr, "\nSpecified %d total cuts but provided %d processes, \
and the maximum number of cuts in a block is %d\n\n", tot_max[0], num_ranks, tot_max[1]);
    return 0;
  }

  struct gkyl_gyrokinetic_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));

  strcpy(mbapp->name, mbinp->name);
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);  
  mbapp->use_gpu = mbinp->use_gpu;
  
  mbapp->gk_block_geom = gkyl_gk_block_geom_acquire(mbinp->gk_block_geom);
  mbapp->block_topo = gkyl_gk_block_geom_topo(mbinp->gk_block_geom);
  
  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  // Construct round-robin decomposition.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    branks[i] = calc_cuts(cdim, bgi->cuts);
  }
  mbapp->round_robin = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int num_local_blocks = 0;
  mbapp->local_blocks = gkyl_malloc(sizeof(int[num_blocks]));

  int lidx = 0;
  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // This is larger than needed.

  mbapp->decomp = gkyl_malloc(num_blocks*sizeof(struct gkyl_rect_decomp*));
  
  // Construct list of block communicators: there are as many
  // communicators as blocks. Not all communicators are valid on each
  // rank. The total number of valid communicators is
  // num_local_blocks.
  mbapp->block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(mbapp->round_robin, i, rank_list);

    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);

    if (is_my_rank_in_decomp) {
      mbapp->local_blocks[lidx++] = i;
      num_local_blocks += 1;      
    }

    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    struct gkyl_range block_global_range;
    gkyl_create_global_range(cdim, bgi->cells, &block_global_range);

    mbapp->decomp[i] = gkyl_rect_decomp_new_from_cuts(
      cdim, bgi->cuts, &block_global_range);

    bool status;
    mbapp->block_comms[i] = gkyl_comm_create_comm_from_ranks(mbinp->comm,
      branks[i], rank_list, mbapp->decomp[i], &status);
  }

  mbapp->num_local_blocks = num_local_blocks;  

  printf("Rank %d handles %d Apps\n", my_rank, num_local_blocks);
  for (int i=0; i<num_local_blocks; ++i)
    printf("  Rank %d handles block %d\n", my_rank, mbapp->local_blocks[i]);

  mbapp->num_species = 0;
  mbapp->num_neut_species = 0;
  mbapp->update_field = 0;
  mbapp->singleb_apps = 0;

  if (num_local_blocks > 0) {
    mbapp->num_species = mbinp->num_species;
    mbapp->num_neut_species = mbinp->num_neut_species;
    mbapp->update_field = !mbinp->skip_field; // Note inversion of truth value (default: update field).


    mbapp->singleb_apps = gkyl_malloc(num_local_blocks*sizeof(struct gkyl_gyrokinetic_app*));
  }

  for (int i=0; i<mbinp->num_species; ++i)
    strcpy(mbapp->species_name[i], mbinp->species[i].name);

  for (int i=0; i<mbinp->num_neut_species; ++i)
    strcpy(mbapp->neut_species_name[i], mbinp->neut_species[i].name);  

  // Write the block topo file.
  gkyl_gyrokinetic_multib_app_write_topo(mbapp);

  // Create single-block grids and geometries.
  for (int i=0; i<num_local_blocks; ++i)
    mbapp->singleb_apps[i] = singleb_app_new_geom(mbinp, mbapp->local_blocks[i], mbapp);


  return mbapp;
}

gkyl_gyrokinetic_multib_app* gkyl_gyrokinetic_multib_app_new(const struct gkyl_gyrokinetic_multib *mbinp)
{
  int my_rank, num_ranks;
  gkyl_comm_get_rank(mbinp->comm, &my_rank);
  gkyl_comm_get_size(mbinp->comm, &num_ranks);

  int tot_max[2];
  calc_tot_and_max_cuts(mbinp->gk_block_geom, tot_max);
  if ((num_ranks > tot_max[0]) || (num_ranks < tot_max[1])) {
    fprintf(stderr, "\nSpecified %d total cuts but provided %d processes, \
and the maximum number of cuts in a block is %d\n\n", tot_max[0], num_ranks, tot_max[1]);
    return 0;
  }

  struct gkyl_gyrokinetic_multib_app *mbapp = gkyl_malloc(sizeof(*mbapp));

  strcpy(mbapp->name, mbinp->name);
  mbapp->comm = gkyl_comm_acquire(mbinp->comm);  
  mbapp->use_gpu = mbinp->use_gpu;
  
  mbapp->gk_block_geom = gkyl_gk_block_geom_acquire(mbinp->gk_block_geom);
  mbapp->block_topo = gkyl_gk_block_geom_topo(mbinp->gk_block_geom);
  
  int cdim = gkyl_gk_block_geom_ndim(mbapp->gk_block_geom);
  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  // Construct round-robin decomposition.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    branks[i] = calc_cuts(cdim, bgi->cuts);
  }
  mbapp->round_robin = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int num_local_blocks = 0;
  mbapp->local_blocks = gkyl_malloc(sizeof(int[num_blocks]));

  int lidx = 0;
  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // This is larger than needed.

  mbapp->decomp = gkyl_malloc(num_blocks*sizeof(struct gkyl_rect_decomp*));
  
  // Construct list of block communicators: there are as many
  // communicators as blocks. Not all communicators are valid on each
  // rank. The total number of valid communicators is
  // num_local_blocks.
  mbapp->block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(mbapp->round_robin, i, rank_list);

    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);

    if (is_my_rank_in_decomp) {
      mbapp->local_blocks[lidx++] = i;
      num_local_blocks += 1;      
    }

    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    struct gkyl_range block_global_range;
    gkyl_create_global_range(cdim, bgi->cells, &block_global_range);

    mbapp->decomp[i] = gkyl_rect_decomp_new_from_cuts(
      cdim, bgi->cuts, &block_global_range);

    bool status;
    mbapp->block_comms[i] = gkyl_comm_create_comm_from_ranks(mbinp->comm,
      branks[i], rank_list, mbapp->decomp[i], &status);
  }

  mbapp->num_local_blocks = num_local_blocks;  

  // Create multiblock metadata to pass to each block, and combine with user's metadata.
  const char *fmt_btopo = "%s-block_topo.gkyl";
  int sz = gkyl_calc_strlen(fmt_btopo, mbapp->name);
  char fileNm_btopo[sz+1]; // ensures no buffer overflow
  snprintf(fileNm_btopo, sizeof fileNm_btopo, fmt_btopo, mbapp->name);

  // Basic metadata for I/O (including metadata optional from user).
  const char* build_id = GIT_COMMIT_ID;
  const char* build_date = GKYL_BUILD_DATE;
  struct gkyl_msgpack_map_elem io_meta_default[] = {
    { .key = "changeset", .elem_type = GKYL_MP_STRING, .cval = (char *)build_id },
    { .key = "builddate", .elem_type = GKYL_MP_STRING, .cval = (char *)build_date },
    { .key = "is_multib", .elem_type = GKYL_MP_UNSIGNED_INT, .uval = 1 },
    { .key = "topo_file", .elem_type = GKYL_MP_STRING, .cval = fileNm_btopo },
  };
  const struct gkyl_msgpack_map_elem *io_meta_union[] = {io_meta_default, mbinp->metadata.attributes};
  int io_meta_union_len[] = {sizeof(io_meta_default)/sizeof(io_meta_default[0]), mbinp->metadata.num_attributes};

  mbapp->io_meta_basic = gkyl_msgpack_map_elem_union(sizeof(io_meta_union)/sizeof(io_meta_union[0]),
    io_meta_union_len, io_meta_union, &mbapp->io_meta_basic_len);

  // Write the block topo file.
  gkyl_gyrokinetic_multib_app_write_topo(mbapp);

  printf("Rank %d handles %d Apps\n", my_rank, num_local_blocks);
  for (int i=0; i<num_local_blocks; ++i)
    printf("  Rank %d handles block %d\n", my_rank, mbapp->local_blocks[i]);

  mbapp->num_species = 0;
  mbapp->num_neut_species = 0;
  mbapp->update_field = 0;
  mbapp->singleb_apps = 0;

  if (num_local_blocks > 0) {
    mbapp->num_species = mbinp->num_species;
    mbapp->num_neut_species = mbinp->num_neut_species;
    mbapp->update_field = !mbinp->skip_field; // Note inversion of truth value (default: update field).

    mbapp->singleb_apps = gkyl_malloc(num_local_blocks*sizeof(struct gkyl_gyrokinetic_app*));
  }

  for (int i=0; i<mbinp->num_species; ++i)
    strcpy(mbapp->species_name[i], mbinp->species[i].name);

  for (int i=0; i<mbinp->num_neut_species; ++i)
    strcpy(mbapp->neut_species_name[i], mbinp->neut_species[i].name);  

  // Create single-block grids and geometries.
  for (int i=0; i<num_local_blocks; ++i)
    mbapp->singleb_apps[i] = singleb_app_new_geom(mbinp, mbapp->local_blocks[i], mbapp);

  // Set bmag_ref.
  double bmag_min_local = DBL_MAX;
  double bmag_min_global;
  for (int i=0; i<num_local_blocks; ++i) {
    double bmag_min = gkyl_gk_geometry_reduce_bmag(mbapp->singleb_apps[i]->gk_geom, GKYL_MIN);
    bmag_min_local = GKYL_MIN2(bmag_min_local, bmag_min);
  }
  gkyl_comm_allreduce_host(mbapp->comm, GKYL_DOUBLE, GKYL_MIN, 1, &bmag_min_local, &bmag_min_global);

  double bmag_max_local = -DBL_MAX;
  double bmag_max_global;
  for (int i=0; i<num_local_blocks; ++i) {
    double bmag_max = gkyl_gk_geometry_reduce_bmag(mbapp->singleb_apps[i]->gk_geom, GKYL_MAX);
    bmag_max_local = GKYL_MAX2(bmag_max_local, bmag_max);
  }
  gkyl_comm_allreduce_host(mbapp->comm, GKYL_DOUBLE, GKYL_MAX, 1, &bmag_max_local, &bmag_max_global);

  mbapp->bmag_ref = (bmag_max_global + bmag_min_global)/2.0;

  // Create connections needed for conf-space syncs.
  int ghost[] = { 1, 1, 1 };
  mbapp->mbcc_sync_conf = gkyl_malloc(sizeof(struct gkyl_mbcc_sr));
  mbapp->mbcc_sync_conf->send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  mbapp->mbcc_sync_conf->recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));

  for (int bI=0; bI<num_local_blocks; ++bI) {
    int bid = mbapp->local_blocks[bI];

    gkyl_rrobin_decomp_getranks(mbapp->round_robin, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbapp->mbcc_sync_conf->recv[bI] = gkyl_multib_comm_conn_new_recv(bid, brank, 
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);
    mbapp->mbcc_sync_conf->send[bI] = gkyl_multib_comm_conn_new_send(bid, brank,
      ghost, &mbapp->block_topo->conn[bid], mbapp->decomp);

    struct gkyl_multib_comm_conn *mbcc_s = mbapp->mbcc_sync_conf->send[bI], *mbcc_r = mbapp->mbcc_sync_conf->recv[bI];
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];

    for (int ns=0; ns<mbcc_s->num_comm_conn; ++ns) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      struct gkyl_comm_conn *ccs = &mbcc_s->comm_conn[ns];
      int rankIdx = ccs->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccs->block_id, rank_list);
      ccs->rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&ccs->range, &sbapp->local_ext, ccs->range.lower, ccs->range.upper);
    }
    for (int nr=0; nr<mbcc_r->num_comm_conn; ++nr) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      struct gkyl_comm_conn *ccr = &mbcc_r->comm_conn[nr];
      int rankIdx = ccr->rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, ccr->block_id, rank_list);
      ccr->rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&ccr->range, &sbapp->local_ext, ccr->range.lower, ccr->range.upper);
    }

    // Sort connections according to rank and block ID (needed by NCCL).
    gkyl_multib_comm_conn_sort(mbcc_r);
    gkyl_multib_comm_conn_sort(mbcc_s);
  }

  // Sync the conf-space volume Jacobian needed for syncing quantities that include a
  // jacobgeo factor in them.
  struct gkyl_array *jacs_vol[mbapp->num_local_blocks];
  for (int b=0; b<mbapp->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
    jacs_vol[b] = sbapp->gk_geom->geo_int.jacobgeo_ghost;
  }
  // Sync across blocks.
  gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
    mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs_vol, jacs_vol);

  // Sync the surface conf-space Jacobian, compute its reciprocal, and 
  // store its product with the Jacobian of this block (in the ghost cell).
  for (int d = 0; d<cdim; d++) {
    // Sync jacobgeo.
    struct gkyl_array *jacs[mbapp->num_local_blocks];
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
      jacs[b] = geo_surf.jacobgeo_ratio;
      gkyl_array_copy_range(jacs[b], geo_surf.jacobgeo, &sbapp->local_lower_skin[d]);
      gkyl_array_copy_range_to_range(jacs[b], geo_surf.jacobgeo, &sbapp->local_upper_skin[d], &sbapp->local_upper_ghost[d]);
    }
    gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
      mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs, jacs);

    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
      struct gkyl_array *jacgeo = mkarr(mbapp->use_gpu, geo_surf.jacobgeo_ratio->ncomp, geo_surf.jacobgeo_ratio->size);

      // Compute 1/jacobgeo in ghost cells.
      gkyl_array_set_range(jacgeo, 1.0, geo_surf.jacobgeo_ratio, &sbapp->local_lower_ghost[d]);
      gkyl_array_set_range(jacgeo, 1.0, geo_surf.jacobgeo_ratio, &sbapp->local_upper_ghost[d]);
      gkyl_dg_inv_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio, 0, jacgeo, &sbapp->local_lower_ghost[d]);
      gkyl_dg_inv_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio, 0, jacgeo, &sbapp->local_upper_ghost[d]);
      // Multiply by the Jacobian of this block.
      gkyl_array_copy_range_to_range(jacgeo, geo_surf.jacobgeo, &sbapp->local_lower_ghost[d], &sbapp->local_lower_skin[d]);
      gkyl_array_copy_range(jacgeo, geo_surf.jacobgeo, &sbapp->local_upper_ghost[d]);
      gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio,
        0, jacgeo, 0, geo_surf.jacobgeo_ratio, &sbapp->local_lower_ghost[d]);
      gkyl_dg_mul_op_range(&sbapp->gk_geom->surf_basis, 0, geo_surf.jacobgeo_ratio,
        0, jacgeo, 0, geo_surf.jacobgeo_ratio, &sbapp->local_upper_ghost[d]);
      // Set the ratio to 1 in the interior (shouldn't be in use).
      gkyl_array_clear_range(geo_surf.jacobgeo_ratio, 0.0, &sbapp->local);
      gkyl_array_shiftc_range(geo_surf.jacobgeo_ratio, pow(sqrt(2.0),sbapp->cdim), 0, &sbapp->local);

      gkyl_array_release(jacgeo);
    }

  }

  const struct gkyl_gk_block_geom_info *bgi0 = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, 0);

  if (cdim > 1 && bgi0->geometry.geometry_id == GKYL_GEOMETRY_TOKAMAK) {
    // Sync the surface deltats. Need to copy the upper ghost into the upper skin before syncing.
    for (int d = 0; d<cdim; d++) {
      struct gkyl_array *deltats[mbapp->num_local_blocks];
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        struct gk_geom_surf geo_surf = sbapp->gk_geom->geo_surf[d];
        deltats[b] = geo_surf.deltats;
        gkyl_array_copy_range_to_range(deltats[b], deltats[b], &sbapp->local_upper_skin[d], &sbapp->local_upper_ghost[d]);
      }
      gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
        mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, deltats, deltats);
    }
    // Accumulate the appropriate shift
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      int par_dir = cdim-1;
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      struct gkyl_array *delta_ts = sbapp->gk_geom->geo_surf[par_dir].deltats;
      struct gkyl_array *buffer = mkarr(sbapp->use_gpu, delta_ts->ncomp, delta_ts->size);

      gkyl_array_copy_range_to_range(buffer, delta_ts, &sbapp->local_upper_skin[par_dir], &sbapp->local_upper_ghost[par_dir]);
      gkyl_array_accumulate_range(delta_ts, -1.0, buffer, &sbapp->local_upper_skin[par_dir]);

      gkyl_array_copy_range_to_range(buffer, delta_ts, &sbapp->local_lower_skin[par_dir], &sbapp->local_lower_ghost[par_dir]);
      gkyl_array_accumulate_range(delta_ts, -1.0, buffer, &sbapp->local_lower_skin[par_dir]);

      gkyl_array_release(buffer);

      // Deflate delta_ts.
      gyrokinetic_deflate_delta_ts(sbapp, delta_ts);

      // Write the twistshift shift.
      gyrokinetic_app_write_ts_shift(sbapp);
    }
  }


  // Create the rest of the single-block solvers.
  for (int i=0; i<num_local_blocks; ++i)
    singleb_app_new_solver(mbinp, mbapp->local_blocks[i], mbapp, mbapp->singleb_apps[i]);

  // Create the MB field app.
  mbapp->field = gk_multib_field_new(mbinp, mbapp);


  // Create connections needed for syncing charged species phase-space quantities.
  mbapp->mbcc_sync_charged = gkyl_malloc(mbapp->num_species * sizeof(struct gkyl_mbcc_sr));
  for (int i=0; i<mbinp->num_species; ++i) {
    mbapp->mbcc_sync_charged[i].send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
    mbapp->mbcc_sync_charged[i].recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  }
  for (int bI=0; bI<num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    struct gkyl_multib_comm_conn *mbcc_cs = mbapp->mbcc_sync_conf->send[bI], *mbcc_cr = mbapp->mbcc_sync_conf->recv[bI];
    for (int i=0; i<mbinp->num_species; ++i) {
      mbapp->mbcc_sync_charged[i].send[bI] = gkyl_multib_comm_conn_new(mbcc_cs->num_comm_conn, mbcc_cs->comm_conn);
      mbapp->mbcc_sync_charged[i].recv[bI] = gkyl_multib_comm_conn_new(mbcc_cr->num_comm_conn, mbcc_cr->comm_conn);
      struct gkyl_multib_comm_conn *mbcc_ps = mbapp->mbcc_sync_charged[i].send[bI],
                                   *mbcc_pr = mbapp->mbcc_sync_charged[i].recv[bI];
      // Extend ranges to include velocity space, and make them a sub range.
      struct gk_species *gks = &sbapp->species[i];
      for (int ns=0; ns<mbcc_cs->num_comm_conn; ++ns) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cs->comm_conn[ns];
        struct gkyl_comm_conn *ccs_phase = &mbcc_ps->comm_conn[ns];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gks->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gks->local_ext, phase_r.lower, phase_r.upper);
      }
      for (int nr=0; nr<mbcc_cr->num_comm_conn; ++nr) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cr->comm_conn[nr];
        struct gkyl_comm_conn *ccs_phase = &mbcc_pr->comm_conn[nr];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gks->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gks->local_ext, phase_r.lower, phase_r.upper);
      }
    }
  }

  // Create connections needed for syncing neutral species phase-space quantities.
  mbapp->mbcc_sync_neut = gkyl_malloc(mbapp->num_neut_species * sizeof(struct gkyl_mbcc_sr));
  for (int i=0; i<mbinp->num_neut_species; ++i) {
    mbapp->mbcc_sync_neut[i].send = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
    mbapp->mbcc_sync_neut[i].recv = gkyl_malloc(mbapp->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  }
  for (int bI=0; bI<num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    struct gkyl_multib_comm_conn *mbcc_cs = mbapp->mbcc_sync_conf->send[bI], *mbcc_cr = mbapp->mbcc_sync_conf->recv[bI];
    for (int i=0; i<mbinp->num_neut_species; ++i) {
      mbapp->mbcc_sync_neut[i].send[bI] = gkyl_multib_comm_conn_new(mbcc_cs->num_comm_conn, mbcc_cs->comm_conn);
      mbapp->mbcc_sync_neut[i].recv[bI] = gkyl_multib_comm_conn_new(mbcc_cr->num_comm_conn, mbcc_cr->comm_conn);
      struct gkyl_multib_comm_conn *mbcc_ps = mbapp->mbcc_sync_neut[i].send[bI],
                                   *mbcc_pr = mbapp->mbcc_sync_neut[i].recv[bI];
      // Extend ranges to include velocity space, and make them a sub range.
      struct gk_neut_species *gkns = &sbapp->neut_species[i];
      for (int ns=0; ns<mbcc_cs->num_comm_conn; ++ns) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cs->comm_conn[ns];
        struct gkyl_comm_conn *ccs_phase = &mbcc_ps->comm_conn[ns];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gkns->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gkns->local_ext, phase_r.lower, phase_r.upper);
      }
      for (int nr=0; nr<mbcc_cr->num_comm_conn; ++nr) {
        struct gkyl_comm_conn *ccs_conf = &mbcc_cr->comm_conn[nr];
        struct gkyl_comm_conn *ccs_phase = &mbcc_pr->comm_conn[nr];
        struct gkyl_range phase_r;
        gkyl_range_ten_prod(&phase_r, &ccs_conf->range, &gkns->local_vel);
        gkyl_sub_range_init(&ccs_phase->range, &gkns->local_ext, phase_r.lower, phase_r.upper);
      }
    }
  }

  // Sync the effective diffusivity of the anomalous diffusion operator.
  // Assume they either all have anomalous diffusion or none of them do.
  bool any_anomalous_diff = false;
  for (int i=0; i<mbapp->num_species; ++i) {
    bool has_anomalous_diff = true;
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      has_anomalous_diff = has_anomalous_diff && sbapp->species[i].anom_diff.anom_diff_id;
    }

    if (has_anomalous_diff) {
      any_anomalous_diff = true;
      // Divide diffD by dz before the transfer
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        int d  = 0;
        gkyl_array_scale_range(sbapp->species[i].anom_diff.diffD, sbapp->grid.dx[sbapp->cdim-1], &sbapp->global);
      }

      // Sync
      struct gkyl_array *gkad_nu[mbapp->num_local_blocks];
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        gkad_nu[b] = sbapp->species[i].anom_diff.diffD;
      }
      gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
        mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, gkad_nu, gkad_nu);

      // Multiply by diffD dz after the transfer to achieve rescaling
      for (int b=0; b<mbapp->num_local_blocks; ++b) {
        struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
        int d  = 0;
        gkyl_array_scale_range(sbapp->species[i].anom_diff.diffD, 1.0/sbapp->grid.dx[sbapp->cdim-1], &sbapp->global_ext);
      }

    }
  }

  if (any_anomalous_diff) {
    // Sync the interior conf-space reciprocal Jacobian.
    struct gkyl_array *jacs_inv_vol[mbapp->num_local_blocks];
    for (int b=0; b<mbapp->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[b];
      jacs_inv_vol[b] = sbapp->gk_geom->geo_int.jacobgeo_inv;
    }
    gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
      mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, jacs_inv_vol, jacs_inv_vol);
  }

  mbapp->stat = (struct gkyl_gyrokinetic_stat) {};

  mbapp->dts = gkyl_dynvec_new(GKYL_DOUBLE, 1); // Dynvector to store time steps.
  mbapp->is_first_dt_write_call = true;

  gkyl_free(rank_list);
  gkyl_free(branks);

  return mbapp;
}

void
gyrokinetic_multib_calc_field(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  struct timespec wtm = gkyl_wall_clock();
  // Compute fields.
  if (app->update_field) {
    // Solve the field equation.
    gk_multib_field_rhs(app, app->field, fin, bflux);
  }
  app->stat.field_tm += gkyl_time_diff_now_sec(wtm);
}

static void
gyrokinetic_multib_apply_bc(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array *distf_neut[])
{
  // Apply boundary conditions in each block (including intrablock sync).
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    for (int i=0; i<app->num_species; ++i) {
      gk_species_apply_bc(sbapp, &sbapp->species[i], distf[li_charged+i]);
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      gk_neut_species_apply_bc(sbapp, &sbapp->neut_species[i], distf_neut[li_neut+i]);
    }
  }

  // Sync blocks.
  struct timespec wst = gkyl_wall_clock();
  int cdim = gkyl_gk_block_geom_ndim(app->gk_block_geom);
  for (int i=0; i<app->num_species; ++i) {
    // Sync charged species.
    struct gkyl_array *fs[app->num_local_blocks];
    for (int b=0; b<app->num_local_blocks; ++b) {
      int li_charged = b * app->num_species;
      fs[b] = distf[li_charged+i];
    }
    gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
      app->mbcc_sync_charged[i].send, app->mbcc_sync_charged[i].recv, fs, fs);
  }
  app->stat.species_bc_tm += gkyl_time_diff_now_sec(wst);

  struct timespec wst_neut = gkyl_wall_clock();
  struct gkyl_gyrokinetic_app *sbapp0 = app->singleb_apps[0];
  for (int i=0; i<app->num_neut_species; ++i) {
    // Sync neutral species.
    if (!sbapp0->neut_species[i].info.is_static) {
      struct gkyl_array *fs[app->num_local_blocks];
      for (int b=0; b<app->num_local_blocks; ++b) {
        int li_neut = b * app->num_neut_species;
        fs[b] = distf_neut[li_neut+i];
      }
      gkyl_multib_comm_conn_array_transfer(app->comm, app->num_local_blocks, app->local_blocks,
        app->mbcc_sync_neut[i].send, app->mbcc_sync_neut[i].recv, fs, fs);
    }
  }
  app->stat.neut_species_bc_tm += gkyl_time_diff_now_sec(wst_neut);

}

void
gyrokinetic_multib_calc_field_and_apply_bc(struct gkyl_gyrokinetic_multib_app* app, double tcurr,
  struct gkyl_array *distf[], struct gkyl_array **bflux[], struct gkyl_array *distf_neut[])
{
  // Compute fields and apply BCs.

  // Compute the field.
  // MF 2024/09/27/: Need the cast here for consistency. Fixing
  // this may require removing 'const' from a lot of places.
  gyrokinetic_multib_calc_field(app, tcurr, (const struct gkyl_array **) distf, bflux);

  // Apply boundary conditions.
  struct timespec wst = gkyl_wall_clock();
  gyrokinetic_multib_apply_bc(app, tcurr, distf, distf_neut);
  app->stat.bc_tm += gkyl_time_diff_now_sec(wst);
}

void
gkyl_gyrokinetic_multib_app_apply_ic(gkyl_gyrokinetic_multib_app* app, double t0)
{
  app->tcurr = t0;
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_apply_ic_species(app, i, t0);
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_apply_ic_neut_species(app, i, t0);
  }  

  // Compute the fields and apply BCs.
  struct gkyl_array *distf[app->num_species * app->num_local_blocks];
  struct gkyl_array **bflux[app->num_species * app->num_local_blocks];
  struct gkyl_array *distf_neut[app->num_neut_species * app->num_local_blocks];
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    int li_charged = b * app->num_species;
    int li_neut = b * app->num_neut_species;
    for (int i=0; i<app->num_species; ++i) {
      distf[li_charged+i] = sbapp->species[i].f;
      bflux[li_charged+i] = sbapp->species[i].bflux.f1;
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      distf_neut[li_neut+i] = sbapp->neut_species[i].f;
    }
  }
  gyrokinetic_multib_calc_field_and_apply_bc(app, t0, distf, bflux, distf_neut);
}

void
gkyl_gyrokinetic_multib_app_apply_ic_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  app->tcurr = t0;
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_apply_ic_species(app->singleb_apps[b], sidx, t0);
  }
  gkyl_comm_barrier(app->comm);
}

void
gkyl_gyrokinetic_multib_app_apply_ic_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double t0)
{
  app->tcurr = t0;
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_apply_ic_neut_species(app->singleb_apps[b], sidx, t0);
  }
  gkyl_comm_barrier(app->comm);
}


struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_file_field(gkyl_gyrokinetic_multib_app *app, const char *fname)
{
 // TO DO
  return (struct gkyl_app_restart_status) { };

}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status 
gkyl_gyrokinetic_multib_app_from_file_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx,
  const char *fname)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_read_from_frame(gkyl_gyrokinetic_multib_app *app, int frame)
{
  struct gkyl_app_restart_status rstat;
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    for (int i=0; i<app->num_neut_species; i++) {
      int neut_frame = frame;
      if (sbapp->neut_species[i].info.is_static) {
        neut_frame = 0;
      }
      rstat = gkyl_gyrokinetic_app_from_frame_neut_species(sbapp, i, neut_frame);
    }
    for (int i=0; i<app->num_species; i++) {
      rstat = gkyl_gyrokinetic_app_from_frame_species(sbapp, i, frame);
    }
  }
  
  if (rstat.io_status == GKYL_ARRAY_RIO_SUCCESS) {
    // Compute the fields and apply BCs.
    struct gkyl_array *distf[app->num_species * app->num_local_blocks];
    struct gkyl_array **bflux[app->num_species * app->num_local_blocks];
    struct gkyl_array *distf_neut[app->num_neut_species * app->num_local_blocks];
    for (int b=0; b<app->num_local_blocks; ++b) {
      struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
      int li_charged = b * app->num_species;
      int li_neut = b * app->num_neut_species;
      for (int i=0; i<app->num_species; ++i) {
        distf[li_charged+i] = sbapp->species[i].f;
        bflux[li_charged+i] = sbapp->species[i].bflux.f;
      }
      for (int i=0; i<app->num_neut_species; ++i) {
        distf_neut[li_neut+i] = sbapp->neut_species[i].f;
      }
    }
//    if (app->update_field && app->field->gkfield_id == GKYL_GK_FIELD_BOLTZMANN) {
//      for (int i=0; i<app->num_species; ++i) {
//        struct gk_species *s = &app->species[i];
//
//        // Compute advection speeds so we can compute the initial boundary flux.
//        gkyl_gk_collisionless_flux_alpha_surf(s->calc_gk_vars, 
//          &app->local, &s->local, &s->local_ext, app->field->phi_smooth,
//          s->alpha_surf, s->sgn_alpha_surf, s->const_sgn_alpha);
//
//        // Compute and store (in the ghost cell of of out) the boundary fluxes.
//        // NOTE: this overwrites ghost cells that may be used for sourcing.
//        gk_species_bflux_rhs(app, s, &s->bflux, distf[i], distf[i]);
//      }
//    }
    gyrokinetic_multib_calc_field_and_apply_bc(app, rstat.stime, distf, bflux, distf_neut);
  }

  struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[0];
  sbapp->field->is_first_energy_write_call = false; // Append to existing diagnostic.
                                                    //
  app->is_first_dt_write_call = false;

  return rstat;
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_field(gkyl_gyrokinetic_multib_app *app, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

struct gkyl_app_restart_status
gkyl_gyrokinetic_multib_app_from_frame_neut_species(gkyl_gyrokinetic_multib_app *app, int sidx, int frame)
{
  // TO DO
  return (struct gkyl_app_restart_status) { };
}

// private function to handle variable argument list for printing
static void
v_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, va_list argp)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if ((rank == 0) && fp)
    vfprintf(fp, fmt, argp);
}

void
gkyl_gyrokinetic_multib_app_cout(const gkyl_gyrokinetic_multib_app* app, FILE *fp, const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  v_gyrokinetic_multib_app_cout(app, fp, fmt, argp);
  va_end(argp);
}

void
gkyl_gyrokinetic_multib_app_write_topo(const gkyl_gyrokinetic_multib_app* app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    cstr file_name = cstr_from_fmt("%s-block_topo.gkyl", app->name);
    gkyl_block_topo_write(app->block_topo, file_name.str);
    cstr_drop(&file_name);
  }
}

void
gkyl_gyrokinetic_multib_app_write_field(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_field(app->singleb_apps[b], tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_field_energy(gkyl_gyrokinetic_multib_app* app, double tm)
{
  if (app->update_field) {
    for (int b=0; b<app->num_local_blocks; ++b) {
      gkyl_gyrokinetic_app_calc_field_energy(app->singleb_apps[b], tm);
    }
  }
}

void
gkyl_gyrokinetic_multib_app_write_field_energy(gkyl_gyrokinetic_multib_app* app)
{
  if (app->update_field) {
    for (int b=0; b<app->num_local_blocks; ++b) {
      gkyl_gyrokinetic_app_write_field_energy(app->singleb_apps[b]);
    }
  }
}

void
gkyl_gyrokinetic_multib_app_write_eirene(gkyl_gyrokinetic_multib_app *app, double tm, int frame)
{

  for (int b=0; b<app->num_local_blocks; ++b)
    gkyl_gyrokinetic_app_write_eirene_diagnostics(app->singleb_apps[b], tm, frame);

  struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[0];
  cstr fileNm = cstr_from_fmt("%snew_data_flag", sbapp->eirene->info.output_data_path);
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);
  if (0 == rank) {
    FILE *fp = fopen(fileNm.str, "w");
    if (fp == NULL)
        return;
    fprintf(fp, "%d\n", frame);
    fclose(fp);
  }
  cstr_drop(&fileNm);
}

void
gkyl_gyrokinetic_multib_app_calc_eirene_integrated_diagnostics(gkyl_gyrokinetic_multib_app *app, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b)
    gkyl_gyrokinetic_app_calc_eirene_integrated_diagnostics(app->singleb_apps[b], tm);
}

void
gkyl_gyrokinetic_multib_app_write_eirene_integrated_diagnostics(gkyl_gyrokinetic_multib_app *app)
{
  for (int b=0; b<app->num_local_blocks; ++b)
    gkyl_gyrokinetic_app_write_eirene_integrated_diagnostics(app->singleb_apps[b]);
}

//
// ............. Species outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_neut_species_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_boundary_flux_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_integrated_mom(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_integrated_mom(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_boundary_flux_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_boundary_flux_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. Source outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_source_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_source_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_neut_species_source_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source_integrated_mom(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_source_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_source_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. BGK Source outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_calc_species_source_bgk_integrated_diagnostics(gkyl_gyrokinetic_multib_app* app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_source_bgk_integrated_diagnostics(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_source_bgk_integrated_diagnostics(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_source_bgk_integrated_diagnostics(app->singleb_apps[b], sidx);
  }
}

//
// ............. LTE outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_lte_max_corr_status(gkyl_gyrokinetic_multib_app* app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_lte_max_corr_status(app->singleb_apps[b], sidx);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_lte_max_corr_status(gkyl_gyrokinetic_multib_app* app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_lte_max_corr_status(app->singleb_apps[b], sidx);
  }
}

//
// ............. Collision outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_lbo_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_lbo_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_bgk_mom(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_bgk_mom(app->singleb_apps[b], sidx, tm, frame);
  }
}

//
// ............. Radiation outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_species_rad_drag(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_drag(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_rad_emissivity(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_emissivity(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx, double tm)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_calc_species_rad_integrated_mom(app->singleb_apps[b], sidx, tm);
  }
}

void
gkyl_gyrokinetic_multib_app_write_species_rad_integrated_mom(gkyl_gyrokinetic_multib_app *app, int sidx)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_rad_integrated_mom(app->singleb_apps[b], sidx);
  }
}

//
// ............. Functions that group several outputs for a single species ............... //
//
void
gkyl_gyrokinetic_multib_app_write_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_species_source(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_species_rad_drag(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_phase(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_neut_species(app, sidx, tm, frame);

  gkyl_gyrokinetic_multib_app_write_neut_species_source(app, sidx, tm, frame);
}

void
gkyl_gyrokinetic_multib_app_write_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_species_conf(app->singleb_apps[b], sidx, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write_neut_species_conf(gkyl_gyrokinetic_multib_app* app, int sidx, double tm, int frame)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    gkyl_gyrokinetic_app_write_neut_species_conf(app->singleb_apps[b], sidx, tm, frame);
  }
}

//
// ............. Functions that group several species outputs ............... //
// 
void
gkyl_gyrokinetic_multib_app_write_mom(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_source_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_lbo_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_bgk_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_species_rad_emissivity(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_mom(app, i, tm, frame);
    gkyl_gyrokinetic_multib_app_write_neut_species_source_mom(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_calc_integrated_mom(gkyl_gyrokinetic_multib_app* app, double tm)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_calc_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_source_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_source_bgk_integrated_diagnostics(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_rad_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_species_boundary_flux_integrated_mom(app, i, tm);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_calc_neut_species_integrated_mom(app, i, tm);
    gkyl_gyrokinetic_multib_app_calc_neut_species_source_integrated_mom(app, i, tm);
  }

  gkyl_gyrokinetic_multib_app_calc_eirene_integrated_diagnostics(app, tm);
}

void
gkyl_gyrokinetic_multib_app_write_integrated_mom(gkyl_gyrokinetic_multib_app *app)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_species_source_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_species_source_bgk_integrated_diagnostics(app, i);
    gkyl_gyrokinetic_multib_app_write_species_lte_max_corr_status(app, i);
    gkyl_gyrokinetic_multib_app_write_species_rad_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_species_boundary_flux_integrated_mom(app, i);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_neut_species_source_integrated_mom(app, i);
    gkyl_gyrokinetic_multib_app_write_neut_species_lte_max_corr_status(app, i);
  }

  gkyl_gyrokinetic_multib_app_write_eirene_integrated_diagnostics(app);
}

void
gkyl_gyrokinetic_multib_app_write_conf(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_field(app, tm, frame);

  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_conf(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_conf(app, i, tm, frame);
  }

  gkyl_gyrokinetic_multib_app_write_eirene(app, tm, frame);
}

void
gkyl_gyrokinetic_multib_app_write_phase(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  for (int i=0; i<app->num_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_species_phase(app, i, tm, frame);
  }

  for (int i=0; i<app->num_neut_species; ++i) {
    gkyl_gyrokinetic_multib_app_write_neut_species_phase(app, i, tm, frame);
  }
}

void
gkyl_gyrokinetic_multib_app_write(gkyl_gyrokinetic_multib_app* app, double tm, int frame)
{
  gkyl_gyrokinetic_multib_app_write_phase(app, tm, frame);

  gkyl_gyrokinetic_multib_app_write_conf(app, tm, frame);
}

//
// ............. End of write functions ............... //
// 

struct gkyl_update_status
gkyl_gyrokinetic_multib_update(gkyl_gyrokinetic_multib_app* app, double dt)
{
  app->stat.nup += 1;
  struct timespec wst = gkyl_wall_clock();

  struct gkyl_update_status status = gyrokinetic_multib_update_ssp_rk3(app, dt);
  app->tcurr += status.dt_actual;

  app->stat.time_loop_tm += gkyl_time_diff_now_sec(wst);

  // Check for any CUDA errors during time step
  if (app->use_gpu)
    checkCuda(cudaGetLastError());

  return status;
}

struct gkyl_gyrokinetic_stat
gkyl_gyrokinetic_multib_app_stat(gkyl_gyrokinetic_multib_app* app)
{
  for (int i=0; i<app->num_species; ++i) {
    app->stat.n_iter_corr[i] = 0;
    app->stat.num_corr[i] = 0;
  }
  for (int i=0; i<app->num_neut_species; ++i) {
    app->stat.neut_n_iter_corr[i] = 0;
    app->stat.neut_num_corr[i] = 0;
  }

  for (int b=0; b<app->num_local_blocks; ++b) {
    // Add time spent on various operations for each local block.
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    struct gkyl_gyrokinetic_stat sb_stat = gkyl_gyrokinetic_app_stat(sbapp);

    app->stat.field_phi_rhs_tm += sb_stat.field_phi_rhs_tm;

    app->stat.species_lte_tm += sb_stat.species_lte_tm;
    app->stat.species_coll_mom_tm += sb_stat.species_coll_mom_tm;
    app->stat.species_coll_tm += sb_stat.species_coll_tm;
    app->stat.species_rad_mom_tm += sb_stat.species_rad_mom_tm;
    app->stat.species_rad_tm += sb_stat.species_rad_tm;
    app->stat.species_react_mom_tm += sb_stat.species_react_mom_tm;
    app->stat.species_react_tm += sb_stat.species_react_tm;
    app->stat.species_omega_cfl_tm += sb_stat.species_omega_cfl_tm;

    app->stat.neut_species_lte_tm += sb_stat.neut_species_lte_tm;
    app->stat.neut_species_coll_mom_tm += sb_stat.neut_species_coll_mom_tm;
    app->stat.neut_species_coll_tm += sb_stat.neut_species_coll_tm;
    app->stat.neut_species_react_mom_tm += sb_stat.neut_species_react_mom_tm;
    app->stat.neut_species_react_tm += sb_stat.neut_species_react_tm;
    app->stat.neut_species_omega_cfl_tm += sb_stat.neut_species_omega_cfl_tm;

    app->stat.species_bc_tm += sb_stat.species_bc_tm;
    app->stat.neut_species_bc_tm += sb_stat.neut_species_bc_tm;

    app->stat.species_io_tm += sb_stat.species_io_tm;
    app->stat.species_diag_calc_tm += sb_stat.species_diag_calc_tm;
    app->stat.species_diag_io_tm += sb_stat.species_diag_io_tm;

    app->stat.neut_species_io_tm += sb_stat.neut_species_io_tm;
    app->stat.neut_species_diag_calc_tm += sb_stat.neut_species_diag_calc_tm;
    app->stat.neut_species_diag_io_tm += sb_stat.neut_species_diag_io_tm;

    app->stat.n_species_omega_cfl += sb_stat.n_species_omega_cfl;
    app->stat.n_mom += sb_stat.n_mom;    
    app->stat.n_diag += sb_stat.n_diag;
    app->stat.n_io += sb_stat.n_io;
    app->stat.n_diag_io += sb_stat.n_diag_io;

    app->stat.n_neut_species_omega_cfl += sb_stat.n_neut_species_omega_cfl;
    app->stat.n_neut_mom += sb_stat.n_neut_mom;  
    app->stat.n_neut_diag += sb_stat.n_neut_diag;
    app->stat.n_neut_io += sb_stat.n_neut_io;
    app->stat.n_neut_diag_io += sb_stat.n_neut_diag_io;

    for (int i=0; i<app->num_species; ++i) {
      app->stat.n_iter_corr[i] += sb_stat.n_iter_corr[i];
      app->stat.num_corr[i] += sb_stat.num_corr[i];
    }
    for (int i=0; i<app->num_neut_species; ++i) {
      app->stat.neut_n_iter_corr[i] += sb_stat.neut_n_iter_corr[i];
      app->stat.neut_num_corr[i] += sb_stat.neut_num_corr[i];
    }

    // Timers not yet computed in app directly.
    app->stat.time_rate_diags_tm += sb_stat.fdot_tm + sb_stat.phidot_tm;
    app->stat.pos_shift_tm += sb_stat.species_pos_shift_tm + sb_stat.neut_species_pos_shift_tm + sb_stat.pos_shift_quasineut_tm;
    app->stat.io_tm += sb_stat.species_io_tm + sb_stat.species_diag_calc_tm + sb_stat.species_diag_io_tm + sb_stat.neut_species_io_tm
      + sb_stat.neut_species_diag_calc_tm + sb_stat.neut_species_diag_io_tm + sb_stat.field_io_tm + sb_stat.field_diag_calc_tm
      + sb_stat.field_diag_io_tm + sb_stat.app_io_tm;
    
    // Additions of several timers.
    app->stat.fwd_euler_sum_tm += sb_stat.species_coll_mom_tm + sb_stat.species_react_mom_tm + sb_stat.neut_species_coll_mom_tm
      + sb_stat.neut_species_react_mom_tm + sb_stat.species_rad_mom_tm + sb_stat.species_gyroavg_tm + sb_stat.species_collisionless_tm
      + sb_stat.species_coll_tm + sb_stat.species_diffusion_tm + sb_stat.species_rad_tm + sb_stat.species_react_tm
      + sb_stat.species_bflux_calc_tm+sb_stat.species_bflux_moms_tm + sb_stat.species_omega_cfl_tm + sb_stat.species_src_tm
      + sb_stat.neut_species_collisionless_tm + sb_stat.neut_species_coll_tm + sb_stat.neut_species_react_tm
      + sb_stat.neut_species_omega_cfl_tm + sb_stat.neut_species_src_tm + sb_stat.dfdt_dt_reduce_tm;

    app->stat.time_rate_diags_sum_tm += sb_stat.fdot_tm + sb_stat.phidot_tm;

    app->stat.pos_shift_sum_tm += sb_stat.species_pos_shift_tm + sb_stat.neut_species_pos_shift_tm + sb_stat.pos_shift_quasineut_tm;

    app->stat.io_sum_tm += sb_stat.species_io_tm + sb_stat.species_diag_calc_tm + sb_stat.species_diag_io_tm + sb_stat.neut_species_io_tm
    + sb_stat.neut_species_diag_calc_tm + sb_stat.neut_species_diag_io_tm + sb_stat.field_io_tm + sb_stat.field_diag_calc_tm
    + sb_stat.field_diag_io_tm + sb_stat.app_io_tm;

  }

  app->stat.fwd_euler_sum_tm += app->stat.fwd_euler_step_f_tm + app->stat.dfdt_dt_reduce_tm;
  app->stat.bc_sum_tm = app->stat.species_bc_tm + app->stat.neut_species_bc_tm;
  app->stat.field_sum_tm = app->stat.field_phi_rhs_tm  + app->stat.field_phi_solve_tm;
  app->stat.time_stepper_sum_tm = app->stat.fwd_euler_tm  + app->stat.field_tm + app->stat.bc_tm
    + app->stat.time_rate_diags_tm + app->stat.pos_shift_tm + app->stat.time_stepper_arithmetic_tm;

  return app->stat;
}

static inline
double
ratio_to_percent(double num, double den, double alt)
{
  return den > 1e-12 ? 100.*num/den : alt;
}

void
gkyl_gyrokinetic_multib_app_print_timings(gkyl_gyrokinetic_multib_app* app, FILE *iostream)
{
  struct gkyl_gyrokinetic_stat *stat = &app->stat;

  double bflux_tm = stat->species_bflux_calc_tm+stat->species_bflux_moms_tm;

  gkyl_gyrokinetic_multib_app_cout(app, iostream, "Timing:\n");
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "  - Time loop:                         %.4e sec.\n", stat->time_loop_tm);
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Forward Euler:                   %.4e sec. / %4.2f %%.\n", stat->fwd_euler_tm, ratio_to_percent(stat->fwd_euler_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision moments (charged):   %.4e sec. / %4.2f %%.\n", stat->species_coll_mom_tm          , ratio_to_percent(stat->species_coll_mom_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction moments (charged):    %.4e sec. / %4.2f %%.\n", stat->species_react_mom_tm         , ratio_to_percent(stat->species_react_mom_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision moments (neutral):   %.4e sec. / %4.2f %%.\n", stat->neut_species_coll_mom_tm     , ratio_to_percent(stat->neut_species_coll_mom_tm     ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction moments (neutral):    %.4e sec. / %4.2f %%.\n", stat->neut_species_react_mom_tm    , ratio_to_percent(stat->neut_species_react_mom_tm    ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Radiation moments:             %.4e sec. / %4.2f %%.\n", stat->species_rad_mom_tm           , ratio_to_percent(stat->species_rad_mom_tm           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species gyroaverage:           %.4e sec. / %4.2f %%.\n", stat->species_gyroavg_tm           , ratio_to_percent(stat->species_gyroavg_tm           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collisionless terms (charged): %.4e sec. / %4.2f %%.\n", stat->species_collisionless_tm     , ratio_to_percent(stat->species_collisionless_tm     ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision terms (charged):     %.4e sec. / %4.2f %%.\n", stat->species_coll_tm              , ratio_to_percent(stat->species_coll_tm              ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Diffusion (charged):           %.4e sec. / %4.2f %%.\n", stat->species_diffusion_tm         , ratio_to_percent(stat->species_diffusion_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Radiation terms:               %.4e sec. / %4.2f %%.\n", stat->species_rad_tm               , ratio_to_percent(stat->species_rad_tm               ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction terms (charged):      %.4e sec. / %4.2f %%.\n", stat->species_react_tm             , ratio_to_percent(stat->species_react_tm             ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Boundary fluxes (charged):     %.4e sec. / %4.2f %%.\n", bflux_tm                           , ratio_to_percent(bflux_tm                           ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ omega_cfl (charged):           %.4e sec. / %4.2f %%.\n", stat->species_omega_cfl_tm         , ratio_to_percent(stat->species_omega_cfl_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Sources (charged):             %.4e sec. / %4.2f %%.\n", stat->species_src_tm               , ratio_to_percent(stat->species_src_tm               ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collisionless terms (neutral): %.4e sec. / %4.2f %%.\n", stat->neut_species_collisionless_tm, ratio_to_percent(stat->neut_species_collisionless_tm,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Collision terms (neutral):     %.4e sec. / %4.2f %%.\n", stat->neut_species_coll_tm         , ratio_to_percent(stat->neut_species_coll_tm         ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Reaction terms (neutral):      %.4e sec. / %4.2f %%.\n", stat->neut_species_react_tm        , ratio_to_percent(stat->neut_species_react_tm        ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ omega_cfl (neutral):           %.4e sec. / %4.2f %%.\n", stat->neut_species_omega_cfl_tm    , ratio_to_percent(stat->neut_species_omega_cfl_tm    ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Sources (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_src_tm          , ratio_to_percent(stat->neut_species_src_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Time step reduction:           %.4e sec. / %4.2f %%.\n", stat->dfdt_dt_reduce_tm            , ratio_to_percent(stat->dfdt_dt_reduce_tm            ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Step f:                        %.4e sec. / %4.2f %%.\n", stat->fwd_euler_step_f_tm          , ratio_to_percent(stat->fwd_euler_step_f_tm          ,stat->fwd_euler_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->fwd_euler_sum_tm, stat->fwd_euler_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field solves:                    %.4e sec. / %4.2f %%.\n", stat->field_tm, ratio_to_percent(stat->field_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Phi eqn RHS:                   %.4e sec. / %4.2f %%.\n", stat->field_phi_rhs_tm  , ratio_to_percent(stat->field_phi_rhs_tm  ,stat->field_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Phi eqn solve:                 %.4e sec. / %4.2f %%.\n", stat->field_phi_solve_tm, ratio_to_percent(stat->field_phi_solve_tm,stat->field_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->field_sum_tm, stat->field_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Boundary conditions::            %.4e sec. / %4.2f %%.\n", stat->bc_tm, ratio_to_percent(stat->bc_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (charged):             %.4e sec. / %4.2f %%.\n", stat->species_bc_tm     , ratio_to_percent(stat->species_bc_tm     ,stat->bc_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_bc_tm, ratio_to_percent(stat->neut_species_bc_tm,stat->bc_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->bc_sum_tm, stat->bc_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Time rate diagnostics:           %.4e sec. / %4.2f %%.\n", stat->time_rate_diags_tm, ratio_to_percent(stat->time_rate_diags_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Charged species:               %.4e sec. / %4.2f %%.\n", stat->fdot_tm  , ratio_to_percent(stat->fdot_tm  ,stat->time_rate_diags_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Phi:                           %.4e sec. / %4.2f %%.\n", stat->phidot_tm, ratio_to_percent(stat->phidot_tm,stat->time_rate_diags_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->time_rate_diags_sum_tm, stat->time_rate_diags_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Positivity shift:                %.4e sec. / %4.2f %%.\n", stat->pos_shift_tm, ratio_to_percent(stat->pos_shift_tm,stat->time_loop_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (charged):             %.4e sec. / %4.2f %%.\n", stat->species_pos_shift_tm     , ratio_to_percent(stat->species_pos_shift_tm     ,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Species (neutral):             %.4e sec. / %4.2f %%.\n", stat->neut_species_pos_shift_tm, ratio_to_percent(stat->neut_species_pos_shift_tm,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Quasineutrality:               %.4e sec. / %4.2f %%.\n", stat->pos_shift_quasineut_tm   , ratio_to_percent(stat->pos_shift_quasineut_tm   ,stat->pos_shift_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "      ^ Accounted for:                 %4.2f %%.\n", ratio_to_percent(stat->pos_shift_sum_tm, stat->pos_shift_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Time stepper arithmetic:         %.4e sec. / %4.2f %%.\n", stat->time_stepper_arithmetic_tm, ratio_to_percent(stat->time_stepper_arithmetic_tm,stat->time_loop_tm, 0.0));

  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Accounted for:                   %4.2f %%.\n", ratio_to_percent(stat->time_stepper_sum_tm, stat->time_loop_tm, 100.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "  - I/O:                               %.4e sec.\n", stat->io_tm);
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * f write (charged):               %.4e sec. / %4.2f %%.\n", stat->species_io_tm            , ratio_to_percent(stat->species_io_tm            , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag calc (charged):     %.4e sec. / %4.2f %%.\n", stat->species_diag_calc_tm     , ratio_to_percent(stat->species_diag_calc_tm     , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag write (charged):    %.4e sec. / %4.2f %%.\n", stat->species_diag_io_tm       , ratio_to_percent(stat->species_diag_io_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * f write (neutral):               %.4e sec. / %4.2f %%.\n", stat->neut_species_io_tm       , ratio_to_percent(stat->neut_species_io_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag calc (neutral):     %.4e sec. / %4.2f %%.\n", stat->neut_species_diag_calc_tm, ratio_to_percent(stat->neut_species_diag_calc_tm, stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Species diag write (neutral):    %.4e sec. / %4.2f %%.\n", stat->neut_species_diag_io_tm  , ratio_to_percent(stat->neut_species_diag_io_tm  , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field write:                     %.4e sec. / %4.2f %%.\n", stat->field_io_tm              , ratio_to_percent(stat->field_io_tm              , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field diag calc:                 %.4e sec. / %4.2f %%.\n", stat->field_diag_calc_tm       , ratio_to_percent(stat->field_diag_calc_tm       , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Field diag write:                %.4e sec. / %4.2f %%.\n", stat->field_diag_io_tm         , ratio_to_percent(stat->field_diag_io_tm         , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Common write:                    %.4e sec. / %4.2f %%.\n", stat->app_io_tm                , ratio_to_percent(stat->app_io_tm                , stat->io_tm, 0.0));
  gkyl_gyrokinetic_multib_app_cout(app, iostream, "    * Accounted for:                   %4.2f %%.\n", ratio_to_percent(stat->io_sum_tm, stat->io_tm, 100.0));
}

void
gkyl_gyrokinetic_multib_app_stat_write(gkyl_gyrokinetic_multib_app* app)
{
  for (int b=0; b<app->num_local_blocks; ++b) {
    struct gkyl_gyrokinetic_app *sbapp = app->singleb_apps[b];
    gkyl_gyrokinetic_app_stat_write(sbapp);
  }
}

void
gkyl_gyrokinetic_multib_app_write_dt(gkyl_gyrokinetic_multib_app* app)
{
  int rank;
  gkyl_comm_get_rank(app->comm, &rank);

  if (rank == 0) {
    // Write integrated diagnostic moments.
    const char *fmt = "%s-%s.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, "dt");
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, "dt");

    struct timespec wtm = gkyl_wall_clock();
    if (app->is_first_dt_write_call) {
      struct gkyl_msgpack_map_elem io_meta_phi[] = {
        { .key = "Description", .elem_type = GKYL_MP_STRING, .cval = "Time step size." }
      };
      int io_meta_len[] = {app->io_meta_basic_len, 1};
      const struct gkyl_msgpack_map_elem* io_meta[] = {app->io_meta_basic, io_meta_phi};
      struct gkyl_msgpack_data *mt = gkyl_msgpack_create_union(sizeof(io_meta_len)/sizeof(int), io_meta_len, io_meta);

      gkyl_dynvec_write_wmeta(app->dts, fileNm, mt);
      app->is_first_dt_write_call = false;
      gkyl_msgpack_data_release(mt);
    }
    else {
      gkyl_dynvec_awrite(app->dts, fileNm);
    }
    app->stat.app_io_tm += gkyl_time_diff_now_sec(wtm);
    app->stat.n_diag_io += 1;
  }
  gkyl_dynvec_clear(app->dts);
}

void
gkyl_gyrokinetic_multib_app_save_dt(gkyl_gyrokinetic_multib_app* app, double tm, double dt)
{
  gkyl_dynvec_append(app->dts, tm, &dt);
}

void
gkyl_gyrokinetic_multib_app_release_geom(gkyl_gyrokinetic_multib_app* mbapp)
{
  if (mbapp->singleb_apps) {
    for (int i=0; i<mbapp->num_local_blocks; ++i)
      gkyl_gyrokinetic_app_release_geom(mbapp->singleb_apps[i]);
    gkyl_free(mbapp->singleb_apps);
  }  

  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(mbapp->block_comms[i]);
  gkyl_free(mbapp->block_comms);

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(mbapp->decomp[i]);
  gkyl_free(mbapp->decomp);

  gkyl_free(mbapp->local_blocks);    

  gkyl_rrobin_decomp_release(mbapp->round_robin);
  
  gkyl_gk_block_geom_release(mbapp->gk_block_geom);
  gkyl_block_topo_release(mbapp->block_topo);
  
  gkyl_comm_release(mbapp->comm);

  gkyl_msgpack_map_elem_release(mbapp->io_meta_basic_len, mbapp->io_meta_basic);

  gkyl_free(mbapp);
}

void
gkyl_gyrokinetic_multib_app_release(gkyl_gyrokinetic_multib_app* mbapp)
{
  for (int i=0; i<mbapp->num_neut_species; ++i) {
    for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_neut[i].send[bI]);
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_neut[i].recv[bI]);
    }
    gkyl_free(mbapp->mbcc_sync_neut[i].send);
    gkyl_free(mbapp->mbcc_sync_neut[i].recv);
  }
  gkyl_free(mbapp->mbcc_sync_neut);

  for (int i=0; i<mbapp->num_species; ++i) {
    for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_charged[i].send[bI]);
      gkyl_multib_comm_conn_release(mbapp->mbcc_sync_charged[i].recv[bI]);
    }
    gkyl_free(mbapp->mbcc_sync_charged[i].send);
    gkyl_free(mbapp->mbcc_sync_charged[i].recv);
  }
  gkyl_free(mbapp->mbcc_sync_charged);

  for (int bI=0; bI<mbapp->num_local_blocks; ++bI) {
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->send[bI]);
    gkyl_multib_comm_conn_release(mbapp->mbcc_sync_conf->recv[bI]);
  }
  gkyl_free(mbapp->mbcc_sync_conf->send);
  gkyl_free(mbapp->mbcc_sync_conf->recv);
  gkyl_free(mbapp->mbcc_sync_conf);

  if (mbapp->singleb_apps) {
    for (int i=0; i<mbapp->num_local_blocks; ++i)
      gkyl_gyrokinetic_app_release(mbapp->singleb_apps[i]);
    gkyl_free(mbapp->singleb_apps);
  }  

  gk_multib_field_release(mbapp->field);

  int num_blocks = gkyl_gk_block_geom_num_blocks(mbapp->gk_block_geom);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(mbapp->block_comms[i]);
  gkyl_free(mbapp->block_comms);

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(mbapp->decomp[i]);
  gkyl_free(mbapp->decomp);

  gkyl_free(mbapp->local_blocks);    

  gkyl_rrobin_decomp_release(mbapp->round_robin);
  
  gkyl_gk_block_geom_release(mbapp->gk_block_geom);
  gkyl_block_topo_release(mbapp->block_topo);
  
  gkyl_comm_release(mbapp->comm);

  gkyl_dynvec_release(mbapp->dts);

  gkyl_free(mbapp);
}
