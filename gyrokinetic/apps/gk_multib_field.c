#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>
#include <gkyl_multib_conn.h>
#include <gkyl_rrobin_decomp.h>
#include <assert.h>
#include <float.h>
#include <time.h>

// compute total number of ranges specified by cuts
static inline int
calc_cuts(int ndim, const int *cuts)
{
  int tc = 1;
  for (int d=0; d<ndim; ++d) tc *= cuts[d];
  return tc;
}

static struct gkyl_array**
gk_multib_field_mkarr(bool on_gpu, long nc, struct gkyl_range **ranges, int num_arr)
{
  // Allocate an array of double arrays (filled with zeros), each array with the
  // same number of components, but possibly a different size given by the volume
  // of a range.
  struct gkyl_array** arr = gkyl_malloc(num_arr * sizeof(struct gkyl_array*));;
  for (int i=0; i<num_arr; i++) {
    if (on_gpu)
      arr[i] = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, ranges[i]->volume);
    else
      arr[i] = gkyl_array_new(GKYL_DOUBLE, nc, ranges[i]->volume);
  }
  return arr;
}

static int**
gk_multib_field_new_connected_list(struct gkyl_gyrokinetic_multib_app *mbapp, int dir, int *nconnected)
{
  // Obtain a list of blocks connected along the specified direction.
  int num_blocks = mbapp->block_topo->num_blocks;
  int **block_list = gkyl_malloc(num_blocks*sizeof(int*));
  for (int bidx=0; bidx<num_blocks; ++bidx) {
    nconnected[bidx] = gkyl_multib_conn_get_num_connected(mbapp->block_topo, bidx, dir, 0, GKYL_CONN_ALL);
    block_list[bidx] = gkyl_malloc(nconnected[bidx]*sizeof(int));
    gkyl_multib_conn_get_connection(mbapp->block_topo, bidx, dir, 0, GKYL_CONN_ALL, block_list[bidx]);
  }

  return block_list;
}

static void
gk_multib_field_release_connected_list(struct gkyl_gyrokinetic_multib_app *mbapp, int **block_list)
{
  // Release the list of connected blocks.
  int num_blocks = mbapp->block_topo->num_blocks;
  for (int bidx=0; bidx<num_blocks; ++bidx)
    gkyl_free(block_list[bidx]);
  gkyl_free(block_list);
}

static void 
gk_multib_field_new_allgather_ranges(struct gk_multib_field *mbf, struct gkyl_gyrokinetic_multib_app *mbapp, int dir,
  struct gkyl_range **multibz_ranges, struct gkyl_range **multibz_ranges_ext)
{
  // Construct the local and global ranges for the allgather along a given
  // direction 'dir'. This function allocates 'multib_ranges' and
  // 'multib_ranges_ext', which must be freed when releasing mbf.

  // Get blocks connected along the specified direction.
  int nconnected[mbapp->block_topo->num_blocks];
  int **block_list = gk_multib_field_new_connected_list(mbapp, dir, nconnected);

  // Construct the local and global ranges for the allgather
  int nghost[] = {1, 1 ,1};
  int *local_blocks = mbapp->local_blocks;
  for (int bI= 0; bI<mbf->num_local_blocks; bI++) {
    int bid = local_blocks[bI];
    multibz_ranges[bI] = gkyl_malloc(sizeof(struct gkyl_range));
    multibz_ranges_ext[bI] = gkyl_malloc(sizeof(struct gkyl_range));
    gkyl_multib_comm_conn_create_multib_ranges_in_dir(multibz_ranges_ext[bI],
      multibz_ranges[bI], nghost, nconnected[bid], block_list[bid], dir, mbapp->decomp);
  }

  gk_multib_field_release_connected_list(mbapp, block_list);
}

static void
gk_multib_field_new_allgather_comm_conns(struct gk_multib_field *mbf, struct gkyl_gyrokinetic_multib_app *mbapp, int dir,
  struct gkyl_range **multib_ranges_ext, struct gkyl_multib_comm_conn **mbcc_allgather_send,
  struct gkyl_multib_comm_conn **mbcc_allgather_recv)
{
  // Construct the comm_conns for the allgather in a given direction.

  int my_rank, num_ranks;
  gkyl_comm_get_rank(mbapp->comm, &my_rank);
  gkyl_comm_get_size(mbapp->comm, &num_ranks);

  // Get branks (number of cuts per block)
  int num_blocks = mbapp->block_topo->num_blocks;
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, i);
    branks[i] = calc_cuts(mbf->cdim, bgi->cuts);
  }
  
  // Get blocks connected along the specified direction.
  int nconnected[num_blocks];
  int **block_list = gk_multib_field_new_connected_list(mbapp, dir, nconnected);

  // Construct the comm_conns for the allgather
  int nghost[] = {1, 1 ,1};
  int rank_list[num_ranks];
  for (int bI= 0; bI<mbf->num_local_blocks; bI++) {
    int bid = mbapp->local_blocks[bI];
    gkyl_rrobin_decomp_getranks(mbapp->round_robin, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbcc_allgather_send[bI] = gkyl_multib_comm_conn_new_send_from_connections(bid, brank, 
      nghost, nconnected[bid], block_list[bid], dir, mbapp->decomp);
    mbcc_allgather_recv[bI] = gkyl_multib_comm_conn_new_recv_from_connections(bid, brank,
      nghost, nconnected[bid], block_list[bid], dir, mbapp->decomp);

    for (int ns=0; ns<mbcc_allgather_send[bI]->num_comm_conn; ++ns) {
      // Need to get the actual rank that owns this cut.
      int rank_idx = mbcc_allgather_send[bI]->comm_conn[ns].rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, 
        mbcc_allgather_send[bI]->comm_conn[ns].block_id, rank_list);
      mbcc_allgather_send[bI]->comm_conn[ns].rank = rank_list[rank_idx];
      // Make range the local range (a subrange of local_ext).
      mbcc_allgather_send[bI]->comm_conn[ns].range = mbapp->singleb_apps[bI]->local;
    }
    for (int nr=0; nr<mbcc_allgather_recv[bI]->num_comm_conn; ++nr) {
      // Need to get the actual rank that owns this cut.
      int rank_idx = mbcc_allgather_recv[bI]->comm_conn[nr].rank;
      gkyl_rrobin_decomp_getranks(mbapp->round_robin, 
        mbcc_allgather_recv[bI]->comm_conn[nr].block_id, rank_list);
      mbcc_allgather_recv[bI]->comm_conn[nr].rank = rank_list[rank_idx];
      // Make range a subrange.
      gkyl_sub_range_init(&mbcc_allgather_recv[bI]->comm_conn[nr].range,
        multib_ranges_ext[bI], mbcc_allgather_recv[bI]->comm_conn[nr].range.lower,
        mbcc_allgather_recv[bI]->comm_conn[nr].range.upper);
    }

    // Sort connections according to rank and block ID (needed by NCCL).
    gkyl_multib_comm_conn_sort(mbcc_allgather_send[bI]);
    gkyl_multib_comm_conn_sort(mbcc_allgather_recv[bI]);
  }

  gk_multib_field_release_connected_list(mbapp, block_list);
  gkyl_free(branks);
}

static struct gkyl_range **
gk_multib_field_new_multib_to_global_ranges(struct gk_multib_field *mbf,
  struct gkyl_gyrokinetic_multib_app *mbapp, int dir, struct gkyl_range **multib_ranges)
{
  // Create ranges for copying smoothed quantity from multib to global after smoothing.

  // Get blocks connected along the specified direction.
  int nconnected[mbapp->block_topo->num_blocks];
  int **block_list = gk_multib_field_new_connected_list(mbapp, dir, nconnected);

  struct gkyl_range **parent_subranges = gkyl_malloc(mbf->num_local_blocks * sizeof(struct gkyl_range *));
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    int bid = mbapp->local_blocks[bI];
    int shift[GKYL_MAX_DIM] = {0};
    for (int i=0; i<nconnected[bid]; i++) {
      if (block_list[bid][i] == bid)
        break;
      else 
        shift[dir] += gkyl_range_shape(&mbapp->decomp[block_list[bid][i]]->parent_range, dir);
    }
    struct gkyl_range shifted_parent_range;
    gkyl_range_shift(&shifted_parent_range, &mbapp->singleb_apps[bI]->global, shift);
    parent_subranges[bI] = gkyl_malloc(sizeof(struct gkyl_range));
    int inter = gkyl_sub_range_intersect(parent_subranges[bI],
      multib_ranges[bI], &shifted_parent_range);
  }

  gk_multib_field_release_connected_list(mbapp, block_list);

  return parent_subranges;
}

static struct gkyl_range **
gk_multib_field_new_multib_to_local_ranges(struct gk_multib_field *mbf,
  struct gkyl_gyrokinetic_multib_app *mbapp, int dir, struct gkyl_range **multib_ranges)
{
  // Create ranges for copying smoothed quantity from multib to local after smoothing.

  // Get blocks connected along the specified direction.
  int nconnected[mbapp->block_topo->num_blocks];
  int **block_list = gk_multib_field_new_connected_list(mbapp, dir, nconnected);

  struct gkyl_range **block_subranges = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_range *));
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    int bid = mbapp->local_blocks[bI];
    int shift[GKYL_MAX_DIM] = {0};
    for (int i=0; i<nconnected[bid]; i++) {
      if (block_list[bid][i] == bid)
        break;
      else 
        shift[dir] += gkyl_range_shape(&mbapp->decomp[block_list[bid][i]]->parent_range, dir);
    }
    struct gkyl_range shifted_block_range;
    gkyl_range_shift(&shifted_block_range, &mbapp->singleb_apps[bI]->local, shift);
    block_subranges[bI] = gkyl_malloc(sizeof(struct gkyl_range));
    int inter = gkyl_sub_range_intersect(block_subranges[bI],
      multib_ranges[bI], &shifted_block_range);
  }

  gk_multib_field_release_connected_list(mbapp, block_list);

  return block_subranges; 
}

static void
gk_multib_field_new_par_smooth(const struct gkyl_gyrokinetic_multib *mbinp,
  struct gkyl_gyrokinetic_multib_app *mbapp, struct gk_multib_field *mbf)
{
  // Initialize objects needed for the multiblock parallel smoothing.
  int dir = mbf->cdim-1;
 
  // Construct the local and global ranges for the allgather along z.
  mbf->multibz_ranges = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_range *));
  mbf->multibz_ranges_ext = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_range *));
  gk_multib_field_new_allgather_ranges(mbf, mbapp, dir, mbf->multibz_ranges, mbf->multibz_ranges_ext);

  // Allocate global-in-z arrays for charge density and potential.
  int num_basis = mbapp->singleb_apps[0]->basis.num_basis;
  mbf->phi_multibz_dg = gk_multib_field_mkarr(mbapp->use_gpu, num_basis, mbf->multibz_ranges_ext, mbf->num_local_blocks);
  mbf->phi_multibz_smooth = gk_multib_field_mkarr(mbapp->use_gpu, num_basis, mbf->multibz_ranges_ext, mbf->num_local_blocks);
  mbf->rho_c_multibz_dg = gk_multib_field_mkarr(mbapp->use_gpu, num_basis, mbf->multibz_ranges_ext, mbf->num_local_blocks);
  mbf->rho_c_multibz_smooth = gk_multib_field_mkarr(mbapp->use_gpu, num_basis, mbf->multibz_ranges_ext, mbf->num_local_blocks);

  // Construct the comm_conns for the allgather along z.
  mbf->mbcc_allgatherz_send = gkyl_malloc(mbf->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  mbf->mbcc_allgatherz_recv = gkyl_malloc(mbf->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  gk_multib_field_new_allgather_comm_conns(mbf, mbapp, dir, mbf->multibz_ranges_ext,
    mbf->mbcc_allgatherz_send, mbf->mbcc_allgatherz_recv);
  
  // Create ranges for copying smoothed quantity from multib to global after smoothing.
  mbf->parent_subrangesz = gk_multib_field_new_multib_to_global_ranges(mbf, mbapp, dir, mbf->multibz_ranges);

  // Create ranges for copying smoothed quantity from multib to local after smoothing.
  mbf->block_subrangesz = gk_multib_field_new_multib_to_local_ranges(mbf, mbapp, dir, mbf->multibz_ranges);

  // Allocate the weights used for the parallel smoother.
  struct gkyl_array **lhs_weight_local = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  struct gkyl_array **rhs_weight_local = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    rhs_weight_local[bI] = mkarr(mbapp->use_gpu, sbapp->basis.num_basis, sbapp->local_ext.volume);
    if (mbf->cdim == 1) {
      lhs_weight_local[bI] = gkyl_array_acquire(sbapp->field->epsilon);
    }
    else {
      lhs_weight_local[bI] = mkarr(mbapp->use_gpu, sbapp->basis.num_basis, sbapp->local_ext.volume);
    }
  }

  // Set the smoothing local RHS weight to 1, and the LHS weight to 1 if cdim>1,
  // and gather them along the magnetic field.
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    gkyl_array_shiftc(rhs_weight_local[bI], sqrt(pow(2,mbf->cdim)), 0);
    if (mbf->cdim > 1)
      gkyl_array_shiftc(lhs_weight_local[bI], sqrt(pow(2,mbf->cdim)), 0);
  }
  mbf->lhs_weight_multibz = gk_multib_field_mkarr(mbapp->use_gpu, num_basis, mbf->multibz_ranges_ext, mbf->num_local_blocks);
  mbf->rhs_weight_multibz = gk_multib_field_mkarr(mbapp->use_gpu, num_basis, mbf->multibz_ranges_ext, mbf->num_local_blocks);
  int stat;
  stat = gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbf->num_local_blocks, 
    mbapp->local_blocks, mbf->mbcc_allgatherz_send, mbf->mbcc_allgatherz_recv, lhs_weight_local, 
    mbf->lhs_weight_multibz);
  stat = gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbf->num_local_blocks, 
    mbapp->local_blocks, mbf->mbcc_allgatherz_send, mbf->mbcc_allgatherz_recv, rhs_weight_local, 
    mbf->rhs_weight_multibz);

  for (int bI= 0; bI<mbf->num_local_blocks; bI++) {
    gkyl_array_release(lhs_weight_local[bI]);
    gkyl_array_release(rhs_weight_local[bI]);
  }
  gkyl_free(lhs_weight_local);
  gkyl_free(rhs_weight_local);

  // Create the parallel smoother.
  mbf->fem_parproj = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_fem_parproj*));
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    int bid = mbapp->local_blocks[bI];
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];

    // Choose no BC for the parallel smoother, unless we are in the core in 2x
    // in which case periodic BCs are needed.
    enum gkyl_fem_parproj_bc_type fem_parbc = GKYL_FEM_PARPROJ_NONE;
    const struct gkyl_gk_block_geom_info *bgi = gkyl_gk_block_geom_get_block(mbapp->gk_block_geom, bid);
    enum gkyl_tok_geo_type ftype = bgi->geometry.tok_grid_info.ftype;
    if (mbf->cdim == 2 && (ftype == GKYL_GEOMETRY_TOKAMAK_CORE || ftype == GKYL_GEOMETRY_TOKAMAK_CORE_R || ftype == GKYL_GEOMETRY_TOKAMAK_CORE_L) && !mbf->info.half_domain)
      fem_parbc = GKYL_FEM_PARPROJ_PERIODIC;

    mbf->fem_parproj[bI] = gkyl_fem_parproj_new(mbf->multibz_ranges[bI], &sbapp->grid,
      &sbapp->basis, fem_parbc, 0, mbf->lhs_weight_multibz[bI], mbf->rhs_weight_multibz[bI], mbapp->use_gpu);
  }
}

static bool
in_array_int(int inp, const int *arr, int num_elements)
{
  // Check if 'inp' is in the array 'arr' which has 'num_elements'.
  bool found = false;
  for (int i=0; i<num_elements; i++) {
    if (arr[i] == inp) {
      found = true;
      break;
    }
  }
  return found;
}

static bool
gk_multib_is_bid_connected_in_dir(int bidx, struct gkyl_gyrokinetic_multib_app *mbapp, int my_bidx, int dir)
{
  // Get blocks connected along the specified direction.
  int nconnected[mbapp->block_topo->num_blocks];
  int **block_list = gk_multib_field_new_connected_list(mbapp, dir, nconnected);

  bool is_in_conn_dir = in_array_int(bidx, block_list[my_bidx], nconnected[my_bidx]);

  gk_multib_field_release_connected_list(mbapp, block_list);
  return is_in_conn_dir;
}

static void
gk_multib_field_new_perp_solve(const struct gkyl_gyrokinetic_multib *mbinp,
  struct gkyl_gyrokinetic_multib_app *mbapp, struct gk_multib_field *mbf)
{
  // Initialize objects needed for the multiblock perpendicular Poisson solve.
  int dir = 0; // Note that for cdim=3 we here assume there are is everywhere a
               // single block along y.
 
  // Construct the local and global ranges for the perpendicular allgather.
  mbf->multib_perp_ranges = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_range *));
  mbf->multib_perp_ranges_ext = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_range *));
  gk_multib_field_new_allgather_ranges(mbf, mbapp, dir, mbf->multib_perp_ranges, mbf->multib_perp_ranges_ext);

  // Allocate global-in-x arrays for charge density and potential.
  int num_basis = mbapp->singleb_apps[0]->basis.num_basis;
  mbf->phi_multib_perp = gk_multib_field_mkarr(mbapp->use_gpu, num_basis, mbf->multib_perp_ranges_ext, mbf->num_local_blocks);
  mbf->rho_c_multib_perp = gk_multib_field_mkarr(mbapp->use_gpu, num_basis, mbf->multib_perp_ranges_ext, mbf->num_local_blocks);

  // Construct the comm_conns for the allgather along z.
  mbf->mbcc_allgather_perp_send = gkyl_malloc(mbf->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  mbf->mbcc_allgather_perp_recv = gkyl_malloc(mbf->num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  gk_multib_field_new_allgather_comm_conns(mbf, mbapp, dir, mbf->multib_perp_ranges_ext,
    mbf->mbcc_allgather_perp_send, mbf->mbcc_allgather_perp_recv);
  
  // Create ranges for copying smoothed quantity from multib to global.
  mbf->parent_subranges_perp = gk_multib_field_new_multib_to_global_ranges(mbf, mbapp, dir, mbf->multib_perp_ranges);

  // Create ranges for copying smoothed quantity from multib to local.
  mbf->block_subranges_perp = gk_multib_field_new_multib_to_local_ranges(mbf, mbapp, dir, mbf->multib_perp_ranges);

  // Allocate and gather the polarization weights.
  struct gkyl_array **epsilon_local = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    epsilon_local[bI] = gkyl_array_acquire(sbapp->field->epsilon);
  }
  mbf->epsilon_multib_perp = gk_multib_field_mkarr(mbapp->use_gpu, mbapp->singleb_apps[0]->field->epsilon->ncomp,
    mbf->multib_perp_ranges_ext, mbf->num_local_blocks);
  int stat = gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbf->num_local_blocks, 
    mbapp->local_blocks, mbf->mbcc_allgather_perp_send, mbf->mbcc_allgather_perp_recv, epsilon_local, 
    mbf->epsilon_multib_perp);
  for (int bI= 0; bI<mbf->num_local_blocks; bI++)
    gkyl_array_release(epsilon_local[bI]);
  gkyl_free(epsilon_local);

  // Create the perpendicular solve.
  mbf->fem_poisson = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_fem_poisson_perp*));
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    int bid = mbapp->local_blocks[bI];
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];

    // Obtain the BCs for connected blocks.
    // Loop through input BCs and find those for the current connnected blocks.
    struct gkyl_poisson_bc bcs = { };
    for (int i=0; i<mbinp->field.num_physical_bcs; i++) {
      const struct gkyl_gyrokinetic_bc *bc_curr = &mbinp->field.bcs[i];
      if (gk_multib_is_bid_connected_in_dir(bc_curr->bidx, mbapp, bid, dir)) {
        if (bc_curr->edge == GKYL_LOWER_EDGE) {
          bcs.lo_type[bc_curr->dir] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_curr->type);
          for (int k=0; k<3; k++)
            bcs.lo_value[bc_curr->dir].v[k] = bc_curr->value[k];
        }
        else {
          bcs.up_type[bc_curr->dir] = gkyl_gyrokinetic_translate_poisson_bc_type(bc_curr->type);
          for (int k=0; k<3; k++)
            bcs.up_value[bc_curr->dir].v[k] = bc_curr->value[k];
        }
      }
    }

    mbf->fem_poisson[bI] = gkyl_fem_poisson_perp_new(mbf->multib_perp_ranges[bI], &sbapp->grid, sbapp->basis,
      &bcs, mbf->info.bias_line_list, mbf->epsilon_multib_perp[bI], NULL, mbapp->use_gpu);
  }

}

// Compute the electrostatic potential.
static void
gk_multib_field_1x_rhs(gkyl_gyrokinetic_multib_app *mbapp, struct gk_multib_field *mbf,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  // Every local block calculates its charge density.
  for (int bI=0; bI<mbf->num_local_blocks; bI++) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    // Construct fin for the local block.
    const struct gkyl_array *fin_local_block[mbapp->num_species];
    struct gkyl_array **bflux_local_block[mbapp->num_species];
    int lin_idx = bI * mbapp->num_species;
    for (int i=0; i<mbapp->num_species; ++i) {
      fin_local_block[i] = fin[lin_idx+i];
      bflux_local_block[i] = bflux[lin_idx+i];
    }
    // Accumulate rho_c in local block.
    gk_field_accumulate_rho_c(sbapp, sbapp->field, fin_local_block, bflux_local_block);
  }

  struct timespec wst = gkyl_wall_clock();

  // Gather the charge density along the magnetic field.
  int stat_par_rho = gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbf->num_local_blocks, mbapp->local_blocks,
    mbf->mbcc_allgatherz_send, mbf->mbcc_allgatherz_recv, mbf->rho_c_local, mbf->rho_c_multibz_dg);
  // Make charge density continuous on the multibz range.
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    gkyl_fem_parproj_set_rhs(mbf->fem_parproj[bI], mbf->rho_c_multibz_dg[bI], mbf->rho_c_multibz_dg[bI]);
    gkyl_fem_parproj_solve(mbf->fem_parproj[bI], mbf->rho_c_multibz_smooth[bI]);
  }

  // Copy continuous charge density back to apps.
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    gkyl_array_copy_range_to_range(mbapp->singleb_apps[bI]->field->phi_smooth, mbf->rho_c_multibz_smooth[bI],
      &mbapp->singleb_apps[bI]->local, mbf->block_subrangesz[bI]);
  }

  mbapp->stat.field_phi_solve_tm += gkyl_time_diff_now_sec(wst);
}

// Compute the electrostatic potential.
static void
gk_multib_field_2x_rhs(gkyl_gyrokinetic_multib_app *mbapp, struct gk_multib_field *mbf,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  // Every local block calculates its charge density.
  for (int bI=0; bI<mbf->num_local_blocks; bI++) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    // Construct fin for the local block.
    const struct gkyl_array *fin_local_block[mbapp->num_species];
    struct gkyl_array **bflux_local_block[mbapp->num_species];
    int lin_idx = bI * mbapp->num_species;
    for (int i=0; i<mbapp->num_species; ++i) {
      fin_local_block[i] = fin[lin_idx+i];
      bflux_local_block[i] = bflux[lin_idx+i];
    }
    // Accumulate rho_c in local block.
    gk_field_accumulate_rho_c(sbapp, sbapp->field, fin_local_block, bflux_local_block);
  }

  struct timespec wst = gkyl_wall_clock();

  // Gather the charge density along the magnetic field.
  int stat_par_rho = gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbf->num_local_blocks, mbapp->local_blocks,
    mbf->mbcc_allgatherz_send, mbf->mbcc_allgatherz_recv, mbf->rho_c_local, mbf->rho_c_multibz_dg);
  // Make charge density continuous on the multibz range.
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    gkyl_fem_parproj_set_rhs(mbf->fem_parproj[bI], mbf->rho_c_multibz_dg[bI], mbf->rho_c_multibz_dg[bI]);
    gkyl_fem_parproj_solve(mbf->fem_parproj[bI], mbf->rho_c_multibz_smooth[bI]);
  }

  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    // Copy continuous charge density back to apps.
    gkyl_array_copy_range_to_range(sbapp->field->rho_c_global_smooth,
      mbf->rho_c_multibz_smooth[bI], &sbapp->global, mbf->parent_subrangesz[bI]);
    // Copy from block-global to block local.
    gkyl_array_copy_range_to_range(sbapp->field->rho_c, sbapp->field->rho_c_global_smooth,
      &sbapp->local, &sbapp->field->global_sub_range);
  }

  //
  // Solve the perpendicular Poisson problem.
  //
  // Gather the charge density in the perpendicular direction.
  int stat_perp = gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbf->num_local_blocks, mbapp->local_blocks,
    mbf->mbcc_allgather_perp_send, mbf->mbcc_allgather_perp_recv, mbf->rho_c_local, mbf->rho_c_multib_perp);
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    // Solve the perp problem.
    gkyl_fem_poisson_perp_set_rhs(mbf->fem_poisson[bI], mbf->rho_c_multib_perp[bI]);
    gkyl_fem_poisson_perp_solve(mbf->fem_poisson[bI], mbf->phi_multib_perp[bI]);
    // Copy the potential from the mulib range to local.
    gkyl_array_copy_range_to_range(mbapp->singleb_apps[bI]->field->phi_smooth, mbf->phi_multib_perp[bI],
      &mbapp->singleb_apps[bI]->local, mbf->block_subranges_perp[bI]);
  }
  //
  // Finished solving the perpendicular Poisson problem.
  //

  // Gather the potential along the magnetic field.
  int stat_par_phi = gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbf->num_local_blocks, mbapp->local_blocks,
    mbf->mbcc_allgatherz_send, mbf->mbcc_allgatherz_recv, mbf->phi_local, mbf->phi_multibz_dg);
  // Make the potential continuous along B on the multibz range.
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    gkyl_fem_parproj_set_rhs(mbf->fem_parproj[bI], mbf->phi_multibz_dg[bI], mbf->phi_multibz_dg[bI]);
    gkyl_fem_parproj_solve(mbf->fem_parproj[bI], mbf->phi_multibz_smooth[bI]);
  }
  // Copy continuous potential back to apps.
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    // Copy from multib to block global.
    gkyl_array_copy_range_to_range(sbapp->field->rho_c_global_smooth,
      mbf->phi_multibz_smooth[bI], &sbapp->global, mbf->parent_subrangesz[bI]);
    // Copy from block-global to block local.
    gkyl_array_copy_range_to_range(sbapp->field->phi_smooth, sbapp->field->rho_c_global_smooth,
      &sbapp->local, &sbapp->field->global_sub_range);
  }

  mbapp->stat.field_phi_solve_tm += gkyl_time_diff_now_sec(wst);
}


void
gk_multib_field_twistshift_and_smooth(gkyl_gyrokinetic_multib_app *mbapp, struct gk_multib_field *mbf,
  struct gkyl_array **arr_local, struct gkyl_array **arr_global_dg, struct gkyl_array **arr_global_smooth)
{
  // There are three arrays:
  // 1) The local one (field->rho_c = rho_c_local[bI])
  // 2) The global non-smooth one (field->rho_c_global_dg) needs an mb equivalent
  // 3) The global smooth one
  // Now apply TS in place on rho for each block
  int par_dir = mbf->cdim-1; // Parallel direction index.
  for (int bI=0; bI<mbapp->num_local_blocks; bI++) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    struct gk_field *field = sbapp->field;
    gkyl_comm_array_allgather(sbapp->comm, &sbapp->local, &sbapp->global, 
      arr_local[bI], arr_global_dg[bI]);

    // Apply TS BC in the lower parallel boundary, and
    // fill upper parallel boundary ghost with skin boundary value.
    if (mbf->num_blocks_below[bI] > 0) {
      gkyl_bc_twistshift_advance(field->bc_ts_lo, arr_global_dg[bI], arr_global_dg[bI]);
    }
    else{
      gkyl_bc_basic_gyrokinetic_advance(field->gfss_bc_op_core_lo, field->bc_buffer, arr_global_dg[bI]);
    }
   
    // Apply TS BC in the upper parallel boundary, and
    // fill lower parallel boundary ghost with skin boundary value.
    if (mbf->num_blocks_above[bI] > 0) {
      gkyl_bc_twistshift_advance(field->bc_ts_up, arr_global_dg[bI], arr_global_dg[bI]);
    }
    else {
      gkyl_bc_basic_gyrokinetic_advance(field->gfss_bc_op_core_up, field->bc_buffer, arr_global_dg[bI]);
    }

    // Smooth the the DG array.
    gkyl_fem_parproj_set_rhs(field->fem_parproj_phi_core, arr_global_dg[bI], arr_global_dg[bI]);
    gkyl_fem_parproj_solve(field->fem_parproj_phi_core, arr_global_smooth[bI]);

    // Copy global, continuous FEM array to a local array.
    gkyl_array_copy_range_to_range(arr_local[bI], arr_global_smooth[bI], &sbapp->local, &field->global_sub_range);
  }
}

// Compute the electrostatic potential.
static void
gk_multib_field_3x_rhs(gkyl_gyrokinetic_multib_app *mbapp, struct gk_multib_field *mbf,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  // Every local block calculates its charge density.
  for (int bI=0; bI<mbf->num_local_blocks; bI++) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    // Construct fin for the local block.
    const struct gkyl_array *fin_local_block[mbapp->num_species];
    struct gkyl_array **bflux_local_block[mbapp->num_species];
    int lin_idx = bI * mbapp->num_species;
    for (int i=0; i<mbapp->num_species; ++i) {
      fin_local_block[i] = fin[lin_idx+i];
      bflux_local_block[i] = bflux[lin_idx+i];
    }
    // Accumulate rho_c in local block.
    gk_field_accumulate_rho_c(sbapp, sbapp->field, fin_local_block, bflux_local_block);
  }

  // Now sync the charge density
  struct gkyl_array *rhos[mbapp->num_local_blocks];
  for (int bI=0; bI<mbapp->num_local_blocks; bI++) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    rhos[bI] = sbapp->field->rho_c;
  }
  gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
    mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, rhos, rhos);

  // Now apply TS in place on rho for each block
  gk_multib_field_twistshift_and_smooth(mbapp, mbf, mbf->rho_c_local, mbf->rho_c_global_dg, mbf->rho_c_global_smooth);
  // Done applying TS now we need to solve the field

  //
  // Solve the perpendicular Poisson problem.
  //
  // Gather the charge density in the perpendicular direction.
  int stat_perp = gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbf->num_local_blocks, mbapp->local_blocks,
    mbf->mbcc_allgather_perp_send, mbf->mbcc_allgather_perp_recv, mbf->rho_c_local, mbf->rho_c_multib_perp);
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    // Solve the perp problem.
    gkyl_fem_poisson_perp_set_rhs(mbf->fem_poisson[bI], mbf->rho_c_multib_perp[bI]);
    gkyl_fem_poisson_perp_solve(mbf->fem_poisson[bI], mbf->phi_multib_perp[bI]);
    // Copy the potential from the mulib range to local.
    gkyl_array_copy_range_to_range(mbapp->singleb_apps[bI]->field->phi_smooth, mbf->phi_multib_perp[bI],
      &mbapp->singleb_apps[bI]->local, mbf->block_subranges_perp[bI]);
  }
  //
  // Finished solving the perpendicular Poisson problem.
  //

  // Now sync the charge density
  struct gkyl_array *phis[mbapp->num_local_blocks];
  for (int bI=0; bI<mbapp->num_local_blocks; bI++) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    phis[bI] = sbapp->field->phi_smooth;
  }
  gkyl_multib_comm_conn_array_transfer(mbapp->comm, mbapp->num_local_blocks, mbapp->local_blocks,
    mbapp->mbcc_sync_conf->send, mbapp->mbcc_sync_conf->recv, phis, phis);

  // Now apply TS in place on phi for each block
  gk_multib_field_twistshift_and_smooth(mbapp, mbf, mbf->phi_local, mbf->rho_c_global_dg, mbf->phi_global_smooth);
  // Done with TS and smooth of phi
}

// Initialize multib field object
struct gk_multib_field* 
gk_multib_field_new(const struct gkyl_gyrokinetic_multib *mbinp, struct gkyl_gyrokinetic_multib_app *mbapp)
{
  struct gk_multib_field *mbf = gkyl_malloc(sizeof(struct gk_multib_field));

  mbf->info = mbinp->field;
  mbf->gkfield_id = mbf->info.gkfield_id? mbf->info.gkfield_id : GKYL_GK_FIELD;
  mbf->num_local_blocks = mbapp->num_local_blocks;
  mbf->cdim = mbapp->block_topo->ndim;
  mbf->half_domain = mbf->info.half_domain ? mbf->info.half_domain : false;

  // Allocate local arrays for charge density and potential.
  mbf->phi_local = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  mbf->rho_c_local = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  mbf->rho_c_global_dg = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  mbf->phi_global_smooth = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  mbf->rho_c_global_smooth = gkyl_malloc(mbf->num_local_blocks* sizeof(struct gkyl_array*));
  for (int bI=0; bI<mbf->num_local_blocks; ++bI) {
    struct gkyl_gyrokinetic_app *sbapp = mbapp->singleb_apps[bI];
    mbf->phi_local[bI] = gkyl_array_acquire(sbapp->field->phi_smooth);
    mbf->rho_c_local[bI] = gkyl_array_acquire(sbapp->field->rho_c);
    mbf->rho_c_global_dg[bI] = gkyl_array_acquire(sbapp->field->rho_c_global_dg);
    mbf->phi_global_smooth[bI] = gkyl_array_acquire(sbapp->field->phi_fem);
    mbf->rho_c_global_smooth[bI] = gkyl_array_acquire(sbapp->field->rho_c_global_smooth);
  }

  // Initialize objects needed for the multiblock parallel smoothing.
  gk_multib_field_new_par_smooth(mbinp, mbapp, mbf);

  if (mbf->cdim > 1) {
    // Initialize objects needed for the multiblock perpendicular solve.
    gk_multib_field_new_perp_solve(mbinp, mbapp, mbf);
  }

  for (int bI = 0; bI < mbapp->num_local_blocks; bI++) {
    int dir = mbf->cdim-1;
    mbf->num_blocks_below[bI] = gkyl_multib_conn_get_num_connected(mbapp->block_topo,
      mbapp->local_blocks[bI], dir, 0, GKYL_CONN_BELOW);
    mbf->num_blocks_above[bI] = gkyl_multib_conn_get_num_connected(mbapp->block_topo,
      mbapp->local_blocks[bI], dir, 0, GKYL_CONN_ABOVE);
  }

  if (mbf->cdim == 1)
    mbf->rhs_func = gk_multib_field_1x_rhs;
  if (mbf->cdim == 2)
    mbf->rhs_func = gk_multib_field_2x_rhs;
  if (mbf->cdim == 3)
    mbf->rhs_func = gk_multib_field_3x_rhs;

  return mbf;
}

void
gk_multib_field_rhs(gkyl_gyrokinetic_multib_app *mbapp, struct gk_multib_field *mbf,
  const struct gkyl_array *fin[], struct gkyl_array **bflux[])
{
  // Compute the electrostatic potential.
  mbf->rhs_func(mbapp, mbf, fin, bflux);
}

// Release resources for multib field.
void
gk_multib_field_release(struct gk_multib_field *mbf)
{
  for (int bI= 0; bI<mbf->num_local_blocks; bI++) {
    gkyl_array_release(mbf->phi_local[bI]);
    gkyl_array_release(mbf->rho_c_local[bI]);
    gkyl_array_release(mbf->rho_c_global_dg[bI]);
    gkyl_array_release(mbf->phi_global_smooth[bI]);
    gkyl_array_release(mbf->rho_c_global_smooth[bI]);
  }
  gkyl_free(mbf->phi_local);
  gkyl_free(mbf->rho_c_local);
  gkyl_free(mbf->rho_c_global_dg);
  gkyl_free(mbf->phi_global_smooth);
  gkyl_free(mbf->rho_c_global_smooth);

  // Free memory allocated for parallel smoothing.
  for (int bI= 0; bI<mbf->num_local_blocks; bI++) {
    gkyl_free(mbf->multibz_ranges[bI]);
    gkyl_free(mbf->multibz_ranges_ext[bI]);
    gkyl_free(mbf->parent_subrangesz[bI]);
    gkyl_free(mbf->block_subrangesz[bI]);
    gkyl_multib_comm_conn_release(mbf->mbcc_allgatherz_send[bI]);
    gkyl_multib_comm_conn_release(mbf->mbcc_allgatherz_recv[bI]);
    gkyl_array_release(mbf->phi_multibz_dg[bI]);
    gkyl_array_release(mbf->phi_multibz_smooth[bI]);
    gkyl_array_release(mbf->rho_c_multibz_dg[bI]);
    gkyl_array_release(mbf->rho_c_multibz_smooth[bI]);
    gkyl_array_release(mbf->lhs_weight_multibz[bI]);
    gkyl_array_release(mbf->rhs_weight_multibz[bI]);
    gkyl_fem_parproj_release(mbf->fem_parproj[bI]);
  }
  gkyl_free(mbf->multibz_ranges);
  gkyl_free(mbf->multibz_ranges_ext);
  gkyl_free(mbf->parent_subrangesz);
  gkyl_free(mbf->block_subrangesz);
  gkyl_free(mbf->mbcc_allgatherz_send);
  gkyl_free(mbf->mbcc_allgatherz_recv);
  gkyl_free(mbf->phi_multibz_dg);
  gkyl_free(mbf->phi_multibz_smooth);
  gkyl_free(mbf->rho_c_multibz_dg);
  gkyl_free(mbf->rho_c_multibz_smooth);
  gkyl_free(mbf->lhs_weight_multibz);
  gkyl_free(mbf->rhs_weight_multibz);
  gkyl_free(mbf->fem_parproj);

  if (mbf->cdim > 1) {
    // Free memory allocated for perp solve.
    for (int bI= 0; bI<mbf->num_local_blocks; bI++) {
      gkyl_free(mbf->multib_perp_ranges[bI]);
      gkyl_free(mbf->multib_perp_ranges_ext[bI]);
      gkyl_free(mbf->parent_subranges_perp[bI]);
      gkyl_free(mbf->block_subranges_perp[bI]);
      gkyl_multib_comm_conn_release(mbf->mbcc_allgather_perp_send[bI]);
      gkyl_multib_comm_conn_release(mbf->mbcc_allgather_perp_recv[bI]);
      gkyl_array_release(mbf->phi_multib_perp[bI]);
      gkyl_array_release(mbf->rho_c_multib_perp[bI]);
      gkyl_array_release(mbf->epsilon_multib_perp[bI]);
      gkyl_fem_poisson_perp_release(mbf->fem_poisson[bI]);
    }
    gkyl_free(mbf->multib_perp_ranges);
    gkyl_free(mbf->multib_perp_ranges_ext);
    gkyl_free(mbf->parent_subranges_perp);
    gkyl_free(mbf->block_subranges_perp);
    gkyl_free(mbf->mbcc_allgather_perp_send);
    gkyl_free(mbf->mbcc_allgather_perp_recv);
    gkyl_free(mbf->phi_multib_perp);
    gkyl_free(mbf->rho_c_multib_perp);
    gkyl_free(mbf->epsilon_multib_perp);
    gkyl_free(mbf->fem_poisson);
  }

  gkyl_free(mbf);
}
