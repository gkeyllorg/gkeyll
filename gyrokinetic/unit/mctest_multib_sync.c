#include <acutest.h>

#include <assert.h>

#include <gkyl_block_geom.h>
#include <gkyl_multib_comm_conn.h>
#include <gkyl_null_comm.h>
#include <gkyl_rrobin_decomp.h>

#ifdef GKYL_HAVE_MPI
#include <gkyl_mpi_comm.h>
#endif

#ifdef GKYL_HAVE_CUDA
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

struct gkyl_comm*
comm_new(bool use_mpi, bool use_gpu, FILE *iostream)
{
  // Construct communicator for use in app.
  struct gkyl_comm *comm = 0;

#ifdef GKYL_HAVE_MPI
  if (use_gpu && use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#else
    fprintf(iostream, " Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = use_gpu
    }
  );
#endif

  return comm;
}

static struct gkyl_block_geom *
create_L_domain_block_geom(int **cuts)
{
  // 2D with 3 blocks
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 3);
  
  /* Block layout

     +------+
     |0     |
     |      |
     +------+-----+
     |1     |2    |
     |      |     |
     +------+-----+
    
  */

  // block 0
  int *cuts0 = cuts[0];
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 300, 300 },
      .cuts = { cuts0[0], cuts0[1] },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
      },
      .connections[1] = { // y-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
      }
    }
  );
  
  // block 1
  int *cuts1 = cuts[1];
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 300, 300 },
      .cuts = { cuts1[0], cuts1[1] },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 2, .dir = 0, .edge = GKYL_LOWER_POSITIVE }
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE }
      }
    }
  );

  // block 2
  int *cuts2 = cuts[2];
  gkyl_block_geom_set_block(bgeom, 2, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 300, 300 },
      .cuts = { cuts2[0], cuts2[1] },
      
      .connections[0] = { // x-direction connections
        { .bid = 1, .dir = 0, .edge = GKYL_UPPER_POSITIVE },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL } // physical boundary
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL } // physical boundary
      }
    }
  );

  return bgeom;
}

static struct gkyl_block_geom *
create_cyclic_domain_block_geom(int **cuts)
{
  // 2D with 2 blocks
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 2);
  
  /* Block layout

     +------+
     |0     |
     |      |
     +------+
     |1     |
     |      |
     +------+
     |0     |
     |      |
     +------+
    
  */

  // block 0
  int *cuts0 = cuts[0];
  gkyl_block_geom_set_block(bgeom, 0, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 4, 8 },
      .cuts = { cuts0[0], cuts0[1] },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
      },
      .connections[1] = { // y-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 1, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );
  
  // block 1
  int *cuts1 = cuts[1];
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 4, 6 },
      .cuts = { cuts1[0], cuts1[1] },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  return bgeom;
}


static inline int
prod_of_elements_int(int ndim, int *arr)
{
  int pr = 1;
  for (int d=0; d<ndim; ++d) pr *= arr[d];
  return pr;
}

static bool
has_int(int n, int val, const int *lst)
{
  for (int i=0; i<n; ++i)
    if (val == lst[i])
      return true;
  return false;
}

struct app_L {
  struct gkyl_range global, global_ext;
  struct gkyl_range local, local_ext;
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];
  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];
  struct gkyl_basis basis;
  struct gkyl_array *f;
  struct gkyl_array *f_ho;
};

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static void
test_L_domain_sync(bool use_gpu, bool use_mpi, int **cuts, int poly_order)
{
  // Create world comm.
  struct gkyl_comm* comm = comm_new(use_mpi, use_gpu, stderr);

  int my_rank, num_ranks;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &num_ranks);

  int ndim = 2;

  struct gkyl_block_geom *geom = create_L_domain_block_geom(cuts);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  int num_blocks = topo->num_blocks;
  int nghost[ndim];
  for (int d=0; d<ndim; ++d)
    nghost[d] = 1;

  // Construct decomp objects.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    branks[i] = prod_of_elements_int(ndim, cuts[i]);
  }
  const struct gkyl_rrobin_decomp* round_robin_decomp = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // Allocate enough space.
  struct gkyl_rect_decomp **decomp = gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));

  int num_blocks_local = 0;
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(round_robin_decomp, i, rank_list);

    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    struct gkyl_range block_global_range;
    gkyl_create_global_range(ndim, ginfo->cells, &block_global_range);
    
    decomp[i] = gkyl_rect_decomp_new_from_cuts(ndim, ginfo->cuts, &block_global_range);

    if (my_rank == 0) {
      printf("b%d ranks: ",i);
      for (int d=0; d<branks[i]; ++d)
        printf(" %d",rank_list[d]);
      printf("\n");
    }

    if (has_int(branks[i], my_rank, rank_list)) {
      num_blocks_local += 1;
    }
  }

  // Store the index (for the list of global blocks) of the local blocks.
  int *local_blocks = gkyl_malloc(sizeof(int[num_blocks_local]));
  int lidx = 0;
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(round_robin_decomp, i, rank_list);
    if (has_int(branks[i], my_rank, rank_list)) {
      local_blocks[lidx++] = i;
    }
  }

  // Construct intrablock communicators.
  struct gkyl_comm **block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    bool status;
    block_comms[i] = gkyl_comm_create_comm_from_ranks(comm,
      branks[i], rank_list, decomp[i], &status);
  }

  // Create "apps", with their field and ranges.
  struct app_L **singleb_apps = gkyl_malloc(num_blocks_local*sizeof(struct app_L *));
  for (int bI=0; bI<num_blocks_local; ++bI) {
    singleb_apps[bI] = gkyl_malloc(sizeof(struct app_L));
    struct app_L *app = singleb_apps[bI];

    int bid = local_blocks[bI];
    gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    gkyl_create_ranges(&decomp[bid]->ranges[brank], nghost, &app->local_ext, &app->local);

    for (int dir=0; dir<ndim; ++dir) {
      gkyl_skin_ghost_ranges(&app->lower_skin[dir], &app->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->local_ext, nghost);
      gkyl_skin_ghost_ranges(&app->upper_skin[dir], &app->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->local_ext, nghost);
    }

    gkyl_cart_modal_serendip(&app->basis, ndim, poly_order);

    app->f = mkarr(use_gpu, app->basis.num_basis, app->local_ext.volume);
    app->f_ho = use_gpu? mkarr(false, app->basis.num_basis, app->local_ext.volume)
	               : gkyl_array_acquire(app->f);
    
    // Put some value in f that is distinct in every rank and block.
//    for (int k=0; k<app->basis.num_basis; k++)
    int k=0;
      gkyl_array_shiftc(app->f, bid+100.0*my_rank, k);
  }

  // Communication connections. 
  struct gkyl_multib_comm_conn **mbcc_recv = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_multib_comm_conn *));
  struct gkyl_multib_comm_conn **mbcc_send = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_multib_comm_conn *));

  // Array of local ranges and fields.
  struct gkyl_array **fs = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_array *));
  struct gkyl_range **locals = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_range *));
  struct gkyl_range **local_exts = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_range *));

  for (int bI=0; bI<num_blocks_local; ++bI) {
    int bid = local_blocks[bI];

    gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbcc_recv[bI] = gkyl_multib_comm_conn_new_recv(bid, brank, nghost,
      &topo->conn[bid], decomp);
    mbcc_send[bI] = gkyl_multib_comm_conn_new_send(bid, brank, nghost,
      &topo->conn[bid], decomp);

    struct app_L *app = singleb_apps[bI];

    for (int ns=0; ns<mbcc_recv[bI]->num_comm_conn; ++ns) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      int rankIdx = mbcc_recv[bI]->comm_conn[ns].rank;
      gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc_recv[bI]->comm_conn[ns].block_id, rank_list);
      mbcc_recv[bI]->comm_conn[ns].rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&mbcc_recv[bI]->comm_conn[ns].range, &app->local_ext,
        mbcc_recv[bI]->comm_conn[ns].range.lower, mbcc_recv[bI]->comm_conn[ns].range.upper);
    }
    for (int ns=0; ns<mbcc_send[bI]->num_comm_conn; ++ns) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      int rankIdx = mbcc_send[bI]->comm_conn[ns].rank;
      gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc_send[bI]->comm_conn[ns].block_id, rank_list);
      mbcc_send[bI]->comm_conn[ns].rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&mbcc_send[bI]->comm_conn[ns].range, &app->local_ext,
        mbcc_send[bI]->comm_conn[ns].range.lower, mbcc_send[bI]->comm_conn[ns].range.upper);
    }

    // Sort connections according to rank and block ID.
    gkyl_multib_comm_conn_sort(mbcc_recv[bI]);
    gkyl_multib_comm_conn_sort(mbcc_send[bI]);

    fs[bI] = app->f;
    locals[bI] = &app->local;
    local_exts[bI] = &app->local_ext;
  }

  // Sync blocks.
  gkyl_multib_comm_conn_array_transfer(comm, num_blocks_local, local_blocks,
    mbcc_send, mbcc_recv, fs, fs);

  // Check results.
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct app_L *app = singleb_apps[bI];
    int bid = local_blocks[bI];

    gkyl_array_copy(app->f_ho, app->f);

    for (int ns=0; ns<mbcc_recv[bI]->num_comm_conn; ++ns) {
      struct gkyl_comm_conn *cc = &mbcc_recv[bI]->comm_conn[ns];
      double ref = cc->block_id + 100.0*cc->rank;

      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, &cc->range);
      while (gkyl_range_iter_next(&iter)) {
        long linidx = gkyl_range_idx(&cc->range, iter.idx);
        double *f_c = gkyl_array_fetch(app->f_ho, linidx);
//        for (int k=0; k<app->basis.num_basis; k++) {
        int k=0;
          TEST_CHECK( gkyl_compare(ref, f_c[k], 1e-10) );
          TEST_MSG( "bid:%d | Expected: %.13e | Got: %.13e | Cell:%d,%d\n", bid, ref, f_c[k], iter.idx[0], iter.idx[1]);
//        }
      }
    }
  }

  gkyl_free(fs);
  gkyl_free(locals);
  gkyl_free(local_exts);
  for (int bI=0; bI<num_blocks_local; ++bI) {
    gkyl_multib_comm_conn_release(mbcc_send[bI]);
    gkyl_multib_comm_conn_release(mbcc_recv[bI]);
  }
  gkyl_free(mbcc_send);
  gkyl_free(mbcc_recv);

  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct app_L *app = singleb_apps[bI];
    gkyl_array_release(app->f_ho);
    gkyl_array_release(app->f);
    gkyl_free(app);
  }
  gkyl_free(singleb_apps);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(block_comms[i]);
  gkyl_free(block_comms);
  gkyl_free(local_blocks);
  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(rank_list);
  gkyl_free(decomp);
  gkyl_rrobin_decomp_release(round_robin_decomp);
  gkyl_free(branks);
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
  gkyl_comm_release(comm);
}

static void
test_cyclic_domain_sync(bool use_gpu, bool use_mpi, int **cuts, int poly_order)
{
  // Create world comm.
  struct gkyl_comm* comm = comm_new(use_mpi, use_gpu, stderr);

  int my_rank, num_ranks;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &num_ranks);

  int ndim = 2;

  struct gkyl_block_geom *geom = create_cyclic_domain_block_geom(cuts);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  int num_blocks = topo->num_blocks;
  int nghost[ndim];
  for (int d=0; d<ndim; ++d)
    nghost[d] = 1;

  // Construct decomp objects.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    branks[i] = prod_of_elements_int(ndim, cuts[i]);
  }
  const struct gkyl_rrobin_decomp* round_robin_decomp = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // Allocate enough space.
  struct gkyl_rect_decomp **decomp = gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));

  int num_blocks_local = 0;
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(round_robin_decomp, i, rank_list);

    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    struct gkyl_range block_global_range;
    gkyl_create_global_range(ndim, ginfo->cells, &block_global_range);
    
    decomp[i] = gkyl_rect_decomp_new_from_cuts(ndim, ginfo->cuts, &block_global_range);

    if (my_rank == 0) {
      printf("b%d ranks: ",i);
      for (int d=0; d<branks[i]; ++d)
        printf(" %d",rank_list[d]);
      printf("\n");
    }

    if (has_int(branks[i], my_rank, rank_list)) {
      num_blocks_local += 1;
    }
  }

  // Store the index (for the list of global blocks) of the local blocks.
  int *local_blocks = gkyl_malloc(sizeof(int[num_blocks_local]));
  int lidx = 0;
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(round_robin_decomp, i, rank_list);
    if (has_int(branks[i], my_rank, rank_list)) {
      local_blocks[lidx++] = i;
    }
  }

  // Construct intrablock communicators.
  struct gkyl_comm **block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    bool status;
    block_comms[i] = gkyl_comm_create_comm_from_ranks(comm,
      branks[i], rank_list, decomp[i], &status);
  }

  // Create "apps", with their field and ranges.
  struct app_L **singleb_apps = gkyl_malloc(num_blocks_local*sizeof(struct app_L *));
  for (int bI=0; bI<num_blocks_local; ++bI) {
    singleb_apps[bI] = gkyl_malloc(sizeof(struct app_L));
    struct app_L *app = singleb_apps[bI];

    int bid = local_blocks[bI];
    gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    gkyl_create_ranges(&decomp[bid]->ranges[brank], nghost, &app->local_ext, &app->local);

    for (int dir=0; dir<ndim; ++dir) {
      gkyl_skin_ghost_ranges(&app->lower_skin[dir], &app->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->local_ext, nghost);
      gkyl_skin_ghost_ranges(&app->upper_skin[dir], &app->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->local_ext, nghost);
    }

    gkyl_cart_modal_serendip(&app->basis, ndim, poly_order);

    app->f = mkarr(use_gpu, app->basis.num_basis, app->local_ext.volume);
    app->f_ho = use_gpu? mkarr(false, app->basis.num_basis, app->local_ext.volume)
	               : gkyl_array_acquire(app->f);
    
    // Put some value in f that is distinct in every rank and block.
//    for (int k=0; k<app->basis.num_basis; k++)
    int k=0;
      gkyl_array_shiftc(app->f, bid+1+100.0*(my_rank+1), k);
      printf("bid = %d\n", bid);
  }

  // Communication connections. 
  struct gkyl_multib_comm_conn **mbcc_recv = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_multib_comm_conn *));
  struct gkyl_multib_comm_conn **mbcc_send = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_multib_comm_conn *));

  // Array of local ranges and fields.
  struct gkyl_array **fs = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_array *));
  struct gkyl_range **locals = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_range *));
  struct gkyl_range **local_exts = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_range *));

  for (int bI=0; bI<num_blocks_local; ++bI) {
    int bid = local_blocks[bI];

    gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbcc_recv[bI] = gkyl_multib_comm_conn_new_recv(bid, brank, nghost,
      &topo->conn[bid], decomp);
    mbcc_send[bI] = gkyl_multib_comm_conn_new_send(bid, brank, nghost,
      &topo->conn[bid], decomp);

    struct app_L *app = singleb_apps[bI];

    for (int ns=0; ns<mbcc_recv[bI]->num_comm_conn; ++ns) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      int rankIdx = mbcc_recv[bI]->comm_conn[ns].rank;
      gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc_recv[bI]->comm_conn[ns].block_id, rank_list);
      mbcc_recv[bI]->comm_conn[ns].rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&mbcc_recv[bI]->comm_conn[ns].range, &app->local_ext,
        mbcc_recv[bI]->comm_conn[ns].range.lower, mbcc_recv[bI]->comm_conn[ns].range.upper);
    }
    for (int ns=0; ns<mbcc_send[bI]->num_comm_conn; ++ns) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      int rankIdx = mbcc_send[bI]->comm_conn[ns].rank;
      gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc_send[bI]->comm_conn[ns].block_id, rank_list);
      mbcc_send[bI]->comm_conn[ns].rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&mbcc_send[bI]->comm_conn[ns].range, &app->local_ext,
        mbcc_send[bI]->comm_conn[ns].range.lower, mbcc_send[bI]->comm_conn[ns].range.upper);
    }

    // Sort connections according to rank and block ID.
    gkyl_multib_comm_conn_sort(mbcc_recv[bI]);
    gkyl_multib_comm_conn_sort(mbcc_send[bI]);

    fs[bI] = app->f;
    locals[bI] = &app->local;
    local_exts[bI] = &app->local_ext;
  }

  // Sync blocks.
  gkyl_multib_comm_conn_array_transfer(comm, num_blocks_local, local_blocks,
    mbcc_send, mbcc_recv, fs, fs);

  // Check results.
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct app_L *app = singleb_apps[bI];
    int bid = local_blocks[bI];

    gkyl_array_copy(app->f_ho, app->f);

    for (int ns=0; ns<mbcc_recv[bI]->num_comm_conn; ++ns) {
      struct gkyl_comm_conn *cc = &mbcc_recv[bI]->comm_conn[ns];
      double ref = cc->block_id+1 + 100.0*(cc->rank+1);

      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, &cc->range);
      while (gkyl_range_iter_next(&iter)) {
        long linidx = gkyl_range_idx(&cc->range, iter.idx);
        double *f_c = gkyl_array_fetch(app->f_ho, linidx);
//        for (int k=0; k<app->basis.num_basis; k++) {
        int k=0;
          TEST_CHECK( gkyl_compare(ref, f_c[k], 1e-10) );
          TEST_MSG( "bid:%d | Expected: %.13e | Got: %.13e | Cell:%d,%d\n", bid, ref, f_c[k], iter.idx[0], iter.idx[1]);
//        }
      }
    }
  }

  gkyl_free(fs);
  gkyl_free(locals);
  gkyl_free(local_exts);
  for (int bI=0; bI<num_blocks_local; ++bI) {
    gkyl_multib_comm_conn_release(mbcc_send[bI]);
    gkyl_multib_comm_conn_release(mbcc_recv[bI]);
  }
  gkyl_free(mbcc_send);
  gkyl_free(mbcc_recv);

  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct app_L *app = singleb_apps[bI];
    gkyl_array_release(app->f_ho);
    gkyl_array_release(app->f);
    gkyl_free(app);
  }
  gkyl_free(singleb_apps);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(block_comms[i]);
  gkyl_free(block_comms);
  gkyl_free(local_blocks);
  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(rank_list);
  gkyl_free(decomp);
  gkyl_rrobin_decomp_release(round_robin_decomp);
  gkyl_free(branks);
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
  gkyl_comm_release(comm);
}

static void
test_cyclic_domain_sync_ser(bool use_gpu, bool use_mpi, int **cuts, int poly_order)
{
  // Create world comm.
  struct gkyl_comm* comm = comm_new(use_mpi, use_gpu, stderr);

  int my_rank, num_ranks;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &num_ranks);

  int ndim = 2;

  struct gkyl_block_geom *geom = create_cyclic_domain_block_geom(cuts);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  int num_blocks = topo->num_blocks;
  int nghost[ndim];
  for (int d=0; d<ndim; ++d)
    nghost[d] = 1;

  // Construct decomp objects.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    branks[i] = prod_of_elements_int(ndim, cuts[i]);
  }
  const struct gkyl_rrobin_decomp* round_robin_decomp = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // Allocate enough space.
  struct gkyl_rect_decomp **decomp = gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));

  int num_blocks_local = 0;
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(round_robin_decomp, i, rank_list);

    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    struct gkyl_range block_global_range;
    gkyl_create_global_range(ndim, ginfo->cells, &block_global_range);
    
    decomp[i] = gkyl_rect_decomp_new_from_cuts(ndim, ginfo->cuts, &block_global_range);

    if (my_rank == 0) {
      printf("b%d ranks: ",i);
      for (int d=0; d<branks[i]; ++d)
        printf(" %d",rank_list[d]);
      printf("\n");
    }

    if (has_int(branks[i], my_rank, rank_list)) {
      num_blocks_local += 1;
    }
  }

  // Store the index (for the list of global blocks) of the local blocks.
  int *local_blocks = gkyl_malloc(sizeof(int[num_blocks_local]));
  int lidx = 0;
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(round_robin_decomp, i, rank_list);
    if (has_int(branks[i], my_rank, rank_list)) {
      local_blocks[lidx++] = i;
    }
  }

  // Construct intrablock communicators.
  struct gkyl_comm **block_comms = gkyl_malloc(num_blocks*sizeof(struct gkyl_comm *));
  for (int i=0; i<num_blocks; ++i) {
    bool status;
    block_comms[i] = gkyl_comm_create_comm_from_ranks(comm,
      branks[i], rank_list, decomp[i], &status);
  }

  // Create "apps", with their field and ranges.
  struct app_L **singleb_apps = gkyl_malloc(num_blocks_local*sizeof(struct app_L *));
  for (int bI=0; bI<num_blocks_local; ++bI) {
    singleb_apps[bI] = gkyl_malloc(sizeof(struct app_L));
    struct app_L *app = singleb_apps[bI];

    int bid = local_blocks[bI];
    gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    gkyl_create_ranges(&decomp[bid]->ranges[brank], nghost, &app->local_ext, &app->local);

    for (int dir=0; dir<ndim; ++dir) {
      gkyl_skin_ghost_ranges(&app->lower_skin[dir], &app->lower_ghost[dir], dir, GKYL_LOWER_EDGE, &app->local_ext, nghost);
      gkyl_skin_ghost_ranges(&app->upper_skin[dir], &app->upper_ghost[dir], dir, GKYL_UPPER_EDGE, &app->local_ext, nghost);
    }

    gkyl_cart_modal_serendip(&app->basis, ndim, poly_order);

    app->f = mkarr(use_gpu, app->basis.num_basis, app->local_ext.volume);
    app->f_ho = use_gpu? mkarr(false, app->basis.num_basis, app->local_ext.volume)
	               : gkyl_array_acquire(app->f);
    
    // Put some value in f that is distinct in every rank and block.
//    for (int k=0; k<app->basis.num_basis; k++)
    int k=0;
      //gkyl_array_shiftc(app->f, bid+1+100.0*(my_rank+1), k);
      //printf("bid = %d\n", bid);

    int sublower[2] = {app->local.lower[0], app->local.lower[1]};
    int subupper[2] = {app->local.upper[0], app->local.upper[1]/2};
    struct gkyl_range lower_range;
    gkyl_sub_range_init(&lower_range, &app->local_ext, sublower, subupper);
    sublower[1] = app->local.upper[1]/2;
    subupper[1] = app->local.upper[1];
    struct gkyl_range upper_range;
    gkyl_sub_range_init(&upper_range, &app->local_ext, sublower, subupper);
    gkyl_array_shiftc(app->f, 1.0, k);
    double val_lower = 2*bid +1;
    double val_upper= 2*bid +2;
    gkyl_array_set_range(app->f, val_lower, app->f, &lower_range);
    gkyl_array_set_range(app->f, val_upper, app->f, &upper_range);
  }

  // Communication connections. 
  struct gkyl_multib_comm_conn **mbcc_recv = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_multib_comm_conn *));
  struct gkyl_multib_comm_conn **mbcc_send = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_multib_comm_conn *));

  // Array of local ranges and fields.
  struct gkyl_array **fs = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_array *));
  struct gkyl_range **locals = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_range *));
  struct gkyl_range **local_exts = gkyl_malloc(num_blocks_local * sizeof(struct gkyl_range *));

  for (int bI=0; bI<num_blocks_local; ++bI) {
    int bid = local_blocks[bI];

    gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbcc_recv[bI] = gkyl_multib_comm_conn_new_recv(bid, brank, nghost,
      &topo->conn[bid], decomp);
    mbcc_send[bI] = gkyl_multib_comm_conn_new_send(bid, brank, nghost,
      &topo->conn[bid], decomp);

    struct app_L *app = singleb_apps[bI];

    for (int ns=0; ns<mbcc_recv[bI]->num_comm_conn; ++ns) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      int rankIdx = mbcc_recv[bI]->comm_conn[ns].rank;
      gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc_recv[bI]->comm_conn[ns].block_id, rank_list);
      mbcc_recv[bI]->comm_conn[ns].rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&mbcc_recv[bI]->comm_conn[ns].range, &app->local_ext,
        mbcc_recv[bI]->comm_conn[ns].range.lower, mbcc_recv[bI]->comm_conn[ns].range.upper);
    }
    for (int ns=0; ns<mbcc_send[bI]->num_comm_conn; ++ns) {
      // Translate the "rank" in gkyl_multib_comm_conn (right now it is a rank index).
      int rankIdx = mbcc_send[bI]->comm_conn[ns].rank;
      gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc_send[bI]->comm_conn[ns].block_id, rank_list);
      mbcc_send[bI]->comm_conn[ns].rank = rank_list[rankIdx];
      // Make range a sub range.
      gkyl_sub_range_init(&mbcc_send[bI]->comm_conn[ns].range, &app->local_ext,
        mbcc_send[bI]->comm_conn[ns].range.lower, mbcc_send[bI]->comm_conn[ns].range.upper);
    }

    // Sort connections according to rank and block ID.
    gkyl_multib_comm_conn_sort(mbcc_recv[bI]);
    gkyl_multib_comm_conn_sort(mbcc_send[bI]);

    fs[bI] = app->f;
    locals[bI] = &app->local;
    local_exts[bI] = &app->local_ext;
  }

  // Sync blocks.
  gkyl_multib_comm_conn_array_transfer(comm, num_blocks_local, local_blocks,
    mbcc_send, mbcc_recv, fs, fs);

  // Check results.
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct app_L *app = singleb_apps[bI];
    int bid = local_blocks[bI];

    gkyl_array_copy(app->f_ho, app->f);

    for (int ns=0; ns<mbcc_recv[bI]->num_comm_conn; ++ns) {
      struct gkyl_comm_conn *cc = &mbcc_recv[bI]->comm_conn[ns];
      double ref_lower = 2*cc->block_id+1;
      double ref_upper = 2*cc->block_id+2;

      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, &cc->range);
      while (gkyl_range_iter_next(&iter)) {
        double ref;
        if(iter.idx[1] == 0)
          ref = ref_upper;
        else
          ref = ref_lower;
                  
        long linidx = gkyl_range_idx(&cc->range, iter.idx);
        double *f_c = gkyl_array_fetch(app->f_ho, linidx);
//        for (int k=0; k<app->basis.num_basis; k++) {
        int k=0;
          TEST_CHECK( gkyl_compare(ref, f_c[k], 1e-10) );
          TEST_MSG( "bid:%d | Expected: %.13e | Got: %.13e | Cell:%d,%d\n", bid, ref, f_c[k], iter.idx[0], iter.idx[1]);
//        }
      }
    }
  }

  gkyl_free(fs);
  gkyl_free(locals);
  gkyl_free(local_exts);
  for (int bI=0; bI<num_blocks_local; ++bI) {
    gkyl_multib_comm_conn_release(mbcc_send[bI]);
    gkyl_multib_comm_conn_release(mbcc_recv[bI]);
  }
  gkyl_free(mbcc_send);
  gkyl_free(mbcc_recv);

  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct app_L *app = singleb_apps[bI];
    gkyl_array_release(app->f_ho);
    gkyl_array_release(app->f);
    gkyl_free(app);
  }
  gkyl_free(singleb_apps);

  for (int i=0; i<num_blocks; ++i)
    gkyl_comm_release(block_comms[i]);
  gkyl_free(block_comms);
  gkyl_free(local_blocks);
  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(rank_list);
  gkyl_free(decomp);
  gkyl_rrobin_decomp_release(round_robin_decomp);
  gkyl_free(branks);
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
  gkyl_comm_release(comm);
}


int **
cuts_array_new(int num_blocks, int ndim, int *cuts_all)
{
  // Create an array of cuts from an array with all the cuts listed flat.
  int **cuts_arr = gkyl_malloc(num_blocks * sizeof(int *));
  for (int i=0; i<num_blocks; i++) {
    cuts_arr[i] = gkyl_malloc(ndim * sizeof(int));

    int *cuts = cuts_arr[i];
    for (int d=0; d<ndim; d++)
      cuts[d] = cuts_all[i*ndim + d];
  }
  return cuts_arr;
}

void
cuts_array_release(int num_blocks, int **cuts_arr)
{
  // Release the array of cuts arrays.
  for (int i=0; i<num_blocks; i++)
    gkyl_free(cuts_arr[i]);
  gkyl_free(cuts_arr);
}

static void
test_L_domain_sync_ho(void)
{
  int num_blocks = 3; // L-shaped example.
  int ndim = 2;

//  int cuts_flat0[] = {
//    1, 1, // Block 0.
//    1, 1, // Block 1.
//    1, 1, // Block 2.
//  };
//  int **cuts0 = cuts_array_new(num_blocks, ndim, cuts_flat0);
//  test_L_domain_sync(false, true, cuts0, 1);
//  cuts_array_release(num_blocks, cuts0);

  int cuts_flat1[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    1, 1, // Block 2.
  };
  int **cuts1 = cuts_array_new(num_blocks, ndim, cuts_flat1);
  test_L_domain_sync(false, true, cuts1, 1);
  cuts_array_release(num_blocks, cuts1);
}

static void
test_cyclic_domain_sync_ho(void)
{
  int num_blocks = 2; // cyclic-shaped example.
  int ndim = 2;

  int cuts_flat1[] = {
    1, 2, // Block 0.
    1, 2, // Block 1.
  };
  int **cuts1 = cuts_array_new(num_blocks, ndim, cuts_flat1);
  test_cyclic_domain_sync(false, true, cuts1, 1);
  cuts_array_release(num_blocks, cuts1);
}

static void
test_cyclic_domain_sync_ser_ho(void)
{
  int num_blocks = 2; // cyclic-shaped example.
  int ndim = 2;

  int cuts_flat1[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
  };
  int **cuts1 = cuts_array_new(num_blocks, ndim, cuts_flat1);
  // Test using MPI
  test_cyclic_domain_sync_ser(false, true, cuts1, 1);
  // Test without MPI
  test_cyclic_domain_sync_ser(false, false, cuts1, 1);
  cuts_array_release(num_blocks, cuts1);
}

#ifdef GKYL_HAVE_CUDA
static void
test_L_domain_sync_dev(void)
{
  int num_blocks = 3; // L-shaped example.
  int ndim = 2;

  int cuts_flat0[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    1, 1, // Block 2.
  };
  int **cuts0 = cuts_array_new(num_blocks, ndim, cuts_flat0);
  #ifdef GKYL_HAVE_NCCL
    test_L_domain_sync(true, true, cuts0, 1);
  #endif
  cuts_array_release(num_blocks, cuts0);
}

static void
test_cyclic_domain_sync_dev(void)
{
  int num_blocks = 2; // cyclic-shaped example.
  int ndim = 2;

  int cuts_flat0[] = {
    1, 2, // Block 0.
    1, 2, // Block 1.
  };
  int **cuts0 = cuts_array_new(num_blocks, ndim, cuts_flat0);
  #ifdef GKYL_HAVE_NCCL
    test_cyclic_domain_sync(true, true, cuts0, 1);
  #endif
  cuts_array_release(num_blocks, cuts0);
}

static void
test_cyclic_domain_sync_ser_dev(void)
{
  int num_blocks = 2; // cyclic-shaped example.
  int ndim = 2;

  int cuts_flat0[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
  };
  int **cuts0 = cuts_array_new(num_blocks, ndim, cuts_flat0);
  // Test using MPI

  #ifdef GKYL_HAVE_NCCL
    test_cyclic_domain_sync_ser(true, true, cuts0, 1);
  #endif
  // Test without MPI
  test_cyclic_domain_sync_ser(true, false, cuts0, 1);
  cuts_array_release(num_blocks, cuts0);
}
#endif

TEST_LIST = {
  { "test_L_domain_sync_ho" , test_L_domain_sync_ho },
  { "test_cyclic_domain_sync_ho" , test_cyclic_domain_sync_ho },
  { "test_cyclic_domain_sync_ser_ho" , test_cyclic_domain_sync_ser_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_L_domain_sync_dev", test_L_domain_sync_dev },
  { "test_cyclic_domain_sync_dev", test_cyclic_domain_sync_dev },
  { "test_cyclic_domain_sync_ser_dev", test_cyclic_domain_sync_ser_dev },
#endif
  { NULL, NULL },
};
