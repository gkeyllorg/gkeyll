#include <acutest.h>

#include <gkyl_block_geom.h>
#include <gkyl_multib_comm_conn.h>
#include <gkyl_multib_conn.h>
#include <gkyl_rrobin_decomp.h>
#include <string.h>

#ifdef GKYL_HAVE_MPI
#include <gkyl_mpi_comm.h>
#endif

#ifdef GKYL_HAVE_CUDA
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

// allocate array (filled with zeros)
static struct gkyl_array*
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array *a = use_gpu? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size)
                                : gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct gkyl_comm*
comm_new(bool use_mpi, bool use_gpu, FILE *iostream)
{
  // Construct communicator for use in app.
  struct gkyl_comm *comm = 0;

#ifdef GKYL_HAVE_MPI
  if (use_gpu && use_mpi) {
#ifdef GKYL_HAVE_NCCL
    printf("using NCCL\n");
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
create_SOL_domain_block_geom(int **cuts)
{
  // 2D with 3 blocks
  struct gkyl_block_geom *bgeom = gkyl_block_geom_new(2, 3);
  
  /* Block layout

     +------+
     |2     |
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
      .cells = { 4, 8},
      .cuts = { cuts0[0], cuts0[1] },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }, // physical boundary
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }  // physical boundary
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL},
        { .bid = 1, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );
  
  // block 1
  int *cuts1 = cuts[1];
  gkyl_block_geom_set_block(bgeom, 1, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 4, 8},
      .cuts = { cuts1[0], cuts1[1] },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL },
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL}
      },
      .connections[1] = { // y-direction connections
        { .bid = 0, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 2, .dir = 1, .edge = GKYL_LOWER_POSITIVE}
      }
    }
  );

  // block 2
  int *cuts2 = cuts[2];
  gkyl_block_geom_set_block(bgeom, 2, &(struct gkyl_block_geom_info) {
      .lower = { 0, 0 },
      .upper = { 1, 1 },
      .cells = { 4, 8},
      .cuts = { cuts2[0], cuts2[1] },
      
      .connections[0] = { // x-direction connections
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL},
        { .bid = 0, .dir = 0, .edge = GKYL_PHYSICAL }
      },
      .connections[1] = { // y-direction connections
        { .bid = 1, .dir = 1, .edge = GKYL_UPPER_POSITIVE},
        { .bid = 0, .dir = 1, .edge = GKYL_PHYSICAL }
      }
    }
  );

  return bgeom;
}




// simple linear search to check if val occurs in lst
static bool
has_int(int n, int val, const int *lst)
{
  for (int i=0; i<n; ++i)
    if (val == lst[i])
      return true;
  return false;
}

static inline int
prod_of_elements_int(int ndim, int *arr)
{
  int pr = 1;
  for (int d=0; d<ndim; ++d) pr *= arr[d];
  return pr;
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
test_L_domain_send_connections_dir0_cuts1()
{
  int num_blocks = 3; // L-shaped example.
  int ndim = 2;
  int cuts_flat1[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    1, 1, // Block 2.
  };
  int **cuts1 = cuts_array_new(num_blocks, ndim, cuts_flat1);
  struct gkyl_block_geom *geom  = create_L_domain_block_geom(cuts1);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);
  cuts_array_release(num_blocks, cuts1);

  int num_cuts[num_blocks];
  int nghost[] = { 1, 1 };

  // Setup for a gather along x
  int block_list[3][2] = {{0},{1,2},{1,2}};
  int dir = 0;
  int nconnected[3] = {1,2,2};
  
  // construct decomp objects
  struct gkyl_rect_decomp **decomp =
    gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    num_cuts[i] = 1;
    for (int d=0; d<topo->ndim; ++d)
      num_cuts[i] *= ginfo->cuts[d];
    
    struct gkyl_range range;
    gkyl_create_global_range(2, ginfo->cells, &range);
    decomp[i] = gkyl_rect_decomp_new_from_cuts(2, ginfo->cuts, &range);
  }

  // for testing
  int num_send_neigh[] = { 1, 2, 2 };

  // for testing (these hard-coded values depend on how the algorithm
  // is implemented)  
  struct gkyl_comm_conn conn_0[] = {
    { .block_id = 0, .rank = 0 },
  };
  gkyl_range_init(&conn_0[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });

  struct gkyl_comm_conn conn_1[] = {
    { .block_id = 1, .rank = 0 },
    { .block_id = 2, .rank = 0 },    
  };
  gkyl_range_init(&conn_1[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[1].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  
  struct gkyl_comm_conn conn_2[] = {
    { .block_id = 1, .rank = 0 },
    { .block_id = 2, .rank = 0 },
  };
  gkyl_range_init(&conn_2[0].range, 2, (int[]) { 301, 1 }, (int[]) { 600, 300 });
  gkyl_range_init(&conn_2[1].range, 2, (int[]) { 301, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn *block_conn[] = { conn_0, conn_1, conn_2 };
  
  for (int bid=0; bid<num_blocks; ++bid) {
    for (int brank=0; brank<num_cuts[bid]; ++brank) {
      struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new_send_from_connections(bid, brank, nghost, nconnected[bid], block_list[bid], dir, decomp);

      TEST_CHECK( num_send_neigh[bid] == mbcc->num_comm_conn );
      if (mbcc->num_comm_conn > 0) {
        for (int ns=0; ns<mbcc->num_comm_conn; ++ns) {
          TEST_CHECK( block_conn[bid][ns].block_id == mbcc->comm_conn[ns].block_id);
          TEST_CHECK( block_conn[bid][ns].rank == mbcc->comm_conn[ns].rank);
          TEST_CHECK( gkyl_range_compare(&block_conn[bid][ns].range, &mbcc->comm_conn[ns].range) );
        }
      }

      gkyl_multib_comm_conn_release(mbcc);
    }
  }

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(decomp);
  
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
}

static void
test_L_domain_recv_connections_dir0_cuts1()
{
  int num_blocks = 3; // L-shaped example.
  int ndim = 2;
  int cuts_flat1[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    1, 1, // Block 2.
  };
  int **cuts1 = cuts_array_new(num_blocks, ndim, cuts_flat1);
  struct gkyl_block_geom *geom  = create_L_domain_block_geom(cuts1);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);
  cuts_array_release(num_blocks, cuts1);

  int num_cuts[num_blocks];
  int nghost[] = { 1, 1 };

  // Setup for a gather along x
  int block_list[3][2] = {{0},{1,2},{1,2}};
  int dir = 0;
  int nconnected[3] = {1,2,2};
  
  // construct decomp objects
  struct gkyl_rect_decomp **decomp =
    gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    num_cuts[i] = 1;
    for (int d=0; d<topo->ndim; ++d)
      num_cuts[i] *= ginfo->cuts[d];
    
    struct gkyl_range range;
    gkyl_create_global_range(2, ginfo->cells, &range);
    decomp[i] = gkyl_rect_decomp_new_from_cuts(2, ginfo->cuts, &range);
  }

  // for testing
  int num_recv_neigh[] = { 1, 2, 2 };

  // for testing (these hard-coded values depend on how the algorithm
  // is implemented)  
  struct gkyl_comm_conn conn_0[] = { 
    { .block_id = 0, .rank = 0 },
  };
  gkyl_range_init(&conn_0[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });

  struct gkyl_comm_conn conn_1[] = {
    { .block_id = 1, .rank = 0 },    
    { .block_id = 2, .rank = 0 },    
  };
  gkyl_range_init(&conn_1[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[1].range, 2, (int[]) { 301, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn conn_2[] = {
    { .block_id = 1, .rank = 0 },
    { .block_id = 2, .rank = 0 },
  };
  gkyl_range_init(&conn_2[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_2[1].range, 2, (int[]) { 301, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn *block_conn[] = { conn_0, conn_1, conn_2 };
  
  for (int bid=0; bid<num_blocks; ++bid) {
    for (int brank=0; brank<num_cuts[bid]; ++brank) {
      struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new_recv_from_connections(bid, brank, nghost, nconnected[bid], block_list[bid], dir, decomp);

      TEST_CHECK( num_recv_neigh[bid] == mbcc->num_comm_conn );
      if (mbcc->num_comm_conn > 0) {
        for (int ns=0; ns<mbcc->num_comm_conn; ++ns) {
          TEST_CHECK( block_conn[bid][ns].block_id == mbcc->comm_conn[ns].block_id);
          TEST_CHECK( block_conn[bid][ns].rank == mbcc->comm_conn[ns].rank);
          TEST_CHECK( gkyl_range_compare(&block_conn[bid][ns].range, &mbcc->comm_conn[ns].range) );
        }
      }

      gkyl_multib_comm_conn_release(mbcc);
    }
  }

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(decomp);
  
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
}

static void
test_L_domain_send_connections_dir0_cuts2()
{
  int num_blocks = 3; // L-shaped example.
  int ndim = 2;
  int cuts_flat1[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    2, 1, // Block 2.
  };
  int **cuts1 = cuts_array_new(num_blocks, ndim, cuts_flat1);
  struct gkyl_block_geom *geom  = create_L_domain_block_geom(cuts1);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  // Construct decomp objects.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    branks[i] = prod_of_elements_int(ndim, cuts1[i]);
  }
  cuts_array_release(num_blocks, cuts1);
  int num_ranks = 1;
  const struct gkyl_rrobin_decomp* round_robin_decomp = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int num_cuts[num_blocks];
  int nghost[] = { 1, 1 };

  // Setup for a gather along x
  int block_list[3][2] = {{0},{1,2},{1,2}};
  int dir = 0;
  int nconnected[3] = {1,2,2};
  
  // construct decomp objects
  struct gkyl_rect_decomp **decomp =
    gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    num_cuts[i] = 1;
    for (int d=0; d<topo->ndim; ++d)
      num_cuts[i] *= ginfo->cuts[d];
    
    struct gkyl_range range;
    gkyl_create_global_range(2, ginfo->cells, &range);
    decomp[i] = gkyl_rect_decomp_new_from_cuts(2, ginfo->cuts, &range);
  }


  // for testing
  int num_send_neigh[] = { 1, 3, 3 };

  // for testing (these hard-coded values depend on how the algorithm
  // is implemented)  
  struct gkyl_comm_conn conn_0[] = { 
    { .block_id = 0, .rank = 0 },
  };
  gkyl_range_init(&conn_0[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });

  struct gkyl_comm_conn conn_1[] = {
    { .block_id = 1, .rank = 0 },    // from 0th cut
    { .block_id = 2, .rank = 0 },    // from 0th cut
    { .block_id = 2, .rank = 0 },    // from 0th cut
  };
  gkyl_range_init(&conn_1[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[1].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[2].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  
  struct gkyl_comm_conn conn_2[] = {
    { .block_id = 1, .rank = 0 }, // from 0th cut
    { .block_id = 2, .rank = 0 }, // from 0th cut
    { .block_id = 2, .rank = 0 }, // from 0th cut
    { .block_id = 1, .rank = 0 }, // from 1st cut
    { .block_id = 2, .rank = 0 }, // from 1st cut
    { .block_id = 2, .rank = 0 }, // from 1st cut
  };
  gkyl_range_init(&conn_2[0].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[1].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[2].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[3].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  gkyl_range_init(&conn_2[4].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  gkyl_range_init(&conn_2[5].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn *block_conn[] = { conn_0, conn_1, conn_2 };

  for (int bid=0; bid<num_blocks; ++bid) {
    int start_ns= 0;
    for (int brank=0; brank<num_cuts[bid]; ++brank) {
      struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new_send_from_connections(bid, brank, nghost, nconnected[bid], block_list[bid], dir, decomp);
      TEST_CHECK( num_send_neigh[bid] == mbcc->num_comm_conn );
      if (mbcc->num_comm_conn > 0) {
        for (int ns=0; ns<mbcc->num_comm_conn; ++ns) {
          // need to get the actual rank that owns this cut
          int rank_list[decomp[mbcc->comm_conn[ns].block_id]->ndecomp]; // max number of ranks
          int rank_idx = mbcc->comm_conn[ns].rank;
          gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc->comm_conn[ns].block_id, rank_list);
          mbcc->comm_conn[ns].rank = rank_list[rank_idx];
          TEST_CHECK( block_conn[bid][start_ns+ns].block_id == mbcc->comm_conn[ns].block_id);
          TEST_CHECK( block_conn[bid][start_ns+ns].rank == mbcc->comm_conn[ns].rank);
          TEST_CHECK( gkyl_range_compare(&block_conn[bid][start_ns+ns].range, &mbcc->comm_conn[ns].range) );
        }
      }

      start_ns+=mbcc->num_comm_conn;
      gkyl_multib_comm_conn_release(mbcc);
    }
  }

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(decomp);
  
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
}

static void
test_L_domain_send_connections_dir0_cuts2_par()
{
  printf("\n");
  // Create world comm.
  bool use_mpi = true;
  bool use_gpu = false;
  struct gkyl_comm* comm = comm_new(use_mpi, use_gpu, stderr);

  int my_rank, num_ranks;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &num_ranks);
  if (num_ranks!=2) {return;}


  int num_blocks = 3; // L-shaped example.
  int ndim = 2;
  int cuts_flat[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    2, 1, // Block 2.
  };
  int **cuts = cuts_array_new(num_blocks, ndim, cuts_flat);
  struct gkyl_block_geom *geom  = create_L_domain_block_geom(cuts);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  // Construct decomp objects.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    branks[i] = prod_of_elements_int(ndim, cuts[i]);
  }
  cuts_array_release(num_blocks, cuts);
  const struct gkyl_rrobin_decomp* round_robin_decomp = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  // Check rank info. This is what was used to produce the expected results
  //int num_ranks_pb[3] = {-1};
  //num_ranks_pb[0] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 0);
  //num_ranks_pb[1] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 1);
  //num_ranks_pb[2] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 2);

  //int ranks_0[num_ranks_pb[0]];
  //int ranks_1[num_ranks_pb[1]];
  //int ranks_2[num_ranks_pb[2]];
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 0, ranks_0);
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 1, ranks_1);
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 2, ranks_2);

  //printf("num_ranks_pb = [%d, %d, %d]\n", num_ranks_pb[0], num_ranks_pb[1], num_ranks_pb[2]);
  //printf("Block 0 has ranks: ");
  //for (int i=0; i<num_ranks_pb[0]; i++) printf("%d ", ranks_0[i]);
  //printf("\n");
  //printf("Block 1 has ranks: ");
  //for (int i=0; i<num_ranks_pb[1]; i++) printf("%d ", ranks_1[i]);
  //printf("\n");
  //printf("Block 2 has ranks: ");
  //for (int i=0; i<num_ranks_pb[2]; i++) printf("%d ", ranks_2[i]);
  //printf("\n");


  int num_cuts[num_blocks];
  int nghost[] = { 1, 1 };

  // Setup for a gather along x
  int block_list[3][2] = {{0},{1,2},{1,2}};
  int dir = 0;
  int nconnected[3] = {1,2,2};
  
  // construct decomp objects
  struct gkyl_rect_decomp **decomp =
    gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    num_cuts[i] = 1;
    for (int d=0; d<topo->ndim; ++d)
      num_cuts[i] *= ginfo->cuts[d];
    
    struct gkyl_range range;
    gkyl_create_global_range(2, ginfo->cells, &range);
    decomp[i] = gkyl_rect_decomp_new_from_cuts(2, ginfo->cuts, &range);
  }


  // for testing
  int num_send_neigh[] = { 1, 3, 3 };

  // for testing (these hard-coded values depend on how the algorithm
  // is implemented)  
  struct gkyl_comm_conn conn_0[] = { 
    { .block_id = 0, .rank = 0 },
  };
  gkyl_range_init(&conn_0[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });

  struct gkyl_comm_conn conn_1[] = {
    { .block_id = 1, .rank = 1 },    // from 0th cut
    { .block_id = 2, .rank = 0 },    // from 0th cut to 0th cut
    { .block_id = 2, .rank = 1 },    // from 0th cut to 1st cut
  };
  gkyl_range_init(&conn_1[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[1].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[2].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  
  struct gkyl_comm_conn conn_2[] = {
    { .block_id = 1, .rank = 1 }, // from 0th cut
    { .block_id = 2, .rank = 0 }, // from 0th cut to 0th cut
    { .block_id = 2, .rank = 1 }, // from 0th cut to 1st cut
    { .block_id = 1, .rank = 1 }, // from 1st cut to 0th cut
    { .block_id = 2, .rank = 0 }, // from 1st cut to 0th cut
    { .block_id = 2, .rank = 1 }, // from 1st cut to 1st cut
  };
  gkyl_range_init(&conn_2[0].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[1].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[2].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[3].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  gkyl_range_init(&conn_2[4].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  gkyl_range_init(&conn_2[5].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn *block_conn[] = { conn_0, conn_1, conn_2 };

  for (int bid=0; bid<num_blocks; ++bid) {
    int start_ns= 0;
    for (int brank=0; brank<num_cuts[bid]; ++brank) {
      struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new_send_from_connections(bid, brank, nghost, nconnected[bid], block_list[bid], dir, decomp);
      TEST_CHECK( num_send_neigh[bid] == mbcc->num_comm_conn );
      for (int ns=0; ns<mbcc->num_comm_conn; ++ns) {
        // need to get the actual rank that owns this cut
        int rank_list[gkyl_rrobin_decomp_nranks(round_robin_decomp, mbcc->comm_conn[ns].block_id)];
        int rank_idx = mbcc->comm_conn[ns].rank;
        gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc->comm_conn[ns].block_id, rank_list);
        mbcc->comm_conn[ns].rank = rank_list[rank_idx];
        TEST_CHECK( block_conn[bid][start_ns+ns].block_id == mbcc->comm_conn[ns].block_id);
        TEST_CHECK( block_conn[bid][start_ns+ns].rank == mbcc->comm_conn[ns].rank);
        TEST_CHECK( gkyl_range_compare(&block_conn[bid][start_ns+ns].range, &mbcc->comm_conn[ns].range) );
      }

      start_ns+=mbcc->num_comm_conn;
      gkyl_multib_comm_conn_release(mbcc);
    }
  }

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(decomp);
  
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
  gkyl_comm_release(comm);

}

static void
test_L_domain_recv_connections_dir0_cuts2_par()
{
  printf("\n");
  // Create world comm.
  bool use_mpi = true;
  bool use_gpu = false;
  struct gkyl_comm* comm = comm_new(use_mpi, use_gpu, stderr);

  int my_rank, num_ranks;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &num_ranks);
  if (num_ranks!=2) {return;}


  int num_blocks = 3; // L-shaped example.
  int ndim = 2;
  int cuts_flat[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    2, 1, // Block 2.
  };
  int **cuts = cuts_array_new(num_blocks, ndim, cuts_flat);
  struct gkyl_block_geom *geom  = create_L_domain_block_geom(cuts);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  // Construct decomp objects.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    branks[i] = prod_of_elements_int(ndim, cuts[i]);
  }
  cuts_array_release(num_blocks, cuts);
  const struct gkyl_rrobin_decomp* round_robin_decomp = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  // Check rank info. This is what was used to produce the expected results
  //int num_ranks_pb[3] = {-1};
  //num_ranks_pb[0] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 0);
  //num_ranks_pb[1] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 1);
  //num_ranks_pb[2] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 2);

  //int ranks_0[num_ranks_pb[0]];
  //int ranks_1[num_ranks_pb[1]];
  //int ranks_2[num_ranks_pb[2]];
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 0, ranks_0);
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 1, ranks_1);
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 2, ranks_2);

  //printf("num_ranks_pb = [%d, %d, %d]\n", num_ranks_pb[0], num_ranks_pb[1], num_ranks_pb[2]);
  //printf("Block 0 has ranks: ");
  //for (int i=0; i<num_ranks_pb[0]; i++) printf("%d ", ranks_0[i]);
  //printf("\n");
  //printf("Block 1 has ranks: ");
  //for (int i=0; i<num_ranks_pb[1]; i++) printf("%d ", ranks_1[i]);
  //printf("\n");
  //printf("Block 2 has ranks: ");
  //for (int i=0; i<num_ranks_pb[2]; i++) printf("%d ", ranks_2[i]);
  //printf("\n");


  int num_cuts[num_blocks];
  int nghost[] = { 1, 1 };

  // Setup for a gather along x
  int block_list[3][2] = {{0},{1,2},{1,2}};
  int dir = 0;
  int nconnected[3] = {1,2,2};
  
  // construct decomp objects
  struct gkyl_rect_decomp **decomp =
    gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    num_cuts[i] = 1;
    for (int d=0; d<topo->ndim; ++d)
      num_cuts[i] *= ginfo->cuts[d];
    
    struct gkyl_range range;
    gkyl_create_global_range(2, ginfo->cells, &range);
    decomp[i] = gkyl_rect_decomp_new_from_cuts(2, ginfo->cuts, &range);
  }


  // for testing
  int num_recv_neigh[] = { 1, 3, 3 };

  // for testing (these hard-coded values depend on how the algorithm
  // is implemented)  
  struct gkyl_comm_conn conn_0[] = { 
    { .block_id = 0, .rank = 0 },
  };
  gkyl_range_init(&conn_0[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });

  struct gkyl_comm_conn conn_1[] = {
    { .block_id = 1, .rank = 1 },    // into 0th from 0th cut
    { .block_id = 2, .rank = 0 },    // into 0th cut from 0th cut
    { .block_id = 2, .rank = 1 },    // into 0th cut from 1st cut
  };
  gkyl_range_init(&conn_1[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_1[1].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_1[2].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn conn_2[] = {
    { .block_id = 1, .rank = 1 }, // into 0th cut from 0th cut
    { .block_id = 2, .rank = 0 }, // into 0th cut from 0th cut
    { .block_id = 2, .rank = 1 }, // into 0th cut from 1st cut
    { .block_id = 1, .rank = 1 }, // into 1st cut from 0th cut
    { .block_id = 2, .rank = 0 }, // into 1st cut from 0th cut
    { .block_id = 2, .rank = 1 }, // into 1st cut from 1st cut
  };
  gkyl_range_init(&conn_2[0].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_2[1].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[2].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  gkyl_range_init(&conn_2[3].range, 2, (int[]) { 1, 1 }, (int[]) { 300, 300 });
  gkyl_range_init(&conn_2[4].range, 2, (int[]) { 301, 1 }, (int[]) { 450, 300 });
  gkyl_range_init(&conn_2[5].range, 2, (int[]) { 451, 1 }, (int[]) { 600, 300 });
  
  struct gkyl_comm_conn *block_conn[] = { conn_0, conn_1, conn_2 };

  for (int bid=0; bid<num_blocks; ++bid) {
    int start_ns= 0;
    for (int brank=0; brank<num_cuts[bid]; ++brank) {
      struct gkyl_multib_comm_conn *mbcc = gkyl_multib_comm_conn_new_recv_from_connections(bid, brank, nghost, nconnected[bid], block_list[bid], dir, decomp);
      TEST_CHECK( num_recv_neigh[bid] == mbcc->num_comm_conn );
      for (int ns=0; ns<mbcc->num_comm_conn; ++ns) {
        // need to get the actual rank that owns this cut
        int rank_list[gkyl_rrobin_decomp_nranks(round_robin_decomp, mbcc->comm_conn[ns].block_id)];
        int rank_idx = mbcc->comm_conn[ns].rank;
        gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc->comm_conn[ns].block_id, rank_list);
        mbcc->comm_conn[ns].rank = rank_list[rank_idx];
        TEST_CHECK( block_conn[bid][start_ns+ns].block_id == mbcc->comm_conn[ns].block_id);
        TEST_CHECK( block_conn[bid][start_ns+ns].rank == mbcc->comm_conn[ns].rank);
        TEST_CHECK( gkyl_range_compare(&block_conn[bid][start_ns+ns].range, &mbcc->comm_conn[ns].range) );
      }

      start_ns+=mbcc->num_comm_conn;
      gkyl_multib_comm_conn_release(mbcc);
    }
  }

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(decomp);
  
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
  gkyl_comm_release(comm);
}

static void
test_L_domain_allgather_dir0_cuts2_par()
{
  printf("\n");
  // Create world comm.
  bool use_mpi = true;
  bool use_gpu = false;
  struct gkyl_comm* comm = comm_new(use_mpi, use_gpu, stderr);

  int my_rank, num_ranks;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &num_ranks);


  int num_blocks = 3; // L-shaped example.
  int ndim = 2;

  int ncuts_block2 = num_ranks > 1 ? 2 : 1;
  
  int cuts_flat[] = {
    1, 1, // Block 0.
    1, 1, // Block 1.
    ncuts_block2, 1, // Block 2.
  };
  int **cuts = cuts_array_new(num_blocks, ndim, cuts_flat);
  struct gkyl_block_geom *geom  = create_L_domain_block_geom(cuts);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  // Construct decomp objects.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    branks[i] = prod_of_elements_int(ndim, cuts[i]);
  }
  cuts_array_release(num_blocks, cuts);
  const struct gkyl_rrobin_decomp* round_robin_decomp = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // Allocate enough space.

  int num_cuts[num_blocks];
  int nghost[] = { 1, 1 };

  // Setup for a gather along x
  int dir = 0;
  int nconnected[num_blocks];
  int **block_list = gkyl_malloc(num_blocks*sizeof(int*));
  for (int bidx=0; bidx<num_blocks; ++bidx) {
    nconnected[bidx] = gkyl_multib_conn_get_num_connected(topo, bidx, dir, 0, GKYL_CONN_ALL);
    block_list[bidx] = gkyl_malloc(nconnected[bidx]*sizeof(int));
    gkyl_multib_conn_get_connection(topo, bidx, dir, 0, GKYL_CONN_ALL, block_list[bidx]);
  }


  //int block_list[3][2] = {{0},{1,2},{1,2}};
  //int nconnected[3] = {1,2,2};
  
  // construct decomp objects
  struct gkyl_rect_decomp **decomp =
    gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    num_cuts[i] = 1;
    for (int d=0; d<topo->ndim; ++d)
      num_cuts[i] *= ginfo->cuts[d];
    
    struct gkyl_range range;
    gkyl_create_global_range(2, ginfo->cells, &range);
    decomp[i] = gkyl_rect_decomp_new_from_cuts(2, ginfo->cuts, &range);
  }


  int lidx = 0;
  int local_blocks[num_blocks];
  int num_local_blocks = 0;
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(round_robin_decomp, i, rank_list);
    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);
    if (is_my_rank_in_decomp) {
      local_blocks[lidx++] = i;
      num_local_blocks += 1;      
    }
  }


  struct gkyl_array **array_local = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_array*));
  struct gkyl_array **array_global = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_array *));

  struct gkyl_range **local_ranges = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_range *));
  struct gkyl_range **local_ranges_ext = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_range *));
  struct gkyl_range **global_ranges = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_range *));
  struct gkyl_range **global_ranges_ext = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_range *));


  for (int bI= 0; bI<num_local_blocks; bI++) {
    int bid = local_blocks[bI];
    global_ranges[bI] = gkyl_malloc( sizeof(struct gkyl_range));
    global_ranges_ext[bI] = gkyl_malloc( sizeof(struct gkyl_range));
    local_ranges[bI] = gkyl_malloc( sizeof(struct gkyl_range));
    local_ranges_ext[bI] = gkyl_malloc( sizeof(struct gkyl_range));
    gkyl_multib_comm_conn_create_multib_ranges_in_dir(global_ranges_ext[bI], global_ranges[bI], nghost, nconnected[bid], block_list[bid], dir, decomp);

    gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;
    //local_ranges[bI] = &decomp[bid]->ranges[brank];
    gkyl_create_ranges(&decomp[bid]->ranges[brank], nghost, local_ranges_ext[bI], local_ranges[bI]);
  }

  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  // populate locals
  for (int bI=0; bI<num_local_blocks; ++bI) {
    int bid = local_blocks[bI];
    gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;
    array_local[bI] = mkarr(false, basis.num_basis, local_ranges_ext[bI]->volume);
    gkyl_array_shiftc(array_local[bI], sqrt(pow(2,ndim)), 0); // Sets es_energy_fac=1.
    if (num_ranks > 1) 
      gkyl_array_scale(array_local[bI], 0.5*my_rank);
    else
      gkyl_array_scale(array_local[bI], 10.0*bid);
    array_global[bI] = mkarr(false,  basis.num_basis, global_ranges_ext[bI]->volume);
  }


  struct gkyl_multib_comm_conn **mbcc_send = gkyl_malloc(num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  struct gkyl_multib_comm_conn **mbcc_recv = gkyl_malloc(num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  for (int bI= 0; bI<num_local_blocks; bI++) {
      int bid = local_blocks[bI];
      gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
      int brank = -1;
      for (int i=0; i<branks[bid]; ++i)
        if (rank_list[i] == my_rank) brank = i;

      mbcc_send[bI] = gkyl_multib_comm_conn_new_send_from_connections(bid, brank, nghost, nconnected[bid], block_list[bid], dir, decomp);
      mbcc_recv[bI] = gkyl_multib_comm_conn_new_recv_from_connections(bid, brank, nghost, nconnected[bid], block_list[bid], dir, decomp);

      for (int ns=0; ns<mbcc_send[bI]->num_comm_conn; ++ns) {
        // need to get the actual rank that owns this cut
        int rank_idx = mbcc_send[bI]->comm_conn[ns].rank;
        gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc_send[bI]->comm_conn[ns].block_id, rank_list);
        mbcc_send[bI]->comm_conn[ns].rank = rank_list[rank_idx];
        // Make range a subrange
        mbcc_send[bI]->comm_conn[ns].range = *local_ranges[bI];
      }
      for (int nr=0; nr<mbcc_recv[bI]->num_comm_conn; ++nr) {
        // need to get the actual rank that owns this cut
        int rank_idx = mbcc_recv[bI]->comm_conn[nr].rank;
        gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc_recv[bI]->comm_conn[nr].block_id, rank_list);
        mbcc_recv[bI]->comm_conn[nr].rank = rank_list[rank_idx];
        // Make range a subrange
        gkyl_sub_range_init(&mbcc_recv[bI]->comm_conn[nr].range, global_ranges_ext[bI], mbcc_recv[bI]->comm_conn[nr].range.lower, mbcc_recv[bI]->comm_conn[nr].range.upper);
      }

  }

  int stat = gkyl_multib_comm_conn_array_transfer(comm, num_local_blocks, local_blocks, mbcc_send, mbcc_recv, array_local, array_global);


  for (int bI=0; bI<num_local_blocks; ++bI) {
    struct gkyl_rect_grid grid;
    int cells[2];
    cells[0] = gkyl_range_shape(global_ranges[bI],0);
    cells[1] = gkyl_range_shape(global_ranges[bI],1);
    double gridlo[2] = {0.0,0.0};
    double l0 = ((double) cells[0])/300.0;
    double l1  = ((double) cells[1])/300.0;
    double gridup[2] = {l0,l1 };
    gkyl_rect_grid_init(&grid, 2, gridlo, gridup, cells);
    char str[50];
    sprintf(str, "lb%d_r%d.gkyl",bI, my_rank);
    gkyl_grid_sub_array_write(&grid, global_ranges[bI], 0, array_global[bI], str);
  }

  for (int bI = 0; bI<num_local_blocks; ++bI) {
    gkyl_multib_comm_conn_release(mbcc_send[bI]);
    gkyl_multib_comm_conn_release(mbcc_recv[bI]);
  }
  gkyl_free(mbcc_send);
  gkyl_free(mbcc_recv);
  gkyl_free(local_ranges);
  gkyl_free(array_local);
  gkyl_free(array_global);

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(decomp);
  
  for (int bidx=0; bidx<num_blocks; ++bidx) gkyl_free(block_list[bidx]);
  gkyl_free(block_list);
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
  gkyl_comm_release(comm);
}

static void
test_SOL_domain_allgather_dir1_cuts2_par(bool use_gpu)
{
  printf("\n");
  // Create world comm.
  bool use_mpi = true;
  struct gkyl_comm* comm = comm_new(use_mpi, use_gpu, stderr);
  printf("made comms\n");

  int my_rank, num_ranks;
  gkyl_comm_get_rank(comm, &my_rank);
  gkyl_comm_get_size(comm, &num_ranks);

  if (num_ranks!=2) 
    return;

  int num_blocks = 3; // SOL-shaped example.
  int ndim = 2;

  int ncuts_block1 = num_ranks > 1 ? 2 : 1;
  
  int cuts_flat[] = {
    1, 1, // Block 0.
    1, ncuts_block1, // Block 1.
    1, 1, // Block 2.
  };
  int **cuts = cuts_array_new(num_blocks, ndim, cuts_flat);
  struct gkyl_block_geom *geom  = create_SOL_domain_block_geom(cuts);
  struct gkyl_block_topo *topo = gkyl_block_geom_topo(geom);

  // Construct decomp objects.
  int *branks = gkyl_malloc(sizeof(int[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    branks[i] = prod_of_elements_int(ndim, cuts[i]);
  }
  cuts_array_release(num_blocks, cuts);
  const struct gkyl_rrobin_decomp* round_robin_decomp = gkyl_rrobin_decomp_new(num_ranks, num_blocks, branks);

  // Check rank info. This is what was used to produce the expected results
  //int num_ranks_pb[3] = {-1};
  //num_ranks_pb[0] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 0);
  //num_ranks_pb[1] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 1);
  //num_ranks_pb[2] = gkyl_rrobin_decomp_nranks(round_robin_decomp, 2);

  //int ranks_0[num_ranks_pb[0]];
  //int ranks_1[num_ranks_pb[1]];
  //int ranks_2[num_ranks_pb[2]];
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 0, ranks_0);
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 1, ranks_1);
  //gkyl_rrobin_decomp_getranks(round_robin_decomp, 2, ranks_2);

  //printf("num_ranks_pb = [%d, %d, %d]\n", num_ranks_pb[0], num_ranks_pb[1], num_ranks_pb[2]);
  //printf("Block 0 has ranks: ");
  //for (int i=0; i<num_ranks_pb[0]; i++) printf("%d ", ranks_0[i]);
  //printf("\n");
  //printf("Block 1 has ranks: ");
  //for (int i=0; i<num_ranks_pb[1]; i++) printf("%d ", ranks_1[i]);
  //printf("\n");
  //printf("Block 2 has ranks: ");
  //for (int i=0; i<num_ranks_pb[2]; i++) printf("%d ", ranks_2[i]);
  //printf("\n");

  int *rank_list = gkyl_malloc(sizeof(int[num_ranks])); // Allocate enough space.

  int num_cuts[num_blocks];
  int nghost[] = { 1, 1 };

  // Setup for a gather along x
  int dir = 1;
  int nconnected[num_blocks];
  int **block_list = gkyl_malloc(num_blocks*sizeof(int*));
  for (int bidx=0; bidx<num_blocks; ++bidx) {
    nconnected[bidx] = gkyl_multib_conn_get_num_connected(topo, bidx, dir, 0, GKYL_CONN_ALL);
    block_list[bidx] = gkyl_malloc(nconnected[bidx]*sizeof(int));
    gkyl_multib_conn_get_connection(topo, bidx, dir, 0, GKYL_CONN_ALL, block_list[bidx]);
  }


  
  // construct decomp objects
  struct gkyl_rect_decomp **decomp =
    gkyl_malloc(sizeof(struct gkyl_rect_decomp*[num_blocks]));
  for (int i=0; i<num_blocks; ++i) {
    const struct gkyl_block_geom_info *ginfo = gkyl_block_geom_get_block(geom, i);

    num_cuts[i] = 1;
    for (int d=0; d<topo->ndim; ++d)
      num_cuts[i] *= ginfo->cuts[d];
    
    struct gkyl_range range;
    gkyl_create_global_range(2, ginfo->cells, &range);
    decomp[i] = gkyl_rect_decomp_new_from_cuts(2, ginfo->cuts, &range);
  }


  int lidx = 0;
  int local_blocks[num_blocks];
  int num_local_blocks = 0;
  for (int i=0; i<num_blocks; ++i) {
    gkyl_rrobin_decomp_getranks(round_robin_decomp, i, rank_list);
    bool is_my_rank_in_decomp = has_int(branks[i], my_rank, rank_list);
    if (is_my_rank_in_decomp) {
      local_blocks[lidx++] = i;
      num_local_blocks += 1;      
    }
  }


  struct gkyl_array **array_local = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_array*));
  struct gkyl_array **array_global = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_array *));
  struct gkyl_array **array_local_ho = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_array*));
  struct gkyl_array **array_global_ho = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_array *));

  struct gkyl_range **local_ranges = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_range *));
  struct gkyl_range **local_ranges_ext = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_range *));
  struct gkyl_range **global_ranges = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_range *));
  struct gkyl_range **global_ranges_ext = gkyl_malloc(num_local_blocks* sizeof(struct gkyl_range *));


  for (int bI= 0; bI<num_local_blocks; bI++) {
    int bid = local_blocks[bI];
    global_ranges[bI] = gkyl_malloc( sizeof(struct gkyl_range));
    global_ranges_ext[bI] = gkyl_malloc( sizeof(struct gkyl_range));
    local_ranges[bI] = gkyl_malloc( sizeof(struct gkyl_range));
    local_ranges_ext[bI] = gkyl_malloc( sizeof(struct gkyl_range));
    gkyl_multib_comm_conn_create_multib_ranges_in_dir(global_ranges_ext[bI], global_ranges[bI], nghost, nconnected[bid], block_list[bid], dir, decomp);

    gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;
    gkyl_create_ranges(&decomp[bid]->ranges[brank], nghost, local_ranges_ext[bI], local_ranges[bI]);
  }

  printf("making arrays\n");
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  // populate locals
  for (int bI=0; bI<num_local_blocks; ++bI) {
    int bid = local_blocks[bI];
    gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;
    array_local[bI] = mkarr(use_gpu, basis.num_basis, local_ranges_ext[bI]->volume);
    array_local_ho[bI] = use_gpu ? mkarr(false, basis.num_basis, local_ranges_ext[bI]->volume)
	                         : gkyl_array_acquire(array_local[bI]);
    gkyl_array_shiftc(array_local[bI], sqrt(pow(2,ndim)), 0); // Sets es_energy_fac=1.
    if (num_ranks > 1) 
      gkyl_array_scale(array_local[bI], 0.5*my_rank);
    else
      gkyl_array_scale_range(array_local[bI], 100.0*bid, local_ranges[bI]);
    array_global[bI] = mkarr(use_gpu,  basis.num_basis, global_ranges_ext[bI]->volume);
    array_global_ho[bI] = use_gpu ? mkarr(false,  basis.num_basis, global_ranges_ext[bI]->volume) : gkyl_array_acquire(array_global[bI]);
  }
  printf("made arrays\n");


  struct gkyl_multib_comm_conn **mbcc_send = gkyl_malloc(num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  struct gkyl_multib_comm_conn **mbcc_recv = gkyl_malloc(num_local_blocks * sizeof(struct gkyl_multib_comm_conn *));
  for (int bI= 0; bI<num_local_blocks; bI++) {
    int bid = local_blocks[bI];
    gkyl_rrobin_decomp_getranks(round_robin_decomp, bid, rank_list);
    int brank = -1;
    for (int i=0; i<branks[bid]; ++i)
      if (rank_list[i] == my_rank) brank = i;

    mbcc_send[bI] = gkyl_multib_comm_conn_new_send_from_connections(bid, brank, nghost, nconnected[bid], block_list[bid], dir, decomp);
    mbcc_recv[bI] = gkyl_multib_comm_conn_new_recv_from_connections(bid, brank, nghost, nconnected[bid], block_list[bid], dir, decomp);

    for (int ns=0; ns<mbcc_send[bI]->num_comm_conn; ++ns) {
      // need to get the actual rank that owns this cut
      int rank_idx = mbcc_send[bI]->comm_conn[ns].rank;
      gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc_send[bI]->comm_conn[ns].block_id, rank_list);
      mbcc_send[bI]->comm_conn[ns].rank = rank_list[rank_idx];
      // Make range a subrange
      mbcc_send[bI]->comm_conn[ns].range = *local_ranges[bI];
    }
    for (int nr=0; nr<mbcc_recv[bI]->num_comm_conn; ++nr) {
      // need to get the actual rank that owns this cut
      int rank_idx = mbcc_recv[bI]->comm_conn[nr].rank;
      gkyl_rrobin_decomp_getranks(round_robin_decomp, mbcc_recv[bI]->comm_conn[nr].block_id, rank_list);
      mbcc_recv[bI]->comm_conn[nr].rank = rank_list[rank_idx];
      // Make range a subrange
      gkyl_sub_range_init(&mbcc_recv[bI]->comm_conn[nr].range, global_ranges_ext[bI], mbcc_recv[bI]->comm_conn[nr].range.lower, mbcc_recv[bI]->comm_conn[nr].range.upper);
    }

    // Sort connections according to rank and block ID.
    gkyl_multib_comm_conn_sort(mbcc_recv[bI]);
    gkyl_multib_comm_conn_sort(mbcc_send[bI]);

  }


  printf("calling transfer\n");
  int stat = gkyl_multib_comm_conn_array_transfer(comm, num_local_blocks, local_blocks, mbcc_send, mbcc_recv, array_local, array_global);
  printf("did transfer\n");


  for (int bI=0; bI<num_local_blocks; ++bI) {
    gkyl_array_copy(array_global_ho[bI], array_global[bI]);
  }

  for (int bI=0; bI<num_local_blocks; ++bI) {
    struct gkyl_rect_grid grid;
    int cells[2];
    cells[0] = gkyl_range_shape(global_ranges[bI],0);
    cells[1] = gkyl_range_shape(global_ranges[bI],1);
    double gridlo[2] = {0.0,0.0};
    double l0 = ((double) cells[0])/4.0;
    double l1  = ((double) cells[1])/8.0;
    double gridup[2] = {l0,l1 };
    gkyl_rect_grid_init(&grid, 2, gridlo, gridup, cells);
    char str[50];
    sprintf(str, "lb%d_r%d.gkyl",bI, my_rank);
    gkyl_grid_sub_array_write(&grid, global_ranges[bI], 0, array_global_ho[bI], str);
  }

  printf("checking\n");
  for (int bI=0; bI<num_local_blocks; ++bI) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, global_ranges[bI]);
    double fcheck[4] = {0.0};
    while(gkyl_range_iter_next(&iter)) {
      if(iter.idx[1]<= 8) {
        fcheck[0] = 0.0;
        fcheck[1] = 0.0;
        fcheck[2] = 0.0;
        fcheck[3] = 0.0;
      }
      else if(iter.idx[1] <= 12) {
        fcheck[0] = 1.0;
        fcheck[1] = 0.0;
        fcheck[2] = 0.0;
        fcheck[3] = 0.0;
      }
      else if(iter.idx[1]<=16) {
        fcheck[0] = 0.0;
        fcheck[1] = 0.0;
        fcheck[2] = 0.0;
        fcheck[3] = 0.0;
      }
      else if(iter.idx[1] <=24) {
          fcheck[0] = 1.0;
          fcheck[1] = 0.0;
          fcheck[2] = 0.0;
          fcheck[3] = 0.0;
      }
      long loc = gkyl_range_idx(global_ranges[bI], iter.idx);
      const double *f = gkyl_array_cfetch(array_global_ho[bI], loc);
      for (int i=0; i<basis.num_basis; ++i) {
        TEST_CHECK( f[i] == fcheck[i] );
      }
    }
      
  }

  for (int bI = 0; bI<num_local_blocks; ++bI) {
    gkyl_multib_comm_conn_release(mbcc_send[bI]);
    gkyl_multib_comm_conn_release(mbcc_recv[bI]);
    gkyl_array_release(array_local_ho[bI]);
    gkyl_array_release(array_local[bI]);
    gkyl_array_release(array_global_ho[bI]);
    gkyl_array_release(array_global[bI]);
  }
  gkyl_free(mbcc_send);
  gkyl_free(mbcc_recv);
  gkyl_free(local_ranges);
  gkyl_free(array_local);
  gkyl_free(array_global);

  for (int i=0; i<num_blocks; ++i)
    gkyl_rect_decomp_release(decomp[i]);
  gkyl_free(decomp);
  
  for (int bidx=0; bidx<num_blocks; ++bidx) gkyl_free(block_list[bidx]);
  gkyl_free(block_list);
  gkyl_block_topo_release(topo);
  gkyl_block_geom_release(geom);
  gkyl_comm_release(comm);
}

static void
test_SOL_domain_allgather_dir1_cuts2_par_ho(void)
{
  test_SOL_domain_allgather_dir1_cuts2_par(false);
}

#ifdef GKYL_HAVE_NCCL
static void
test_SOL_domain_allgather_dir1_cuts2_par_dev(void)
{
  test_SOL_domain_allgather_dir1_cuts2_par(true);
}
#endif



TEST_LIST = {
  //{ "test_L_domain_send_connections_dir0_cuts1", test_L_domain_send_connections_dir0_cuts1},
  //{ "test_L_domain_recv_connections_dir0_cuts1", test_L_domain_recv_connections_dir0_cuts1},
  //{ "test_L_domain_send_connections_dir0_cuts2", test_L_domain_send_connections_dir0_cuts2},
  //
  //{ "test_L_domain_send_connections_dir0_cuts2_par", test_L_domain_send_connections_dir0_cuts2_par},
  //{ "test_L_domain_recv_connections_dir0_cuts2_par", test_L_domain_recv_connections_dir0_cuts2_par},
  //{ "test_L_domain_allgather_dir0_cuts2_par", test_L_domain_allgather_dir0_cuts2_par},
  { "test_SOL_domain_allgather_dir1_cuts2_par_ho", test_SOL_domain_allgather_dir1_cuts2_par_ho},
#ifdef GKYL_HAVE_NCCL
  { "test_SOL_domain_allgather_dir1_cuts2_par_dev", test_SOL_domain_allgather_dir1_cuts2_par_dev},
#endif
  { NULL, NULL },
};
