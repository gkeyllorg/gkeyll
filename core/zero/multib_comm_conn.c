#include <gkyl_alloc.h>
#include <gkyl_multib_comm_conn.h>
#include <gkyl_util.h>
#include <gkyl_multib_comm_conn_priv.h>

#include <string.h>
#include <assert.h>

struct multib_comm_conn {
  struct gkyl_multib_comm_conn mcc;
};

static void
multib_comm_conn_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_multib_comm_conn *mcc = container_of(ref, struct gkyl_multib_comm_conn, ref_count);
  struct multib_comm_conn *cconn = container_of(mcc, struct multib_comm_conn, mcc);

  if (cconn->mcc.num_comm_conn > 0)
    gkyl_free(cconn->mcc.comm_conn);
  gkyl_free(cconn);
}

struct gkyl_multib_comm_conn *
gkyl_multib_comm_conn_new(int num, const struct gkyl_comm_conn *comm_conn)
{
  struct multib_comm_conn *cconn = gkyl_malloc(sizeof *cconn);
  cconn->mcc.num_comm_conn = num;

  cconn->mcc.comm_conn = 0;
  if (num > 0)
    cconn->mcc.comm_conn = gkyl_malloc(sizeof(struct gkyl_comm_conn[num]));

  for (int i=0; i<num; ++i)
    memcpy(&cconn->mcc.comm_conn[i], &comm_conn[i], sizeof(struct gkyl_comm_conn));

  cconn->mcc.ref_count = gkyl_ref_count_init(multib_comm_conn_free);

  return &cconn->mcc;
}


// private method to compute send/recv connections
static struct gkyl_multib_comm_conn *
multib_comm_conn_new_sr(enum multib_send_recv sr,
  int block_id, int block_rank, const int *nghost,
  const struct gkyl_block_connections *block_conn, struct gkyl_rect_decomp **decomp)
{
  int ndim = decomp[0]->ndim;
  // determine maximum number of ranks we can send/recv data 
  int max_sr_ranks = 0;
  for (int d=0; d<ndim; ++d) {
    for (int e=0; e<2; ++e) {
    
      if (block_conn->connections[d][e].edge != GKYL_PHYSICAL) {
        int nranks = decomp[block_conn->connections[d][e].bid]->ndecomp;
        max_sr_ranks += nranks;
      }
    }
  }

  int comm_conn_idx = 0;
  struct gkyl_comm_conn *comm_conn
    = gkyl_malloc(sizeof(struct gkyl_comm_conn[max_sr_ranks]));

  const struct gkyl_range *src_pr = &decomp[block_id]->parent_range;
  const struct gkyl_range *src_br = &decomp[block_id]->ranges[block_rank];

  for (int dir=0; dir<ndim; ++dir) {
    for (int e=0; e<2; ++e) {
      int tar_bid = block_conn->connections[dir][e].bid;
      int tar_dir = block_conn->connections[dir][e].dir;
      enum gkyl_oriented_edge tar_edge = block_conn->connections[dir][e].edge;
      enum gkyl_oriented_edge src_edge = e == 0? GKYL_LOWER_POSITIVE : GKYL_UPPER_POSITIVE;

      if (block_conn->connections[dir][e].edge != GKYL_PHYSICAL) {

        const struct gkyl_rect_decomp *tar_decomp = decomp[tar_bid];
        // to find intersection with the neighbor we first need to get
        // the neigbor block into the index-space to source block: we
        // do this by reseting the lower indices of the neigbor's
        // parent
        const struct gkyl_range *tar_pr = &tar_decomp->parent_range;

        int new_lower[GKYL_MAX_DIM];
        for (int d=0; d<ndim; ++d)
          new_lower[d] = src_pr->lower[d];

        if (GKYL_UPPER_EDGE == e)
          new_lower[dir] += gkyl_range_shape(src_pr, dir);
        else
          new_lower[dir] -= gkyl_range_shape(tar_pr, dir);
        
        struct gkyl_range reset_tar_pr;
        gkyl_range_reset_lower(&reset_tar_pr, tar_pr, new_lower);

        int delta[GKYL_MAX_DIM] = { 0 };
        for (int d=0; d<ndim; ++d)
          delta[d] = reset_tar_pr.lower[d] - tar_pr->lower[d];

        int minus_delta[GKYL_MAX_DIM] = { 0 };
        for (int d=0; d<ndim; ++d)
          minus_delta[d] = -delta[d];

        struct gkyl_range_dir_edge dir_edge = gkyl_range_edge_match(&reset_tar_pr, src_br);
        if (dir_edge.eloc != GKYL_NO_EDGE) {
          // source block range touches the parent of the neigbor
          // block: find the appropriate sub-blocks it intersects

          for (int nn=0; nn<tar_decomp->ndecomp; ++nn) {
            struct gkyl_range sub_range;
            gkyl_range_shift(&sub_range, &tar_decomp->ranges[nn], delta);

            int elo[GKYL_MAX_DIM] = { 0 }, eup[GKYL_MAX_DIM] = { 0 };
            elo[dir] = eup[dir] = nghost[dir];  // only extend in 1 direction

            struct gkyl_range range_ext;
            if (GKYL_COMM_CONN_SEND == sr)
              gkyl_range_extend(&range_ext, &sub_range, elo, eup);
            else
              gkyl_range_extend(&range_ext, src_br, elo, eup);

            int is_inter;
            struct gkyl_range irng;
            if (GKYL_COMM_CONN_SEND == sr)
              is_inter = gkyl_range_intersect(&irng, &range_ext, src_br);
            else
              is_inter = gkyl_range_intersect(&irng, &range_ext, &sub_range);
            
            if (is_inter) {
              comm_conn[comm_conn_idx].sr = sr;
              comm_conn[comm_conn_idx].rank = nn;
              comm_conn[comm_conn_idx].block_id = tar_bid;
              comm_conn[comm_conn_idx].src_edge = src_edge;
              comm_conn[comm_conn_idx].tar_edge = tar_edge;
              memcpy(&comm_conn[comm_conn_idx].range, &irng, sizeof(struct gkyl_range));

              comm_conn_idx += 1;
            }
          }

        }
      }
    }
  }

  struct gkyl_multib_comm_conn *mbcc =
    gkyl_multib_comm_conn_new(comm_conn_idx, comm_conn);

  gkyl_free(comm_conn);
  
  return mbcc;
}

void
gkyl_multib_comm_conn_create_multib_ranges_in_dir(struct gkyl_range *multib_range_ext, struct gkyl_range *multib_range,
    const int *nghost, int nconnected, int* block_list, int dir, struct gkyl_rect_decomp **decomp)
{
  // Construct the multib range which spans all the parent ranges
  // Block list is always in order in dir
  // Block parent ranges always share range extents in the other directions
  int ndim = decomp[0]->ndim;
  int multib_lower[GKYL_MAX_DIM];
  int multib_upper[GKYL_MAX_DIM];
  for (int i =0; i<ndim; i++) {
    multib_lower[i] = decomp[block_list[0]]->parent_range.lower[i];
    multib_upper[i] = decomp[block_list[0]]->parent_range.upper[i];
  }
  for (int i=1; i<nconnected; i++) {
    multib_upper[dir] += gkyl_range_shape(&decomp[block_list[i]]->parent_range, dir);
  }
  struct gkyl_range multib_rng;
  gkyl_range_init(&multib_rng, ndim, multib_lower, multib_upper);
  gkyl_create_ranges(&multib_rng, nghost, multib_range_ext, multib_range);
}

// public method to compute send connections from block list
struct gkyl_multib_comm_conn *
gkyl_multib_comm_conn_new_send_from_connections(
  int block_id, int block_rank, const int *nghost,
  int nconnected, int* block_list, int dir,
  struct gkyl_rect_decomp **decomp)
{
  int ndim = decomp[0]->ndim;
  // determine maximum number of ranks we can send/recv data 
  int max_sr_ranks = 0;
  for (int i=0; i<nconnected; i++) {
    int nranks = decomp[block_list[i]]->ndecomp;
    max_sr_ranks += nranks;
  }

  int comm_conn_idx = 0;
  struct gkyl_comm_conn *comm_conn
    = gkyl_malloc(sizeof(struct gkyl_comm_conn[max_sr_ranks]));

  const struct gkyl_range *src_parent_range = &decomp[block_id]->parent_range;
  const struct gkyl_range *src_block_range = &decomp[block_id]->ranges[block_rank];

  struct gkyl_range cross_range, cross_range_ext;
  gkyl_multib_comm_conn_create_multib_ranges_in_dir(&cross_range_ext, &cross_range, nghost, nconnected, block_list, dir, decomp);

  // Get block indices into the index space of the cross range
  int new_lower[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d)
    new_lower[d] = src_parent_range->lower[d];

  int source_idx = -1;
  for(int i=0; i<nconnected; i++) { if (block_list[i] == block_id) source_idx = i; }
  for (int i=0; i<source_idx; i++) {
    new_lower[dir] += gkyl_range_shape(&decomp[block_list[i]]->parent_range, dir);
  }
  
  struct gkyl_range reset_src_parent_range;
  gkyl_range_reset_lower(&reset_src_parent_range, src_parent_range, new_lower);

  int delta[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<ndim; ++d)
    delta[d] = reset_src_parent_range.lower[d] - src_parent_range->lower[d];

  int minus_delta[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<ndim; ++d)
    minus_delta[d] = -delta[d];

  // The actual intersection should be between a our rank's range
  // and the cross range. So we will need to shift src_block_range into the cross range
  struct gkyl_range sub_range;
  gkyl_range_shift(&sub_range, src_block_range, delta);

  // Now loop over all the connected ranks and create a connection
  for (int ib=0; ib<nconnected; ib++) {
    int tar_ndecomp = decomp[block_list[ib]]->ndecomp;
    for (int ir=0; ir<tar_ndecomp; ir++) {
      int is_inter;
      struct gkyl_range irng;
      is_inter = gkyl_range_intersect(&irng, &cross_range, &sub_range);
      if (is_inter) {
        comm_conn[comm_conn_idx].sr = GKYL_COMM_CONN_SEND;
        comm_conn[comm_conn_idx].src_edge = 0;
        comm_conn[comm_conn_idx].tar_edge = 0;
        comm_conn[comm_conn_idx].rank = ir;
        comm_conn[comm_conn_idx].block_id = block_list[ib];
        memcpy(&comm_conn[comm_conn_idx].range, &irng, sizeof(struct gkyl_range));
        comm_conn_idx += 1;
      }
    }
  }
  struct gkyl_multib_comm_conn *mbcc =
    gkyl_multib_comm_conn_new(comm_conn_idx, comm_conn);

  gkyl_free(comm_conn);
  
  return mbcc;
}

// public method to compute recv connections from block list
struct gkyl_multib_comm_conn *
gkyl_multib_comm_conn_new_recv_from_connections(
  int block_id, int block_rank, const int *nghost,
  int nconnected, int* block_list, int dir,
  struct gkyl_rect_decomp **decomp)
{
  int ndim = decomp[0]->ndim;
  // determine maximum number of ranks we can send/recv data 
  int max_sr_ranks = 0;
  for (int i=0; i<nconnected; i++) {
    int nranks = decomp[block_list[i]]->ndecomp;
    max_sr_ranks += nranks;
  }

  int comm_conn_idx = 0;
  struct gkyl_comm_conn *comm_conn
    = gkyl_malloc(sizeof(struct gkyl_comm_conn[max_sr_ranks]));

  const struct gkyl_range *src_parent_range = &decomp[block_id]->parent_range;
  const struct gkyl_range *src_block_range = &decomp[block_id]->ranges[block_rank];

  struct gkyl_range cross_range, cross_range_ext;
  gkyl_multib_comm_conn_create_multib_ranges_in_dir(&cross_range_ext, &cross_range, nghost, nconnected, block_list, dir, decomp);

  // Need to get other block indices into the index space of the cross range
  for (int ib=0; ib<nconnected; ib++) {
    int tar_bid = block_list[ib];
    const struct gkyl_rect_decomp *tar_decomp = decomp[tar_bid];
    const struct gkyl_range *tar_parent_range = &tar_decomp->parent_range;
    int new_lower[GKYL_MAX_DIM];
    for (int d=0; d<ndim; ++d)
      new_lower[d] = tar_parent_range->lower[d];

    int tar_idx = -1;
    for(int i=0; i<nconnected; i++) { if (block_list[i] == tar_bid) tar_idx = i; }
    for (int i=0; i<tar_idx; i++) {
      new_lower[dir] += gkyl_range_shape(&decomp[block_list[i]]->parent_range, dir);
    }
    
    struct gkyl_range reset_tar_parent_range;
    gkyl_range_reset_lower(&reset_tar_parent_range, tar_parent_range, new_lower);

    int delta[GKYL_MAX_DIM] = { 0 };
    for (int d=0; d<ndim; ++d)
      delta[d] = reset_tar_parent_range.lower[d] - tar_parent_range->lower[d];

    int minus_delta[GKYL_MAX_DIM] = { 0 };
    for (int d=0; d<ndim; ++d)
      minus_delta[d] = -delta[d];

  // The actual intersection should be between the target rank's range
  // and the cross range. So we will need to shift tar_block_range into the cross range
    int tar_ndecomp = tar_decomp->ndecomp;
    for (int ir=0; ir<tar_ndecomp; ir++) {
      const struct gkyl_range *tar_block_range = &decomp[tar_bid]->ranges[ir];
      struct gkyl_range sub_range;
      gkyl_range_shift(&sub_range, tar_block_range, delta);
      int is_inter;
      struct gkyl_range irng;
      is_inter = gkyl_range_intersect(&irng, &cross_range, &sub_range);
      if (is_inter) {
        comm_conn[comm_conn_idx].sr = GKYL_COMM_CONN_RECV;
        comm_conn[comm_conn_idx].src_edge = 0;
        comm_conn[comm_conn_idx].tar_edge = 0;
        comm_conn[comm_conn_idx].rank = ir;
        comm_conn[comm_conn_idx].block_id = tar_bid;
        memcpy(&comm_conn[comm_conn_idx].range, &irng, sizeof(struct gkyl_range));
        comm_conn_idx += 1;
      }
    }
  }
  struct gkyl_multib_comm_conn *mbcc =
    gkyl_multib_comm_conn_new(comm_conn_idx, comm_conn);

  gkyl_free(comm_conn);
  
  return mbcc;
}

struct gkyl_multib_comm_conn *
gkyl_multib_comm_conn_new_send(
  int block_id, int block_rank, const int *nghost,
  const struct gkyl_block_connections *block_conn, struct gkyl_rect_decomp **decomp)
{
  return multib_comm_conn_new_sr(GKYL_COMM_CONN_SEND, block_id, block_rank, nghost, block_conn, decomp);
}

struct gkyl_multib_comm_conn *
gkyl_multib_comm_conn_new_recv(
  int block_id, int block_rank, const int *nghost,
  const struct gkyl_block_connections *block_conn, struct gkyl_rect_decomp **decomp)
{
  return multib_comm_conn_new_sr(GKYL_COMM_CONN_RECV, block_id, block_rank, nghost, block_conn, decomp);
}

int
gkyl_multib_comm_conn_array_transfer(struct gkyl_comm *comm, int num_blocks_local, const int *blocks_local,
  struct gkyl_multib_comm_conn **mbcc_send, struct gkyl_multib_comm_conn **mbcc_recv,
  struct gkyl_array **arr_send, struct gkyl_array **arr_recv)
{
  int err;
  if (strcmp(comm->id, "null_comm") == 0) {
    err = gkyl_multib_comm_conn_array_transfer_null(comm, num_blocks_local, blocks_local,
      mbcc_send, mbcc_recv, arr_send, arr_recv);
  }
  else if (strcmp(comm->id, "mpi_comm") == 0) {
    err = gkyl_multib_comm_conn_array_transfer_mpi(comm, num_blocks_local, blocks_local,
      mbcc_send, mbcc_recv, arr_send, arr_recv);
  }
  else if (strcmp(comm->id, "nccl_comm") == 0) {
    err = gkyl_multib_comm_conn_array_transfer_nccl(comm, num_blocks_local, blocks_local,
      mbcc_send, mbcc_recv, arr_send, arr_recv);
  }
  else
    assert(false);

  return err;
}

static void
swap_comm_conns(struct gkyl_comm_conn *ccj, struct gkyl_comm_conn *cck)
{
  struct gkyl_comm_conn cc_tmp = *ccj;
  *ccj = *cck;
  *cck = cc_tmp;
}

void
gkyl_multib_comm_conn_sort(struct gkyl_multib_comm_conn *mbcc)
{
  int num_conn = mbcc->num_comm_conn;
  // First sort connections in ascending rank (w/ bubble sort).
  for (int i=0; i<num_conn-1; i++) {
    bool swapped = false;
    for (int j=0; j<num_conn-i-1; j++) {
      struct gkyl_comm_conn *cj = &mbcc->comm_conn[j], *cjp1 = &mbcc->comm_conn[j+1];
      if (cj->rank > cjp1->rank) {
        swap_comm_conns(cj, cjp1);
        swapped = true;
      }
    }
    if (swapped == false)
      break; // Stop if no swaps happened.
  }

  // Now sort connections in ascending block ID.
  for (int i=0; i<num_conn-1; i++) {
    bool swapped = false;
    for (int j=0; j<num_conn-i-1; j++) {
      struct gkyl_comm_conn *cj = &mbcc->comm_conn[j], *cjp1 = &mbcc->comm_conn[j+1];
      if ((cj->rank == cjp1->rank) && (cj->block_id > cjp1->block_id)) {
        swap_comm_conns(cj, cjp1);
        swapped = true;
      }
    }
    if (swapped == false)
      break; // Stop if no swaps happened.
  }

  // Now sort connections in ascending tar/src edge if receiver/sender
  for (int i=0; i<num_conn-1; i++) {
    bool swapped = false;
    for (int j=0; j<num_conn-i-1; j++) {
      struct gkyl_comm_conn *cj = &mbcc->comm_conn[j], *cjp1 = &mbcc->comm_conn[j+1];
      int edge = cj->sr == GKYL_COMM_CONN_SEND ? cj->src_edge : cj->tar_edge;
      int edgep1 = cjp1->sr == GKYL_COMM_CONN_SEND ? cjp1->src_edge : cjp1->tar_edge;
      if ((cj->rank == cjp1->rank) && (cj->block_id == cjp1->block_id) && (edge > edgep1) ) {
        swap_comm_conns(cj, cjp1);
        swapped = true;
      }
    }
    if (swapped == false)
      break; // Stop if no swaps happened.
  }

}

void
gkyl_multib_comm_conn_release(const struct gkyl_multib_comm_conn *cconn)
{
  if (cconn)
    gkyl_ref_count_dec(&cconn->ref_count);
}
