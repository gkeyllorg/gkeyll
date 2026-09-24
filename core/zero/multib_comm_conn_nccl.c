#include <gkyl_multib_comm_conn_priv.h>
#include <gkyl_comm_priv.h>

#ifdef GKYL_HAVE_NCCL

#include <gkyl_nccl_comm_priv.h>

int
gkyl_multib_comm_conn_array_transfer_nccl(struct gkyl_comm *comm, int num_blocks_local, const int *local_blocks,
  struct gkyl_multib_comm_conn **mbcc_send, struct gkyl_multib_comm_conn **mbcc_recv,
  struct gkyl_array **arr_send, struct gkyl_array **arr_recv)
{
  struct nccl_comm *nccl = container_of(comm, struct nccl_comm, priv_comm.pub_comm);

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);

  checkNCCL(ncclGroupStart());

  // post nonblocking recv to get data into ghost-cells  
  int nridx = 0;
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct gkyl_multib_comm_conn *mbcc_r = mbcc_recv[bI];

    for (int nr=0; nr<mbcc_r->num_comm_conn; ++nr) {
      int nid = mbcc_r->comm_conn[nr].rank;

      if (nid != my_rank) {
        int bid = mbcc_r->comm_conn[nr].block_id;
        
        size_t recv_vol = arr_recv[bI]->esznc*mbcc_r->comm_conn[nr].range.volume;
  
        if (recv_vol>0) {
          if (gkyl_mem_buff_size(nccl->recv[nridx].buff) < recv_vol)
            gkyl_mem_buff_resize(nccl->recv[nridx].buff, recv_vol);
  
            checkNCCL(ncclRecv(gkyl_mem_buff_data(nccl->recv[nridx].buff),
              recv_vol, ncclChar, nid, nccl->ncomm, nccl->custream));
  
          nridx += 1;
        }
      }
    }
  }

  // post non-blocking sends of skin-cell data to neighbors
  int nsidx = 0;
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct gkyl_multib_comm_conn *mbcc_s = mbcc_send[bI];

    for (int ns=0; ns<mbcc_s->num_comm_conn; ++ns) {
      int nid = mbcc_s->comm_conn[ns].rank;

      if (nid != my_rank) {
        size_t send_vol = arr_send[bI]->esznc*mbcc_s->comm_conn[ns].range.volume;

        if (send_vol>0) {
          if (gkyl_mem_buff_size(nccl->send[nsidx].buff) < send_vol)
            gkyl_mem_buff_resize(nccl->send[nsidx].buff, send_vol);
          
          gkyl_array_copy_to_buffer(gkyl_mem_buff_data(nccl->send[nsidx].buff),
            arr_send[bI], &mbcc_s->comm_conn[ns].range);

          checkNCCL(ncclSend(gkyl_mem_buff_data(nccl->send[nsidx].buff),
            send_vol, ncclChar, nid, nccl->ncomm, nccl->custream));

          nsidx += 1;
        }
      }
    }
  }

  checkNCCL(ncclGroupEnd());

  // Simply do a copy for the connections whose rank is the same as my rank.
  // This is the same code in gkyl_multib_comm_conn_array_transfer_null.
  for (int bI=0; bI<num_blocks_local; ++bI) {
    int bid = local_blocks[bI];
    struct gkyl_multib_comm_conn *mbcc_r = mbcc_recv[bI];

    for (int nr=0; nr<mbcc_r->num_comm_conn; ++nr) {
      int nid = mbcc_r->comm_conn[nr].rank;

      if (nid == my_rank) {
        int bid_src = mbcc_r->comm_conn[nr].block_id;
        struct gkyl_range *range_dest = &mbcc_r->comm_conn[nr].range;
        int e = mbcc_r->comm_conn[nr].tar_edge;

        // The send connections may not be ordered in the same way as the recvs.
        // Loop over the send connections and find the one from the same block.
        int bid_src_idx = -1;
        for (int cI=0; cI<num_blocks_local; ++cI) {
          if (local_blocks[cI] == bid_src) {
            bid_src_idx = cI;
            break;
          }
        }
        struct gkyl_multib_comm_conn *mbcc_s = mbcc_send[bid_src_idx];
        int conn_src_idx = -1;
        for (int ns=0; ns<mbcc_s->num_comm_conn; ++ns) {
          if (mbcc_s->comm_conn[ns].block_id == bid && mbcc_s->comm_conn[ns].src_edge == e) {
            conn_src_idx = ns;
            break;
          }
        }

        struct gkyl_range *range_src = &mbcc_s->comm_conn[conn_src_idx].range;

        gkyl_array_copy_range_to_range(arr_recv[bI], arr_send[bid_src_idx], range_dest, range_src);
      }
    }
  }

  // Complete sends and recvs.
  ncclResult_t nstat;
  do {
    checkNCCL(ncclCommGetAsyncError(nccl->ncomm, &nstat));
  } while(nstat == ncclInProgress);
  checkCuda(cudaStreamSynchronize(nccl->custream));

  // Copy data into ghost-cells.
  nridx = 0;
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct gkyl_multib_comm_conn *mbcc_r = mbcc_recv[bI];
    for (int nr=0; nr<mbcc_r->num_comm_conn; ++nr) {
      int nid = mbcc_r->comm_conn[nr].rank;
      if (nid != my_rank) {
        int isrecv = mbcc_r->comm_conn[nr].range.volume;
        if (isrecv) {
          gkyl_array_copy_from_buffer(arr_recv[bI],
            gkyl_mem_buff_data(nccl->recv[nridx].buff),
            &(mbcc_r->comm_conn[nr].range)
          );
          nridx += 1;
        }
      }
    }
  }

  return 0;
}

#else

int
gkyl_multib_comm_conn_array_transfer_nccl(struct gkyl_comm *comm, int num_blocks_local, const int *local_blocks,
  struct gkyl_multib_comm_conn **mbcc_send, struct gkyl_multib_comm_conn **mbcc_recv,
  struct gkyl_array **arr_send, struct gkyl_array **arr_recv)
{
  return 1;
}

#endif
