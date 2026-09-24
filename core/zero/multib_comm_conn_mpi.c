#include <gkyl_multib_comm_conn_priv.h>
#include <gkyl_comm_priv.h>

#ifdef GKYL_HAVE_MPI

#include <gkyl_mpi_comm_priv.h>

int
gkyl_multib_comm_conn_array_transfer_mpi(struct gkyl_comm *comm, int num_blocks_local, const int *local_blocks,
  struct gkyl_multib_comm_conn **mbcc_send, struct gkyl_multib_comm_conn **mbcc_recv,
  struct gkyl_array **arr_send, struct gkyl_array **arr_recv)
{
  struct mpi_comm *mpi = container_of(comm, struct mpi_comm, priv_comm.pub_comm);

  int my_rank;
  gkyl_comm_get_rank(comm, &my_rank);

  int tag = MPI_BASE_TAG;

  // post nonblocking recv to get data into ghost-cells  
  int nridx = 0;
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct gkyl_multib_comm_conn *mbcc_r = mbcc_recv[bI];

    for (int n=0; n<mbcc_r->num_comm_conn; ++n) {
      int nid = mbcc_r->comm_conn[n].rank;
      int bid = mbcc_r->comm_conn[n].block_id;
      int e = mbcc_r->comm_conn[n].tar_edge;
      
      size_t recv_vol = arr_recv[bI]->esznc*mbcc_r->comm_conn[n].range.volume;

      if (recv_vol>0) {
        if (gkyl_mem_buff_size(mpi->recv[nridx].buff) < recv_vol)
          gkyl_mem_buff_resize(mpi->recv[nridx].buff, recv_vol);

        int rtag = tag + 1000*nid + 100*e + bid;

        MPI_Irecv(gkyl_mem_buff_data(mpi->recv[nridx].buff),
          recv_vol, MPI_CHAR, nid, rtag, mpi->mcomm, &mpi->recv[nridx].status);

        nridx += 1;
      }
    }
  }

  // post non-blocking sends of skin-cell data to neighbors
  int nsidx = 0;
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct gkyl_multib_comm_conn *mbcc_s = mbcc_send[bI];

    for (int n=0; n<mbcc_s->num_comm_conn; ++n) {
      int nid = mbcc_s->comm_conn[n].rank;
      int e = mbcc_s->comm_conn[n].src_edge;
    
      size_t send_vol = arr_send[bI]->esznc*mbcc_s->comm_conn[n].range.volume;

      if (send_vol>0) {
        if (gkyl_mem_buff_size(mpi->send[nsidx].buff) < send_vol)
          gkyl_mem_buff_resize(mpi->send[nsidx].buff, send_vol);
        
        gkyl_array_copy_to_buffer(gkyl_mem_buff_data(mpi->send[nsidx].buff),
          arr_send[bI], &mbcc_s->comm_conn[n].range);

        int stag = tag + 1000*my_rank + 100*e + local_blocks[bI];

        MPI_Isend(gkyl_mem_buff_data(mpi->send[nsidx].buff),
          send_vol, MPI_CHAR, nid, stag, mpi->mcomm, &mpi->send[nsidx].status);

        nsidx += 1;
      }
    }
  }

  // complete send
  nsidx = 0;
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct gkyl_multib_comm_conn *mbcc_s = mbcc_send[bI];
    for (int n=0; n<mbcc_s->num_comm_conn; ++n) {
      int issend = mbcc_s->comm_conn[n].range.volume;
      if (issend) {
        MPI_Wait(&mpi->send[nsidx].status, MPI_STATUS_IGNORE);
        nsidx += 1;
      }
    }
  }

  // complete recv, copying data into ghost-cells
  nridx = 0;
  for (int bI=0; bI<num_blocks_local; ++bI) {
    struct gkyl_multib_comm_conn *mbcc_r = mbcc_recv[bI];
    for (int n=0; n<mbcc_r->num_comm_conn; ++n) {
      int isrecv = mbcc_r->comm_conn[n].range.volume;
      if (isrecv) {
        MPI_Wait(&mpi->recv[nridx].status, MPI_STATUS_IGNORE);
        
        gkyl_array_copy_from_buffer(arr_recv[bI],
          gkyl_mem_buff_data(mpi->recv[nridx].buff),
          &(mbcc_r->comm_conn[n].range)
        );

        nridx += 1;
      }
    }
  }
  
  return 0;
}

#else 

int
gkyl_multib_comm_conn_array_transfer_mpi(struct gkyl_comm *comm, int num_blocks_local, const int *local_blocks,
  struct gkyl_multib_comm_conn **mbcc_send, struct gkyl_multib_comm_conn **mbcc_recv,
  struct gkyl_array **arr_send, struct gkyl_array **arr_recv)
{
  return 1;
}

#endif
