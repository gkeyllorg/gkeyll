#pragma once

// Private header for mpi_comm. Do not include in user-facing header files.

#include <gkyl_mpi_comm.h>
#include <gkyl_alloc.h>

// Maximum number of recv neighbors: not sure hard-coding this is a
// good idea.
#define MAX_RECV_NEIGH 128

#define MPI_BASE_TAG 4242
#define MPI_BASE_PER_TAG 5252

// Object with a range, a status and a buffer used for send/recv.
struct comm_buff_stat {
  struct gkyl_range range;
  MPI_Request status;
  gkyl_mem_buff buff;
};

// Private struct wrapping MPI-specific code
struct mpi_comm {
  struct gkyl_comm_priv priv_comm; // base communicator

  MPI_Comm mcomm; // MPI communicator to use
  bool is_mcomm_allocated; // Is the mcomm allocated?
  struct gkyl_rect_decomp *decomp; // pre-computed decomposition
  long local_range_offset; // offset of the local region

  bool sync_corners; // should we sync corners?

  struct gkyl_rect_decomp_neigh *neigh; // neighbors of local region
  struct gkyl_rect_decomp_neigh *per_neigh[GKYL_MAX_DIM]; // periodic neighbors

  struct gkyl_range dir_edge; // for use in computing tags
  int is_on_edge[2][GKYL_MAX_DIM]; // flags to indicate if local range is on edge
  bool touches_any_edge; // true if this range touches any edge

  int nrecv; // number of elements in rinfo array
  struct comm_buff_stat recv[MAX_RECV_NEIGH]; // info for recv data

  int nsend; // number of elements in sinfo array
  struct comm_buff_stat send[MAX_RECV_NEIGH]; // info for send data

  // buffers for for allgather
  struct comm_buff_stat allgather_buff_local; 
  struct comm_buff_stat allgather_buff_global; 
};


