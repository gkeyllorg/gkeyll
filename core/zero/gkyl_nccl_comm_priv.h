#pragma once

// Private header for nccl_comm. Do not include in user-facing header files.

#include <gkyl_nccl_comm.h>
#include <gkyl_alloc.h>

// Maximum number of recv neighbors: not sure hard-coding this is a
// good idea.
#define MAX_RECV_NEIGH 32

#define NCCL_BASE_TAG 4343
#define NCCL_BASE_PER_TAG 5353

// Some NCCL calls should return ncclSuccess when done properly,
// but others (e.g. ncclGroupEnd) may return ncclInProgress.
// If a function that should return ncclSuccess returns ncclInProgress
// for some reason, having only one check function may be a problem.
// We could create a separate check function which waits and times out
// after a set amount of time.
#define checkNCCL(cmd) do {                            \
  ncclResult_t res = cmd;                              \
  if (res != ncclSuccess  && res != ncclInProgress) {  \
    fprintf(stderr, "Failed, NCCL error %s:%d '%s'\n", \
        __FILE__,__LINE__,ncclGetErrorString(res));    \
    exit(EXIT_FAILURE);                                \
  }                                                    \
} while(0)

// Object with a range, a status and a buffer used for send/recv.
struct comm_buff_stat {
  struct gkyl_range range;
//  MPI_Request status;
  gkyl_mem_buff buff;
};

// Private struct wrapping NCCL-specific code
struct nccl_comm {
  struct gkyl_comm_priv priv_comm; // base communicator.

  int rank; // Process ID in this communicator.
  int size; // Size of this communicator.

  ncclComm_t ncomm; // NCCL communicator to use.
  MPI_Comm mcomm; // MPI communicator to use
  bool is_mcomm_allocated; // Is the mcomm allocated?
  struct gkyl_comm *mpi_comm; // MPI comm this NCCL comm derives from.
  bool has_decomp; // Whether this comm is associated with a decomposition (e.g. of a range)
  cudaStream_t custream; // Cuda stream for NCCL comms.
  struct gkyl_rect_decomp *decomp; // pre-computed decomposition
  bool sync_corners; // Whether to sync corners.
  long local_range_offset; // Offset of the local region.

  struct gkyl_rect_decomp_neigh *neigh; // neighbors of local region
  struct gkyl_rect_decomp_neigh *per_neigh[GKYL_MAX_DIM]; // periodic neighbors

  int nrecv; // number of elements in rinfo array
  struct comm_buff_stat recv[MAX_RECV_NEIGH]; // info for recv data

  int nsend; // number of elements in sinfo array
  struct comm_buff_stat send[MAX_RECV_NEIGH]; // info for send data

  struct gkyl_range dir_edge; // for use in computing tags
  int is_on_edge[2][GKYL_MAX_DIM]; // flags to indicate if local range is on edge
  bool touches_any_edge; // true if this range touches any edge

  // buffers for for allgather
  struct comm_buff_stat allgather_buff_local; 
  struct comm_buff_stat allgather_buff_global; 
};
