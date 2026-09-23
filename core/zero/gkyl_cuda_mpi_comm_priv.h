#pragma once

#include <gkyl_cuda_mpi_comm.h>

#if defined(GKYL_HAVE_MPI) && defined(GKYL_HAVE_CUDA)
#include <gkyl_alloc.h>
#include <gkyl_comm_priv.h>

struct cuda_mpi_msg {
  struct gkyl_range range;
  struct gkyl_array *array;
  gkyl_mem_buff buffer;
  int peer, tag, count;
};

struct cuda_mpi_comm {
  struct gkyl_comm_priv priv_comm;
  MPI_Comm mcomm;
  struct gkyl_comm *host;
  struct gkyl_rect_decomp *decomp;
  struct gkyl_rect_decomp_neigh *neigh, *periodic[GKYL_MAX_DIM];
  int rank, size, capacity;
  bool sync_corners;
  struct cuda_mpi_msg *send, *recv;
  MPI_Request *requests;
  int *counts, *offsets;
  gkyl_mem_buff gather_local[2], gather_global[2]; // Host and device buffers.
};
#endif
