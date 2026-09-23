#pragma once

#include <gkyl_util.h>

#if defined(GKYL_HAVE_MPI) && defined(GKYL_HAVE_CUDA)

#include <mpi.h>
#include <gkyl_comm.h>

struct gkyl_cuda_mpi_comm_inp {
  MPI_Comm mpi_comm;
  const struct gkyl_rect_decomp *decomp;
  bool sync_corners;
  bool device_set; // Caller has already selected its GPU.
};

/**
 * Construct a CUDA-aware MPI communicator (requires MPI device-buffer support).
 * The MPI communicator is duplicated and owned by the returned object. Select
 * the GPU before allocating arrays, or let this constructor select by shared
 * memory rank (respecting CUDA_VISIBLE_DEVICES). No NCCL dependency is required.
 *
 * Device methods complete before returning. They synchronize CUDA producers
 * before MPI accesses device memory; MPI itself is not assumed stream-aware.
 * Host methods and file I/O accept host arrays. Halo and gather buffers grow
 * on first use and are reused. Messages must fit MPI's int byte counts.
 * Device halo, multiblock and gather methods use Gkeyll's double-only range
 * packing kernels and require GKYL_DOUBLE arrays. Broadcasts accept all array
 * types. Reductions support INT, INT_64, LONG, FLOAT and DOUBLE.
 * Release with gkyl_comm_release before MPI_Finalize.
 */
struct gkyl_comm *gkyl_cuda_mpi_comm_new(const struct gkyl_cuda_mpi_comm_inp *inp);

#endif
