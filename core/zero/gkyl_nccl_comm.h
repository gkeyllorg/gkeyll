#pragma once

#ifdef GKYL_HAVE_NCCL

#include <mpi.h>
#include <nccl.h>
#include <gkyl_comm.h>
#include <gkyl_rect_decomp.h>

// input to create new NCCL communicator
struct gkyl_nccl_comm_inp {
  MPI_Comm mpi_comm; // MPI communicator to use
  const struct gkyl_rect_decomp *decomp; // pre-computed decomposition
  bool sync_corners; // should we sync corners?
  bool device_set; // Whether the device has been set already (should only happen once).
  cudaStream_t custream; // CUDA stream to use.
};

/**
 * Return a new NCCL communicator. The methods in the communicator work
 * on arrays defined on the decomp specified in the input. This
 * contract is not checked (and can't be checked) and so the user of
 * this object must ensure consistency of the arrays and the
 * decomposition.
 *
 * @param inp Input struct to use for initialization
 * @return New NCCL communicator
 */
struct gkyl_comm *gkyl_nccl_comm_new(const struct gkyl_nccl_comm_inp *inp);

/**
 * Select the CUDA device of the calling rank: its rank in mpi_comm modulo the
 * number of devices. gkyl_nccl_comm_new makes the same choice unless
 * device_set is true. Call this before creating device objects that are built
 * before the communicator, so that they live on the communicator's device.
 *
 * @param mpi_comm MPI communicator whose rank selects the device
 */
void gkyl_nccl_comm_set_device(MPI_Comm mpi_comm);

#endif
