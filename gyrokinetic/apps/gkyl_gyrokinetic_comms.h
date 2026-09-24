#pragma once

#include <gkyl_util.h>

#include <gkyl_null_comm.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#include <gkyl_mpi_comm.h>
#ifdef GKYL_HAVE_NCCL
#include <gkyl_nccl_comm.h>
#endif
#endif

/**
 * Create a new communicator for the gyrokinetic app.
 *
 * @param use_mpi Whether using MPI to run a parallel sim.
 * @param use_gpu Whether to run on GPU(s).
 * @param Pointer to place where to put error messages.
 * @return New gkyl_comm communicator object.
 */
struct gkyl_comm* 
gkyl_gyrokinetic_comms_new(bool use_mpi, bool use_gpu, FILE *iostream);

/**
 * Free gyrokinetic app decomp and comm objects.
 *
 * @param decomp Decomposition object.
 * @param comm Communicator object.
 */
void
gkyl_gyrokinetic_comms_release(struct gkyl_comm *comm);
