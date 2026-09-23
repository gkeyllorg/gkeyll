#include <gkyl_gyrokinetic_comms.h>
#include <gkyl_cuda_mpi_comm.h>
#include <assert.h>
#include <string.h>

struct gkyl_comm *gkyl_gyrokinetic_comms_new(bool use_mpi, bool use_gpu, FILE *iostream)
{
  // Construct communicator for use in app.
  struct gkyl_comm *comm = 0;

#ifdef GKYL_HAVE_MPI
  if (use_gpu && use_mpi) {
#ifdef GKYL_HAVE_CUDA
    const char *backend = getenv("GKYL_GPU_COMM");
    if (backend && strcmp(backend, "cuda_mpi") == 0) {
      return gkyl_cuda_mpi_comm_new(&(struct gkyl_cuda_mpi_comm_inp){.mpi_comm = MPI_COMM_WORLD});
    }
    if (backend && strcmp(backend, "nccl") != 0) {
      fprintf(iostream, " Unknown GKYL_GPU_COMM='%s'; choose cuda_mpi or nccl.\n", backend);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
#endif
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new(&(struct gkyl_nccl_comm_inp){.mpi_comm = MPI_COMM_WORLD});
#else
    fprintf(iostream, " Using -g and -M requires NCCL or GKYL_GPU_COMM=cuda_mpi.\n");
    assert(0 == 1);
#endif
  } else if (use_mpi) {
    comm = gkyl_mpi_comm_new(&(struct gkyl_mpi_comm_inp){.mpi_comm = MPI_COMM_WORLD});
  } else {
    comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){.use_gpu = use_gpu});
  }
#else
  comm = gkyl_null_comm_inew(&(struct gkyl_null_comm_inp){.use_gpu = use_gpu});
#endif

  return comm;
}

void gkyl_gyrokinetic_comms_release(struct gkyl_comm *comm)
{
  if (comm != 0) {
    gkyl_comm_release(comm);
  }
}
