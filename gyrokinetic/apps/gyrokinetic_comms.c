#include <gkyl_gyrokinetic_comms.h>
#include <assert.h>

struct gkyl_comm* 
gkyl_gyrokinetic_comms_new(bool use_mpi, bool use_gpu, FILE *iostream)
{
  // Construct communicator for use in app.
  struct gkyl_comm *comm = 0;

#ifdef GKYL_HAVE_MPI
  if (use_gpu && use_mpi) {
#ifdef GKYL_HAVE_NCCL
    comm = gkyl_nccl_comm_new( &(struct gkyl_nccl_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
#else
    fprintf(iostream, " Using -g and -M together requires NCCL.\n");
    assert(0 == 1);
#endif
  }
  else if (use_mpi) {
    comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
        .mpi_comm = MPI_COMM_WORLD,
      }
    );
  }
  else {
    comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
        .use_gpu = use_gpu
      }
    );
  }
#else
  comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = use_gpu
    }
  );
#endif

  return comm;
}

void
gkyl_gyrokinetic_comms_release(struct gkyl_comm *comm)
{
  if (comm != 0)
    gkyl_comm_release(comm);
}
