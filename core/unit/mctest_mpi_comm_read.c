#include <acutest.h>

#ifdef GKYL_HAVE_MPI

#include <math.h>
#include <mpi.h>
#include <stc/cstr.h>
#include <unistd.h>

#include <gkyl_array_rio.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_mpi_comm.h>
#include <gkyl_comm_io.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

void
mpi_read(int nrank, int cuts[2])
{
  int m_sz;
  MPI_Comm_size(MPI_COMM_WORLD, &m_sz);
  if (m_sz != nrank) return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  enum gkyl_array_rio_status status;

  struct gkyl_rect_grid grid;
  struct gkyl_array_header_info hdr;

  status = gkyl_grid_sub_array_header_read(&grid, &hdr,
    "core/data/unit/ser-euler_riem_2d_hllc-euler_1.gkyl");

  TEST_CHECK( GKYL_ARRAY_RIO_SUCCESS == status );

  if (hdr.meta_size > 0)
    free(hdr.meta);

  int nghost[] = { 2, 2 };
  struct gkyl_range global, ext_global;
  gkyl_create_grid_ranges(&grid, nghost, &ext_global, &global);

  size_t nc = hdr.esznc/gkyl_elem_type_size[hdr.etype];
  struct gkyl_array *s_arr = gkyl_array_new(hdr.etype, nc, ext_global.volume);
  gkyl_array_clear(s_arr, 0.0);

  status = gkyl_grid_sub_array_read(&grid, &global, s_arr,
    "core/data/unit/ser-euler_riem_2d_hllc-euler_1.gkyl");

  TEST_CHECK( GKYL_ARRAY_RIO_SUCCESS == status );

  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(global.ndim, cuts, &global);

  struct gkyl_comm *comm = gkyl_mpi_comm_new( &(struct gkyl_mpi_comm_inp) {
      .mpi_comm = MPI_COMM_WORLD,
      .decomp = decomp
    }
  );

  struct gkyl_range local, ext_local;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &ext_local, &local);

  struct gkyl_array *p_arr = gkyl_array_new(hdr.etype, nc, ext_local.volume);
  gkyl_array_clear(p_arr, 0.0);

  status = gkyl_comm_array_read(comm, &grid, &local, p_arr,
    "core/data/unit/euler_riem_2d_hllc-euler_1.gkyl");

  TEST_CHECK( GKYL_ARRAY_RIO_SUCCESS == status );

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    const double *s_dat = gkyl_array_cfetch(s_arr, gkyl_range_idx(&global, iter.idx));
    const double *p_dat = gkyl_array_cfetch(p_arr, gkyl_range_idx(&local, iter.idx));

    for (int c=0; c<nc; ++c)
      TEST_CHECK( gkyl_compare_double(s_dat[c], p_dat[c], 1e-14) );
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(s_arr);
  gkyl_array_release(p_arr);
}

void mpi_n1_read() { mpi_read(1, (int[]){1, 1}); }
void mpi_n2_read() { mpi_read(2, (int[]){2, 1}); }
void mpi_n4_read() { mpi_read(4, (int[]){2, 2}); }

TEST_LIST = {
  {"mpi_n1_read", mpi_n1_read},
  {"mpi_n2_read", mpi_n2_read},
  {"mpi_n4_read", mpi_n4_read},
  {NULL, NULL},
};

#else

// nothing to test if not building with MPI
TEST_LIST = {
  {NULL, NULL},
};

#endif
