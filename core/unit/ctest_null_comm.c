#include <acutest.h>

#include <gkyl_array_rio.h>
#include <gkyl_comm_io.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_null_comm.h>

void
test_1d()
{
  struct gkyl_range range;
  gkyl_range_init(&range, 1, (int[]) { 1 }, (int[]) { 100 });

  int cuts[] = { 1 };
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp
    }
  );

  TEST_CHECK( strcmp(comm->id, "null_comm") == 0 );
  TEST_CHECK( comm->has_decomp );

  int rank;
  gkyl_comm_get_rank(comm, &rank);
  TEST_CHECK( rank == 0 );

  int sz;
  gkyl_comm_get_size(comm, &sz);
  TEST_CHECK( sz == 1 );  

  double out[3], inp[3] = { 2.0, 4.0, 8.0 };
  gkyl_comm_allreduce(comm, GKYL_DOUBLE, GKYL_MIN, 3, inp, out);

  for (int i=0; i<3; ++i)
    TEST_CHECK( out[i] == inp[i] );

  int nghost[] = { 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr, 200005);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr, idx);
    f[0] = iter.idx[0];
  }

  int per_dirs[] = { 0 };
  gkyl_comm_array_per_sync(comm, &local, &local_ext, 1, per_dirs, arr );

  int idx[GKYL_MAX_DIM] = { 0 };
  
  for (int d=0; d<local.ndim; ++d) {
    int ncell = gkyl_range_shape(&local, d);

    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {

      if (!gkyl_range_contains_idx(&local, iter.idx)) {
        long lidx = gkyl_range_idx(&local_ext, iter.idx);
        
        for (int n=0; n<local.ndim; ++n)
          idx[n] = iter.idx[n];
        if (idx[d] > local.upper[d])
          idx[d] = idx[d] - ncell;
        else
          idx[d] = idx[d] + ncell;

        const double  *f = gkyl_array_cfetch(arr, lidx);
        TEST_CHECK( idx[0] == f[0] );
      }
    }
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
}

void
test_allgather_1d()
{
  struct gkyl_range range;
  gkyl_range_init(&range, 1, (int[]) { 1 }, (int[]) { 100 });

  int cuts[] = { 1 };
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp
    }
  );

  int rank, sz;
  gkyl_comm_get_rank(comm, &rank);
  gkyl_comm_get_size(comm, &sz);

  int nghost[] = { 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr, 200005);

  struct gkyl_range_iter iter;

  // Test gkyl_comm_array_allgather
  
  struct gkyl_array *arr_recv = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_comm_array_allgather(comm, &local, &local, arr, arr_recv);

  for (int d=0; d<local.ndim; ++d) {
    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&local_ext, iter.idx);
      const double *arr_c = gkyl_array_cfetch(arr, lidx);
      const double *arr_recv_c = gkyl_array_cfetch(arr_recv, lidx);
      TEST_CHECK( arr_c[d] == arr_recv_c[d] );
    }
  }

  // Test gkyl_comm_array_allgather_host

  gkyl_array_clear(arr_recv, 0);
  gkyl_comm_array_allgather_host(comm, &local, &local, arr, arr_recv);
  for (int d=0; d<local.ndim; ++d) {
    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&local_ext, iter.idx);
      const double *arr_c = gkyl_array_cfetch(arr, lidx);
      const double *arr_recv_c = gkyl_array_cfetch(arr_recv, lidx);
      TEST_CHECK( arr_c[d] == arr_recv_c[d] );
    }
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
  gkyl_array_release(arr_recv);
}

void
test_bcast_1d()
{
  struct gkyl_range range;
  gkyl_range_init(&range, 1, (int[]) { 1 }, (int[]) { 100 });

  int cuts[] = { 1 };
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp
    }
  );

  int rank, sz;
  gkyl_comm_get_rank(comm, &rank);
  gkyl_comm_get_size(comm, &sz);

  int nghost[] = { 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  struct gkyl_array *arr_recv = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);

  gkyl_array_clear(arr, 200005);

  struct gkyl_range_iter iter;
  // Test gkyl_comm_array_bcast 
  int bcast_rank = 0;
  gkyl_array_clear(arr_recv, 0);
  gkyl_comm_array_bcast(comm, arr, arr_recv, bcast_rank);
  for (int d=0; d<local.ndim; ++d) {
    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&local_ext, iter.idx);
      const double *arr_c = gkyl_array_cfetch(arr, lidx);
      const double *arr_recv_c = gkyl_array_cfetch(arr_recv, lidx);
      TEST_CHECK( arr_c[d] == arr_recv_c[d] );
    }
  }

  // Test gkyl_comm_array_bcast_host 

  gkyl_array_clear(arr_recv, 0);
  gkyl_comm_array_bcast_host(comm, arr, arr_recv, bcast_rank);
  for (int d=0; d<local.ndim; ++d) {
    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&local_ext, iter.idx);
      const double *arr_c = gkyl_array_cfetch(arr, lidx);
      const double *arr_recv_c = gkyl_array_cfetch(arr_recv, lidx);
      TEST_CHECK( arr_c[d] == arr_recv_c[d] );
    }
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
  gkyl_array_release(arr_recv);
}

void
test_2d()
{
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 4, 4 });

  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp
    }
  );

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_range local_x[2], local_ext_x[2];
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { nghost[0], 0 },
    &local_ext_x[0], &local_x[0]);
  
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { 0, nghost[1] },
    &local_ext_x[1], &local_x[1]);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr, 200005);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr, idx);

    for (int d=0; d<local.ndim; ++d)
      f[d] = iter.idx[d];
  }

  int per_dirs[] = { 0, 1 };
  gkyl_comm_array_per_sync(comm, &local, &local_ext, 2, per_dirs, arr );

  int idx[GKYL_MAX_DIM] = { 0 };
  int count = 0;
  
  for (int d=0; d<local.ndim; ++d) {
    int ncell = gkyl_range_shape(&local, d);

    gkyl_range_iter_init(&iter, &local_ext_x[d]);
    while (gkyl_range_iter_next(&iter)) {

      if (!gkyl_range_contains_idx(&local, iter.idx)) {
        long lidx = gkyl_range_idx(&local_ext, iter.idx);
        
        for (int n=0; n<local.ndim; ++n)
          idx[n] = iter.idx[n];
        if (idx[d] > local.upper[d])
          idx[d] = idx[d] - ncell;
        else
          idx[d] = idx[d] + ncell;

        const double  *f = gkyl_array_cfetch(arr, lidx);
        for (int n=0; n<local.ndim; ++n)
          TEST_CHECK( idx[n] == f[n] );
      }
    }
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
}

void
test_bcast_2d(){
  struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 4, 4 });

  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp
    }
  );

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_range local_x[2], local_ext_x[2];
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { nghost[0], 0 },
    &local_ext_x[0], &local_x[0]);
  
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { 0, nghost[1] },
    &local_ext_x[1], &local_x[1]);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  struct gkyl_array *arr_recv = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr, 200005);

  struct gkyl_range_iter iter;
  // Test gkyl_comm_array_bcast 
  int bcast_rank = 0;
  gkyl_array_clear(arr_recv, 0);
  gkyl_comm_array_bcast(comm, arr, arr_recv, bcast_rank);
  for (int d=0; d<local.ndim; ++d) {
    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&local_ext, iter.idx);
      const double *arr_c = gkyl_array_cfetch(arr, lidx);
      const double *arr_recv_c = gkyl_array_cfetch(arr_recv, lidx);
      TEST_CHECK( arr_c[d] == arr_recv_c[d] );
    }
  }

  // Test gkyl_comm_array_bcast_host 

  gkyl_array_clear(arr_recv, 0);
  gkyl_comm_array_bcast_host(comm, arr, arr_recv, bcast_rank);
  for (int d=0; d<local.ndim; ++d) {
    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&local_ext, iter.idx);
      const double *arr_c = gkyl_array_cfetch(arr, lidx);
      const double *arr_recv_c = gkyl_array_cfetch(arr_recv, lidx);
      TEST_CHECK( arr_c[d] == arr_recv_c[d] );
    }
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
  gkyl_array_release(arr_recv);
}

void
test_allgather_2d(){

    struct gkyl_range range;
  gkyl_range_init(&range, 2, (int[]) { 1, 1 }, (int[]) { 4, 4 });

  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp
    }
  );

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_range local_x[2], local_ext_x[2];
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { nghost[0], 0 },
    &local_ext_x[0], &local_x[0]);
  
  gkyl_create_ranges(&decomp->ranges[rank], (int[]) { 0, nghost[1] },
    &local_ext_x[1], &local_x[1]);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  struct gkyl_array *arr_recv = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr, 200005);

  struct gkyl_range_iter iter;

  // Test gkyl_comm_array_allgather

  gkyl_comm_array_allgather(comm, &local, &local, arr, arr_recv);

  for (int d=0; d<local.ndim; ++d) {
    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&local_ext, iter.idx);
      const double *arr_c = gkyl_array_cfetch(arr, lidx);
      const double *arr_recv_c = gkyl_array_cfetch(arr_recv, lidx);
      TEST_CHECK( arr_c[d] == arr_recv_c[d] );
    }
  }

  // Test gkyl_comm_array_allgather_host 

  gkyl_array_clear(arr_recv, 0);
  gkyl_comm_array_allgather_host(comm, &local, &local, arr, arr_recv);
  for (int d=0; d<local.ndim; ++d) {
    gkyl_range_iter_init(&iter, &local_ext);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&local_ext, iter.idx);
      const double *arr_c = gkyl_array_cfetch(arr, lidx);
      const double *arr_recv_c = gkyl_array_cfetch(arr_recv, lidx);
      TEST_CHECK( arr_c[d] == arr_recv_c[d] );
    }
  }

  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
  gkyl_array_release(arr_recv);
}

void
test_io_2d()
{
  int cells[] = { 32, 32 };
  struct gkyl_range range;
  gkyl_range_init_from_shape1(&range, 2, cells);

  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);

  struct gkyl_comm *comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp
    }
  );

  double lower[] = {0.0, 0.5}, upper[] = {1.0, 2.5};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int rank;
  gkyl_comm_get_rank(comm, &rank);

  int nghost[] = { 1, 1 };
  struct gkyl_range local, local_ext;
  gkyl_create_ranges(&decomp->ranges[rank], nghost, &local_ext, &local);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);
  gkyl_array_clear(arr, 1.5);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    double  *f = gkyl_array_fetch(arr, idx);

    double xc[GKYL_MAX_DIM] = { 0.0 };
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);
    f[0] = sin(2*M_PI*xc[0])*sin(2*M_PI*xc[1]);
    f[1] = cos(2*M_PI*xc[0])*sin(2*M_PI*xc[1]);
  }  

  int status;
  status = gkyl_comm_array_write(comm, &grid, &local, 0, arr, "ctest_null_comm_io_2d.gkyl");

  struct gkyl_array *arr_rw = gkyl_array_new(GKYL_DOUBLE, range.ndim, local_ext.volume);  
  status = gkyl_comm_array_read(comm, &grid, &local, arr_rw, "ctest_null_comm_io_2d.gkyl");
  TEST_CHECK( status == 0 );

  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long idx = gkyl_range_idx(&local, iter.idx);
    const double *f = gkyl_array_cfetch(arr, idx);
    const double *frw = gkyl_array_cfetch(arr_rw, idx);

    TEST_CHECK( gkyl_compare_double(f[0], frw[0], 1e-15) );
    TEST_CHECK( gkyl_compare_double(f[1], frw[1], 1e-15) );
  }
  
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(arr);
  gkyl_array_release(arr_rw);
}

void
test_io_p1_p4(void)
{
  struct gkyl_rect_grid grid;
  struct gkyl_array_header_info hdr;
  gkyl_grid_sub_array_header_read(&grid, &hdr,
    "core/data/unit/ser-euler_riem_2d_hllc-euler_1.gkyl");

  size_t nc = hdr.esznc/gkyl_elem_type_size[hdr.etype];  
  
  int nghost[] = { 1, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  int cuts[] = { 1, 1 };
  struct gkyl_rect_decomp *decomp =
    gkyl_rect_decomp_new_from_cuts(range.ndim, cuts, &range);
  
  struct gkyl_comm *comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp
    }
  );
  
  struct gkyl_array *s_arr = gkyl_array_new(hdr.etype, nc, ext_range.volume);
  gkyl_array_clear(s_arr, 0.0);

  int status;
  status = gkyl_grid_sub_array_read(&grid, &range, s_arr,
    "core/data/unit/ser-euler_riem_2d_hllc-euler_1.gkyl");

  TEST_CHECK( 0 == status );

  struct gkyl_array *p_arr = gkyl_array_new(hdr.etype, nc, ext_range.volume);
  gkyl_array_clear(p_arr, 0.0);

  status = gkyl_comm_array_read(comm, &grid, &range, p_arr,
    "core/data/unit/euler_riem_2d_hllc-euler_1.gkyl");

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    const double *s_dat = gkyl_array_cfetch(s_arr, gkyl_range_idx(&range, iter.idx));
    const double *p_dat = gkyl_array_cfetch(p_arr, gkyl_range_idx(&range, iter.idx));

    for (int c=0; c<nc; ++c)
      TEST_CHECK( gkyl_compare_double(s_dat[c], p_dat[c], 1e-14) );
  }

  TEST_CHECK( 0 == status );  

  gkyl_array_header_info_release(&hdr);
  gkyl_rect_decomp_release(decomp);
  gkyl_comm_release(comm);
  gkyl_array_release(s_arr);
  gkyl_array_release(p_arr);
}

TEST_LIST = {
  { "test_1d", test_1d },
  { "test_allgather_1d", test_allgather_1d },
  { "test_bcast_1d", test_bcast_1d },
  { "test_2d", test_2d },
  { "test_allgather_2d", test_allgather_2d },
  { "test_bcast_2d", test_bcast_2d },
  { "test_io_2d", test_io_2d },
  { "test_io_p1_p4", test_io_p1_p4 },
  { NULL, NULL },
};
