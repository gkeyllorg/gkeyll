#include <acutest.h>
#include <mpack.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_rio.h>
#include <gkyl_array_rio_format_desc.h>
#include <gkyl_array_rio_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

static void set_array_to_zero_ho(struct gkyl_array *arr)
{
  double *arr_d  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arr_d[i] = 0.0;
}

void test_array_0()
{
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, 200);
  gkyl_array_release(arr);
}

void test_array_base()
{
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, 200);

  TEST_CHECK( gkyl_array_is_using_buffer(arr) == false );

  TEST_CHECK( arr->type = GKYL_DOUBLE );
  TEST_CHECK( arr->elemsz == sizeof(double) );
  TEST_CHECK( arr->ncomp == 1 );
  TEST_CHECK( arr->size == 20*10 );
  TEST_CHECK( arr->ref_count.count == 1 );

  TEST_CHECK( arr->on_dev == arr );

  TEST_CHECK( gkyl_array_is_cu_dev(arr) == false );

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i){
    TEST_CHECK( arrData[i] == 0. );
    arrData[i] = (i+0.5)*0.1;
  }

  // clone array
  struct gkyl_array *brr = gkyl_array_clone(arr);

  TEST_CHECK( brr->elemsz == sizeof(double) );
  TEST_CHECK( arr->ncomp == 1 );  
  TEST_CHECK( brr->size == 20*10 );
  TEST_CHECK( brr->ref_count.count == 1 );

  double *brrData  = brr->data;
  for (unsigned i=0; i<brr->size; ++i)
    TEST_CHECK( brrData[i] == arrData[i] );

  // reset values in brr
  for (unsigned i=0; i<brr->size; ++i)
    brrData[i] = (i-0.5)*0.5;

  gkyl_array_copy(arr, brr);

  for (unsigned i=0; i<arr->size; ++i)
    TEST_CHECK( arrData[i] == brrData[i] );

  // acquire pointer
  struct gkyl_array *crr = gkyl_array_acquire(arr);

  TEST_CHECK( crr->ref_count.count == 2 );
  TEST_CHECK( arr->ref_count.count == 2 );

  struct gkyl_array *drr = gkyl_array_acquire(crr);

  TEST_CHECK( drr->ref_count.count == 3 );
  TEST_CHECK( crr->ref_count.count == 3 );  
  TEST_CHECK( arr->ref_count.count == 3 );  
  
  gkyl_array_release(crr);
  TEST_CHECK( arr->ref_count.count == 2 );
  gkyl_array_release(drr);
  TEST_CHECK( arr->ref_count.count == 1 );
  
  gkyl_array_release(arr);
  gkyl_array_release(brr);
}

void test_array_fetch()
{
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, 20);

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;

  double *arrDataLh = gkyl_array_fetch(arr, 0);
  TEST_CHECK( arrDataLh[0] == 0.05 );

  double *arrDataUh = gkyl_array_fetch(arr, 10);
  TEST_CHECK( arrDataUh[0] == (10+0.5)*0.1);
 
  gkyl_array_release(arr);
}

void test_non_numeric()
{
  struct euler { double rho, u, E; };
  struct gkyl_array *arr = gkyl_array_new(GKYL_USER, sizeof(struct euler), 10);

  for (unsigned i=0; i<arr->size; ++i) {
    struct euler *e = gkyl_array_fetch(arr, i);
    e->rho = 1.0; e->u = 0.0; e->E = 100.5;
  }
  
  struct gkyl_array *brr = gkyl_array_new(GKYL_USER, sizeof(struct euler), 10);
  gkyl_array_copy(brr, arr);

  for (unsigned i=0; i<arr->size; ++i) {
    struct euler *e = gkyl_array_fetch(brr, i);
    TEST_CHECK( e->rho == 1.0 );
    TEST_CHECK( e->u == 0.0 );
    TEST_CHECK( e->E == 100.5 );
  }  

  gkyl_array_release(arr);
  gkyl_array_release(brr);
}

void
test_grid_sub_array_read_1()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 60};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int nghost[] = { 1, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 2, ext_range.volume);
  set_array_to_zero_ho(arr);

  // set some values in array
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    double *d = gkyl_array_fetch(arr, loc);
    for (int k=0; k<2; ++k)
      d[k] = (10.5*iter.idx[0] + 220.5*iter.idx[1])*(k+0.5);
  }

  gkyl_grid_sub_array_write(&grid, &range, 0, arr, "ctest_grid_sub_array_1.gkyl");

  struct gkyl_array *arr2 = gkyl_array_new(GKYL_DOUBLE, 2, ext_range.volume);
  set_array_to_zero_ho(arr2);

  struct gkyl_rect_grid grid2;

  // read just header
  FILE *fp;
  with_file(fp, "ctest_grid_sub_array_1.gkyl", "r") {
    struct gkyl_array_header_info hdr;
    
    int status = gkyl_grid_sub_array_header_read_fp(&grid2, &hdr, fp);
    TEST_CHECK( status == 0 );

    TEST_CHECK( hdr.file_type == 1);
    TEST_CHECK( hdr.etype == GKYL_DOUBLE);

    long tot_cells = 1L;
    TEST_CHECK( grid.ndim == grid2.ndim );
    for (int d=0; d<grid.ndim; ++d) {
      TEST_CHECK( grid.lower[d] == grid2.lower[d] );
      TEST_CHECK( grid.upper[d] == grid2.upper[d] );
      TEST_CHECK( grid.cells[d] == grid2.cells[d] );
      TEST_CHECK( grid.dx[d] == grid2.dx[d] );

      tot_cells *= grid.cells[d];
    }
    TEST_CHECK( grid.cellVolume == grid2.cellVolume );

    TEST_CHECK( hdr.esznc = arr->esznc );
    TEST_CHECK( hdr.tot_cells == tot_cells );

    TEST_CHECK( 0 == hdr.meta_size );
    TEST_CHECK( 1 == hdr.nrange );
  }

  int file_type = gkyl_get_gkyl_file_type("ctest_grid_sub_array_1.gkyl");
  TEST_CHECK( 1 == file_type );
  
  // read back the grid and the array
  int err =
    gkyl_grid_sub_array_read(&grid2, &range, arr2, "ctest_grid_sub_array_1.gkyl");

  TEST_CHECK( err < 1 );
  
  if (err < 1) {

    TEST_CHECK( grid.ndim == grid2.ndim );
    for (int d=0; d<grid.ndim; ++d) {
      TEST_CHECK( grid.lower[d] == grid2.lower[d] );
      TEST_CHECK( grid.upper[d] == grid2.upper[d] );
      TEST_CHECK( grid.cells[d] == grid2.cells[d] );
      TEST_CHECK( grid.dx[d] == grid2.dx[d] );
    }
    TEST_CHECK( grid.cellVolume == grid2.cellVolume );
    
    gkyl_range_iter_init(&iter, &range);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&range, iter.idx);
      
      const double *rhs = gkyl_array_cfetch(arr, loc);
      const double *lhs = gkyl_array_cfetch(arr2, loc);
      for (int k=0; k<2; ++k)
        TEST_CHECK( lhs[k] == rhs[k] );
    }
  }
  
  gkyl_array_release(arr);
  gkyl_array_release(arr2);
}

void
test_grid_sub_array_read_2()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 60};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int nghost[] = { 1, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 2, ext_range.volume);
  set_array_to_zero_ho(arr);

  // set some values in array
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    double *d = gkyl_array_fetch(arr, loc);
    for (int k=0; k<2; ++k)
      d[k] = (10.5*iter.idx[0] + 220.5*iter.idx[1])*(k+0.5);
  }

  gkyl_grid_sub_array_write(&grid, &range, 0, arr, "ctest_grid_sub_array_2.gkyl");

  struct gkyl_rect_grid grid2;

  struct gkyl_range srange;
  gkyl_range_init(&srange, grid.ndim, (int[]) { 5, 5 }, (int[]) { 10, 15 });

  struct gkyl_array *arr2 = gkyl_array_new(GKYL_DOUBLE, 2, srange.volume);
  set_array_to_zero_ho(arr2);

  // read back the grid and the array
  int err =
    gkyl_grid_sub_array_read(&grid2, &srange, arr2, "ctest_grid_sub_array_2.gkyl");

  TEST_CHECK( err < 1 );
  
  if (err < 1) {

    TEST_CHECK( grid.ndim == grid2.ndim );
    for (int d=0; d<grid.ndim; ++d) {
      TEST_CHECK( grid.lower[d] == grid2.lower[d] );
      TEST_CHECK( grid.upper[d] == grid2.upper[d] );
      TEST_CHECK( grid.cells[d] == grid2.cells[d] );
      TEST_CHECK( grid.dx[d] == grid2.dx[d] );
    }
    TEST_CHECK( grid.cellVolume == grid2.cellVolume );
    
    gkyl_range_iter_init(&iter, &srange);
    while (gkyl_range_iter_next(&iter)) {
      const double *rhs = gkyl_array_cfetch(arr, gkyl_range_idx(&range, iter.idx));
      const double *lhs = gkyl_array_cfetch(arr2, gkyl_range_idx(&srange, iter.idx));
      for (int k=0; k<2; ++k)
        TEST_CHECK( lhs[k] == rhs[k] );
    }
  }
  
  gkyl_array_release(arr);
  gkyl_array_release(arr2);
}

void
test_grid_array_new_from_file_1()
{
  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 60};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  int nghost[] = { 1, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 2, ext_range.volume);
  set_array_to_zero_ho(arr);

  // set some values in array
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    double *d = gkyl_array_fetch(arr, loc);
    for (int k=0; k<2; ++k)
      d[k] = (10.5*iter.idx[0] + 220.5*iter.idx[1])*(k+0.5);
  }

  gkyl_grid_sub_array_write(&grid, &range, 0, arr, "ctest_grid_array_new_from_file_1.gkyl");

  struct gkyl_rect_grid grid2;
  struct gkyl_array *arr2 = gkyl_grid_array_new_from_file(&grid2, "ctest_grid_array_new_from_file_1.gkyl");

  TEST_CHECK( arr2->type == GKYL_DOUBLE );

  TEST_CHECK( grid.ndim == grid2.ndim );
  for (int d=0; d<grid.ndim; ++d) {
    TEST_CHECK( grid.lower[d] == grid2.lower[d] );
    TEST_CHECK( grid.upper[d] == grid2.upper[d] );
    TEST_CHECK( grid.cells[d] == grid2.cells[d] );
    TEST_CHECK( grid.dx[d] == grid2.dx[d] );
  }
  TEST_CHECK( grid.cellVolume == grid2.cellVolume );  

  // NOTE: array read from file is not the same shape as "arr". This
  // is because the new_from_file does not read ghost cells.
  long lhs_loc = 0;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    const double *rhs = gkyl_array_cfetch(arr, loc);
    const double *lhs = gkyl_array_cfetch(arr2, lhs_loc++);
    for (int k=0; k<2; ++k)
      TEST_CHECK( lhs[k] == rhs[k] );
  }
  
  gkyl_array_release(arr);
  gkyl_array_release(arr2);
}

void
test_grid_array_read_p1(void)
{
  // read just header
  struct gkyl_rect_grid grid;  
  struct gkyl_array_header_info hdr;
  FILE *fp = 0;  
  with_file(fp, "core/data/unit/ser-euler_riem_2d_hllc-euler_1.gkyl", "r") {
    
    int status = gkyl_grid_sub_array_header_read_fp(&grid, &hdr, fp);
    TEST_CHECK( status == 0 );

    TEST_CHECK( hdr.file_type == gkyl_file_type_int[GKYL_FIELD_DATA_FILE]);
    TEST_CHECK( hdr.etype == GKYL_DOUBLE);

    TEST_CHECK( hdr.esznc = 5*sizeof(double));
    TEST_CHECK( 50*50 == hdr.tot_cells );

    TEST_CHECK( 50 == grid.cells[0] );
    TEST_CHECK( 50 == grid.cells[1] );
    
    TEST_CHECK( 0.0 == grid.lower[0] );
    TEST_CHECK( 0.0 == grid.lower[1] );

    TEST_CHECK( 1.0 == grid.upper[0] );
    TEST_CHECK( 1.0 == grid.upper[1] );

    TEST_CHECK( 1 == hdr.nrange );

    TEST_CHECK( hdr.meta_size > 0 );

    mpack_tree_t tree;
    mpack_tree_init_data(&tree, hdr.meta, hdr.meta_size);
    mpack_tree_parse(&tree);

    mpack_node_t root = mpack_tree_root(&tree);
    TEST_CHECK(mpack_node_type(root) == mpack_type_map);

    mpack_node_t tm_node = mpack_node_map_cstr(root, "time");
    TEST_CHECK( mpack_node_double(tm_node) >  0.80675 );

    mpack_node_t fr_node = mpack_node_map_cstr(root, "frame");
    TEST_CHECK( mpack_node_i64(fr_node) == 1 );

    status = mpack_tree_destroy(&tree);
    TEST_CHECK( mpack_ok == status );

    gkyl_free(hdr.meta);
  }

  size_t nc = hdr.esznc/gkyl_elem_type_size[hdr.etype];

  int nghost[] = { 1, 2 };
  struct gkyl_range range, ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &ext_range, &range);

  // read serial data for comparison
  struct gkyl_rect_grid s_grid;
  struct gkyl_array *s_arr = gkyl_array_new(hdr.etype, nc, ext_range.volume);
  int s_status = gkyl_grid_sub_array_read(&s_grid, &range, s_arr,
    "core/data/unit/ser-euler_riem_2d_hllc-euler_1.gkyl");

  // read parallel data (whole domain)
  do {
    struct gkyl_rect_grid p_grid;  
    struct gkyl_array *p_arr = gkyl_array_new(hdr.etype,  nc, ext_range.volume);
    int p_status = gkyl_grid_sub_array_read(&p_grid, &range, p_arr,
      "core/data/unit/euler_riem_2d_hllc-euler_1.gkyl");

    TEST_CHECK( 0 == p_status );
    
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &range);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&range, iter.idx);
      const double *s_dat = gkyl_array_fetch(s_arr, loc);
      const double *p_dat = gkyl_array_fetch(p_arr, loc);
      
      for (int c=0; c<nc; ++c)
        TEST_CHECK( gkyl_compare_double(s_dat[c], p_dat[c], 1e-15) );
    }
    gkyl_array_release(p_arr);
  } while (0);

  // read parallel data (partial domain)
  do {
    struct gkyl_range prange;
    gkyl_range_init(&prange, 2, (int[]) { 10, 10 }, (int[]) { 30, 40 });
    
    struct gkyl_rect_grid p_grid;
    struct gkyl_array *p_arr = gkyl_array_new(hdr.etype,  nc, prange.volume);
    int p_status = gkyl_grid_sub_array_read(&p_grid, &prange, p_arr,
      "core/data/unit/euler_riem_2d_hllc-euler_1.gkyl");

    TEST_CHECK( 0 == p_status );
    
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &prange);
    while (gkyl_range_iter_next(&iter)) {

      const double *s_dat = gkyl_array_fetch(s_arr, gkyl_range_idx(&range, iter.idx));
      const double *p_dat = gkyl_array_fetch(p_arr, gkyl_range_idx(&prange, iter.idx));
      
      for (int c=0; c<nc; ++c)
        TEST_CHECK( gkyl_compare_double(s_dat[c], p_dat[c], 1e-15) );
    }
    gkyl_array_release(p_arr);
  } while (0);

  // read parallel data (partial domain)
  do {
    struct gkyl_range prange;
    gkyl_range_init(&prange, 2, (int[]) { 4, 5 }, (int[]) { 10, 10 });
    
    struct gkyl_rect_grid p_grid;
    struct gkyl_array *p_arr = gkyl_array_new(hdr.etype,  nc, prange.volume);
    int p_status = gkyl_grid_sub_array_read(&p_grid, &prange, p_arr,
      "core/data/unit/euler_riem_2d_hllc-euler_1.gkyl");

    TEST_CHECK( 0 == p_status );
    
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &prange);
    while (gkyl_range_iter_next(&iter)) {

      const double *s_dat = gkyl_array_fetch(s_arr, gkyl_range_idx(&range, iter.idx));
      const double *p_dat = gkyl_array_fetch(p_arr, gkyl_range_idx(&prange, iter.idx));
      
      for (int c=0; c<nc; ++c)
        TEST_CHECK( gkyl_compare_double(s_dat[c], p_dat[c], 1e-15) );
    }
    gkyl_array_release(p_arr);
  } while (0);

  gkyl_array_release(s_arr);
}

static void
test_array_from_buff(void)
{
  double *buff = gkyl_malloc(sizeof(double[400]));
  
  struct gkyl_array *arr = gkyl_array_new_from_buff(GKYL_DOUBLE, 1, 200, buff);

  TEST_CHECK( gkyl_array_is_using_buffer(arr) == true );

  TEST_CHECK( arr->type = GKYL_DOUBLE );
  TEST_CHECK( arr->elemsz == sizeof(double) );
  TEST_CHECK( arr->ncomp == 1 );
  TEST_CHECK( arr->size == 20*10 );
  TEST_CHECK( arr->ref_count.count == 1 );

  TEST_CHECK( arr->on_dev == arr );

  TEST_CHECK( gkyl_array_is_cu_dev(arr) == false );

  set_array_to_zero_ho(arr);
  
  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i){
    TEST_CHECK( arrData[i] == 0. );
    arrData[i] = (i+0.5)*0.1;
  }

  // clone array
  struct gkyl_array *brr = gkyl_array_clone(arr);

  TEST_CHECK( brr->elemsz == sizeof(double) );
  TEST_CHECK( arr->ncomp == 1 );  
  TEST_CHECK( brr->size == 20*10 );
  TEST_CHECK( brr->ref_count.count == 1 );

  double *brrData  = brr->data;
  for (unsigned i=0; i<brr->size; ++i)
    TEST_CHECK( brrData[i] == arrData[i] );

  // reset values in brr
  for (unsigned i=0; i<brr->size; ++i)
    brrData[i] = (i-0.5)*0.5;

  gkyl_array_copy(arr, brr);

  for (unsigned i=0; i<arr->size; ++i)
    TEST_CHECK( arrData[i] == brrData[i] );

  // acquire pointer
  struct gkyl_array *crr = gkyl_array_acquire(arr);

  TEST_CHECK( crr->ref_count.count == 2 );
  TEST_CHECK( arr->ref_count.count == 2 );

  struct gkyl_array *drr = gkyl_array_acquire(crr);

  TEST_CHECK( drr->ref_count.count == 3 );
  TEST_CHECK( crr->ref_count.count == 3 );  
  TEST_CHECK( arr->ref_count.count == 3 );
  
  gkyl_array_release(crr);
  TEST_CHECK( arr->ref_count.count == 2 );
  gkyl_array_release(drr);
  TEST_CHECK( arr->ref_count.count == 1 );

  gkyl_free(buff);
  
  gkyl_array_release(arr);
  gkyl_array_release(brr);
}

// Cuda specific tests
#ifdef GKYL_HAVE_CUDA

/* Function signatures of kernel calls */
int cu_array_test_and_flip_sign( struct gkyl_array *arr);

void test_cu_array_base()
{
  struct gkyl_array *arr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 200);

  TEST_CHECK( arr_cu->type = GKYL_DOUBLE );
  TEST_CHECK( arr_cu->elemsz == sizeof(double) );
  TEST_CHECK( arr_cu->ncomp == 1 );
  TEST_CHECK( arr_cu->size == 20*10 );
  TEST_CHECK( arr_cu->ref_count.count == 1 );

  TEST_CHECK( gkyl_array_is_cu_dev(arr_cu) == true );

  // create host array and initialize it
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, 200);

  gkyl_array_copy(arr, arr_cu);

  double *arrData = arr->data;
  for (unsigned i=0; i<arr->size; ++i) {
    TEST_CHECK( arrData[i] == 0. );
    arrData[i] = (i+0.5)*0.1;  
  }

  // copy to device
  gkyl_array_copy(arr_cu, arr);

  // reset host array
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = 0.0;

  // copy from device and check if things are ok
  gkyl_array_copy(arr, arr_cu);

  for (unsigned i=0; i<arr->size; ++i)
    TEST_CHECK( arrData[i] == (i+0.5)*0.1 );

  gkyl_array_release(arr);
  gkyl_array_release(arr_cu);
}

void test_cu_array_dev_kernel()
{
  // create a host array struct containing device data
  struct gkyl_array *arr_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, 1, 20);

  TEST_CHECK( arr_cu->type = GKYL_DOUBLE );
  TEST_CHECK( arr_cu->elemsz == sizeof(double) );
  TEST_CHECK( arr_cu->ncomp == 1 );
  TEST_CHECK( arr_cu->size == 20 );
  TEST_CHECK( arr_cu->ref_count.count == 1 );
  TEST_CHECK( gkyl_array_is_cu_dev(arr_cu) == true );

  // create host array and initialize it
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, arr_cu->ncomp, arr_cu->size);

  double *arrData  = arr->data;
  for (unsigned i=0; i<arr->size; ++i)
    arrData[i] = (i+0.5)*0.1;  

  // copy arr data to device data in arr_cu
  gkyl_array_copy(arr_cu, arr);

  // create new struct arr_cu_cl that is a clone of arr_cu
  struct gkyl_array *arr_cu_cl = gkyl_array_clone(arr_cu);

  // check arr_cu on device and flip sign
  int nfail = cu_array_test_and_flip_sign(arr_cu->on_dev);
  TEST_CHECK( nfail == 0 );

  // restore arr_cu by copying from arr_cu_cl, and test again
  gkyl_array_copy(arr_cu, arr_cu_cl);
  nfail = cu_array_test_and_flip_sign(arr_cu->on_dev);
  TEST_CHECK( nfail == 0 );

  // check arr_cu_cl on device and flip sign
  nfail = cu_array_test_and_flip_sign(arr_cu_cl->on_dev);
  TEST_CHECK( nfail == 0 );

  // copy arr_cu back to host and check
  gkyl_array_copy(arr, arr_cu);
  for (unsigned i=0; i<arr->size; ++i)
    TEST_CHECK( arrData[i] == -(i+0.5)*0.1);  

  // copy arr_cu_cl back to host and check
  // zero out arr first (no cheating)
  set_array_to_zero_ho(arr);
  gkyl_array_copy(arr, arr_cu_cl);
  for (unsigned i=0; i<arr->size; ++i)
    TEST_CHECK( arrData[i] == -(i+0.5)*0.1);  

  // release all data
  gkyl_array_release(arr_cu_cl);  
  gkyl_array_release(arr);
  gkyl_array_release(arr_cu);  
}

#endif

TEST_LIST = {
  { "array_0", test_array_0 },  
  { "array_base", test_array_base },
  { "array_fetch", test_array_fetch },
  { "non_numeric", test_non_numeric },
  { "grid_sub_array_read_1", test_grid_sub_array_read_1 },
  { "grid_sub_array_read_2", test_grid_sub_array_read_2 },  
  { "grid_array_new_from_file_1", test_grid_array_new_from_file_1 },
  { "grid_array_read_1", test_grid_array_read_p1 },
  { "array_from_buff", test_array_from_buff },  
#ifdef GKYL_HAVE_CUDA
  { "cu_array_base", test_cu_array_base },
  { "cu_array_dev_kernel", test_cu_array_dev_kernel },
#endif
  { NULL, NULL },
};
