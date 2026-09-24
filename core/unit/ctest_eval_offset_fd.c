#include <acutest.h>

#include <gkyl_eval_offset_fd.h>
#include <gkyl_rect_decomp.h>

void elc_field_1d(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = x+0.1;
  fout[1] = x*x+0.2;
  fout[2] = x*x*x+0.3;
}

void test_1d()
{
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  struct gkyl_offset_descr offsets[] = {
    { 0.0 },
    { -0.5 },
    { -0.5 }    
  };
  
  struct gkyl_eval_offset_fd *ev = gkyl_eval_offset_fd_new(
    &(struct gkyl_eval_offset_fd_inp) {
      .grid = &grid,
      .num_ret_vals = 3,
      .offsets = offsets,
      .eval = elc_field_1d
    }
  );

  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  struct gkyl_array *elc_fld = gkyl_array_new(GKYL_DOUBLE, 3, arr_range.volume);

  gkyl_eval_offset_fd_advance(ev, 0.0, &arr_range, elc_fld);

  double elc_out[3];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &arr_range);

  gkyl_range_iter_next(&iter);
  const double *El = gkyl_array_cfetch(elc_fld, gkyl_range_idx(&arr_range, iter.idx));

  elc_field_1d(0.0, (double[]) { -1.0 }, elc_out, 0);
  TEST_CHECK( El[0] == elc_out[0] );
  
  elc_field_1d(0.0, (double[]) { -2.0 }, elc_out, 0);
  TEST_CHECK( El[1] == elc_out[1] );
  TEST_CHECK( El[2] == elc_out[2] );

  gkyl_range_iter_next(&iter);
  const double *Er = gkyl_array_cfetch(elc_fld, gkyl_range_idx(&arr_range, iter.idx));

  elc_field_1d(0.0, (double[]) { 1.0 }, elc_out, 0);
  TEST_CHECK( Er[0] == elc_out[0] );
  
  elc_field_1d(0.0, (double[]) { 0.0 }, elc_out, 0);
  TEST_CHECK( Er[1] == elc_out[1] );
  TEST_CHECK( Er[2] == elc_out[2] );  

  gkyl_eval_offset_fd_release(ev);
  gkyl_array_release(elc_fld);
}

TEST_LIST = {
  { "test_1d", test_1d },
  { NULL, NULL },  
};
