#include <acutest.h>

#include <gkyl_fv_proj.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

void evalFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = x*x;
}

void
test_1()
{
  double lower[] = {-2.0}, upper[] = {4.0};
  int cells[] = {2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  gkyl_fv_proj *fv_proj = gkyl_fv_proj_new(&grid, 2, 1, evalFunc, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, 1, arr_range.volume);

  // project distribution function on basis
  gkyl_fv_proj_advance(fv_proj, 0.0, &arr_range, distf);

  // left cell
  double *dfl = gkyl_array_fetch(distf, 0);
  TEST_CHECK( gkyl_compare(1.0, dfl[0], 1e-12) );

  // right cell
  double *dfr = gkyl_array_fetch(distf, 1);
  TEST_CHECK( gkyl_compare(7.0, dfr[0], 1e-12) );

  gkyl_fv_proj_release(fv_proj);
  gkyl_array_release(distf);
}

TEST_LIST = {
  { "test_1", test_1 },
  { NULL, NULL },
};
