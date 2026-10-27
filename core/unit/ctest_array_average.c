/* Test gkyl_array_average updater
This unit test performs several weighted averages and integrations from 1x to 3x arrays in various orders
with the gkyl_array_average updater. It compares the final result (full averaging) to the gkyl_array_integrate
updater.
*/

#include <acutest.h>
#include <gkyl_rect_grid.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_average.h>
#include <math.h>
#include <assert.h>
#include <gkyl_dg_bin_ops.h>

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(long nc, long size, bool use_gpu)
{
  return use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) :
                   gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Compare the computed result with the average computed with another updater.
double
solution_array_integrate(
  struct gkyl_rect_grid grid, struct gkyl_basis basis, struct gkyl_range local_ext,
  struct gkyl_range local, struct gkyl_array *win, struct gkyl_array *fin, bool use_gpu
)
{
  double *avgf_ref = use_gpu ? gkyl_cu_malloc(sizeof(double)) : gkyl_malloc(sizeof(double));

  if (win) {
    gkyl_dg_mul_op_range(&basis, 0, fin, 0, win, 0, fin, &local_ext);
  }

  struct gkyl_array_integrate *arr_integ =
    gkyl_array_integrate_new(&grid, &basis, 1, GKYL_ARRAY_INTEGRATE_OP_NONE, use_gpu);

  gkyl_array_integrate_advance(arr_integ, fin, 1.0, fin, &local, &local, avgf_ref);

  gkyl_array_integrate_release(arr_integ);

  double *avgf_ref_ho = gkyl_malloc(sizeof(double));
  if (use_gpu) {
    gkyl_cu_memcpy(avgf_ref_ho, avgf_ref, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(avgf_ref_ho, avgf_ref, sizeof(double));
  }

  double out = avgf_ref_ho[0];
  if (use_gpu) {
    gkyl_cu_free(avgf_ref);
  } else {
    gkyl_free(avgf_ref);
  }
  gkyl_free(avgf_ref_ho);
  return out;
}

// test 1x
void
evalFunc_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double lower[] = {-4.0}, upper[] = {6.0}; // Has to match the test below.
  double Lx = upper[0] - lower[0];
  double k_x = 2. * M_PI / Lx;
  double phi = 0.5;

  fout[0] = x * sin(k_x * x + phi);
}

void
evalWeight_1x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = 1 + x * x;
}
// direct weighted averaging x -> avg
void
test_1x(int poly_order, bool use_gpu)
{
  // define grid and basis
  double lower[] = {-4.0}, upper[] = {6.0};
  int cells[] = {16};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // define the reduced range and basis for averaging
  struct gkyl_range red_local, red_local_ext;
  gkyl_range_init(&red_local, 1, &local.lower[0], &local.lower[0]);
  gkyl_range_init(&red_local_ext, 1, &local_ext.lower[0], &local_ext.lower[0]);
  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  // project the target function and weight
  gkyl_proj_on_basis *projf =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalFunc_1x, NULL);

  struct gkyl_array *fx_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *fx_c_ho = use_gpu ? mkarr(fx_c->ncomp, fx_c->size, false) :
                                         gkyl_array_acquire(fx_c);
  gkyl_proj_on_basis_advance(projf, 0.0, &local, fx_c_ho);
  gkyl_array_copy(fx_c, fx_c_ho);

  gkyl_proj_on_basis_release(projf);

  gkyl_proj_on_basis *projw =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalWeight_1x, NULL);

  struct gkyl_array *wx_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *wx_c_ho = use_gpu ? mkarr(wx_c->ncomp, wx_c->size, false) :
                                         gkyl_array_acquire(wx_c);
  gkyl_proj_on_basis_advance(projw, 0.0, &local_ext, wx_c_ho);
  gkyl_array_copy(wx_c, wx_c_ho);

  gkyl_proj_on_basis_release(projw);

  // compute weighted average
  int avg_dim_x[] = {1, 0, 0};
  struct gkyl_array_average_inp inp_avg_full = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = red_basis,
    .local = &local,
    .local_avg = &red_local,
    .local_avg_ext = &red_local_ext,
    .weight = wx_c,
    .avg_dim = avg_dim_x,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *avg_full = gkyl_array_average_inew(&inp_avg_full);

  struct gkyl_array *avgf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(avg_full, fx_c, avgf_c);

  gkyl_array_average_release(avg_full);

  // fetch and transfer results
  double *avg_c0 = gkyl_array_fetch(avgf_c, 0);
  double *avg_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu) {
    gkyl_cu_memcpy(avg_c0_ho, avg_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(avg_c0_ho, avg_c0, sizeof(double));
  }

  // check results
  double intf_ref = solution_array_integrate(grid, basis, local_ext, local, wx_c, fx_c, use_gpu);
  double intw_ref = solution_array_integrate(grid, basis, local_ext, local, NULL, wx_c, use_gpu);
  double solution = intf_ref / intw_ref;
  double result = avg_c0_ho[0] * 0.5 * sqrt(2);
  double rel_err = fabs(result - solution) / fabs(solution);
  TEST_CHECK(gkyl_compare(result, solution, 1e-12));

  // clean up
  gkyl_array_release(avgf_c);
  gkyl_array_release(fx_c);
  gkyl_array_release(wx_c);
  gkyl_array_release(fx_c_ho);
  gkyl_array_release(wx_c_ho);
  gkyl_free(avg_c0_ho);
}

// tests 2x
void
evalFunc_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  double lower[] = {-4., -3.}, upper[] = {6., 5.};
  double Lx = upper[0] - lower[0], Ly = upper[1] - lower[1];
  double k_x = 2. * M_PI / Lx, k_y = 2. * M_PI / Ly;
  double phi = 0.5;
  fout[0] = 1 + sin(k_x * x + k_y * y);
  fout[0] = x * y * sin(1.5 * k_x * x + 0.75 * k_y * y + phi) * cos(1.42 * k_y * y);
}
void
evalWeight_2x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = 1 + x * x + y * y;
}

// one step weighted averaging x,y -> avg
void
test_2x_1step(int poly_order, bool use_gpu)
{
  double lower[] = {-4.0, -3.0}, upper[] = {6.0, 5.0};
  int cells[] = {16, 8};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // define the reduced range and basis for averaging
  struct gkyl_range red_local, red_local_ext;
  gkyl_range_init(&red_local, 1, &local.lower[0], &local.lower[0]);
  gkyl_range_init(&red_local_ext, 1, &local_ext.lower[0], &local_ext.lower[0]);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  // project the target and weight functions
  gkyl_proj_on_basis *projf =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalFunc_2x, NULL);

  struct gkyl_array *fxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);

  struct gkyl_array *fxy_c_ho = use_gpu ? mkarr(fxy_c->ncomp, fxy_c->size, false) :
                                          gkyl_array_acquire(fxy_c);
  gkyl_proj_on_basis_advance(projf, 0.0, &local_ext, fxy_c_ho);
  gkyl_array_copy(fxy_c, fxy_c_ho);

  gkyl_proj_on_basis_release(projf);

  gkyl_proj_on_basis *projw =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalWeight_2x, NULL);

  struct gkyl_array *wxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);

  struct gkyl_array *wxy_c_ho = use_gpu ? mkarr(wxy_c->ncomp, wxy_c->size, false) :
                                          gkyl_array_acquire(wxy_c);
  gkyl_proj_on_basis_advance(projw, 0.0, &local_ext, wxy_c_ho);
  gkyl_array_copy(wxy_c, wxy_c_ho);

  gkyl_proj_on_basis_release(projw);

  // perform the one step average
  int avg_dim_xy[] = {1, 1, 0};
  struct gkyl_array_average_inp inp_avg_xy = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = red_basis,
    .local = &local,
    .local_avg = &red_local,
    .local_avg_ext = &red_local_ext,
    .weight = wxy_c,
    .avg_dim = avg_dim_xy,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *avg_xy = gkyl_array_average_inew(&inp_avg_xy);

  struct gkyl_array *avgf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(avg_xy, fxy_c, avgf_c);

  gkyl_array_average_release(avg_xy);

  // check results
  double *avgf_c0 = gkyl_array_fetch(avgf_c, 0);
  double avgf_c0_ho[1];
  if (use_gpu) {
    gkyl_cu_memcpy(avgf_c0_ho, avgf_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(avgf_c0_ho, avgf_c0, sizeof(double));
  }
  double result = avgf_c0_ho[0] * sqrt(2) / 2;

  double intwf_ref = solution_array_integrate(grid, basis, local_ext, local, wxy_c, fxy_c, use_gpu);
  double intw_ref = solution_array_integrate(grid, basis, local_ext, local, NULL, wxy_c, use_gpu);
  double solution = intwf_ref / intw_ref;

  double rel_err = fabs(result - solution) / fabs(solution);
  TEST_CHECK(gkyl_compare(rel_err, 0, 1e-12));

  // clean up
  gkyl_array_release(avgf_c);
  gkyl_array_release(fxy_c);
  gkyl_array_release(wxy_c);
  gkyl_array_release(fxy_c_ho);
  gkyl_array_release(wxy_c_ho);
}

// two step integration of the weight x,y -> y -> int
void
test_2x_intx_inty(int poly_order, bool use_gpu)
{
  // define grids and basis
  double lower[] = {-4.0, -3.0}, upper[] = {6.0, 5.0};
  int cells[] = {4, 2};
  int ghost[] = {1, 1};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_rect_grid grid_y;
  gkyl_rect_grid_init(&grid_y, 1, &lower[1], &upper[1], &cells[1]);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  struct gkyl_basis basis_y;
  gkyl_cart_modal_serendip(&basis_y, 1, poly_order);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_range local_y, local_y_ext;
  int ghost_y[] = {ghost[1]};
  gkyl_create_grid_ranges(&grid_y, ghost_y, &local_y_ext, &local_y);

  struct gkyl_range red_local, red_local_ext;
  gkyl_range_init(&red_local, 1, &local.lower[0], &local.lower[0]);
  gkyl_range_init(&red_local_ext, 1, &local_ext.lower[0], &local_ext.lower[0]);

  gkyl_proj_on_basis *projw =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalWeight_2x, NULL);

  struct gkyl_array *wxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);

  struct gkyl_array *wxy_c_ho = use_gpu ? mkarr(wxy_c->ncomp, wxy_c->size, false) :
                                          gkyl_array_acquire(wxy_c);
  gkyl_proj_on_basis_advance(projw, 0.0, &local_ext, wxy_c_ho);
  gkyl_array_copy(wxy_c, wxy_c_ho);

  gkyl_proj_on_basis_release(projw);

  // integration over x only, (x,y) to (y)
  int int_dim_x[] = {1, 0, 0};
  struct gkyl_array_average_inp inp_int_x = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = basis_y,
    .local = &local,
    .local_avg = &local_y,
    .local_avg_ext = &local_y_ext,
    .weight = NULL,
    .avg_dim = int_dim_x,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *int_x = gkyl_array_average_inew(&inp_int_x);

  struct gkyl_array *wy_c = mkarr(basis_y.num_basis, local_y_ext.volume, use_gpu);
  gkyl_array_average_advance(int_x, wxy_c, wy_c);

  gkyl_array_average_release(int_x);

  // integration over remaining dimensions (y)
  int int_dim_y[] = {1, 0, 0};
  struct gkyl_array_average_inp inp_int_y = {
    .grid = &grid_y,
    .basis = basis_y,
    .basis_avg = red_basis,
    .local = &local_y,
    .local_avg = &red_local,
    .local_avg_ext = &red_local_ext,
    .weight = NULL,
    .avg_dim = int_dim_y,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *int_y = gkyl_array_average_inew(&inp_int_y);

  struct gkyl_array *intw_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(int_y, wy_c, intw_c);

  gkyl_array_average_release(int_y);

  // multiply by the volume to get the integral
  double volume = 1;
  for (int d = 0; d < ndim; d++) {
    volume *= grid.upper[d] - grid.lower[d];
  }
  gkyl_array_scale(intw_c, volume);

  // retrieve the computed average from the device (if applicable)
  double *intw_c0 = gkyl_array_fetch(intw_c, 0);
  double intw_c0_ho[1];
  if (use_gpu) {
    gkyl_cu_memcpy(intw_c0_ho, intw_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(intw_c0_ho, intw_c0, sizeof(double));
  }

  double result = intw_c0_ho[0];

  double intw_ref = solution_array_integrate(grid, basis, local_ext, local, NULL, wxy_c, use_gpu);
  double solution = intw_ref;

  double rel_err = fabs(result - solution) / fabs(solution);
  TEST_CHECK(gkyl_compare(rel_err, 0.0, 1e-12));

  // clean up
  gkyl_array_release(wxy_c);
  gkyl_array_release(wy_c);
  gkyl_array_release(intw_c);
  gkyl_array_release(wxy_c_ho);
}

// two steps averaging x,y -> y -> avg
void
test_2x_avgx_avgy(int poly_order, bool use_gpu)
{
  // define grids and basis
  double lower[] = {-4.0, -3.0}, upper[] = {6.0, 5.0};
  int cells[] = {16, 8};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_rect_grid grid_y;
  gkyl_rect_grid_init(&grid_y, 1, &lower[1], &upper[1], &cells[1]);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  struct gkyl_basis basis_y;
  gkyl_cart_modal_serendip(&basis_y, 1, poly_order);

  int ghost[] = {0, 0};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  struct gkyl_range local_y, local_y_ext;
  int ghost_y[] = {ghost[1]};
  gkyl_create_grid_ranges(&grid_y, ghost_y, &local_y_ext, &local_y);

  struct gkyl_range red_local, red_local_ext;
  gkyl_range_init(&red_local, 1, &local.lower[0], &local.lower[0]);
  gkyl_range_init(&red_local_ext, 1, &local_ext.lower[0], &local_ext.lower[0]);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  // project the target function and weight
  gkyl_proj_on_basis *projf =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalFunc_2x, NULL);

  struct gkyl_array *fxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);

  struct gkyl_array *fxy_c_ho = use_gpu ? mkarr(fxy_c->ncomp, fxy_c->size, false) :
                                          gkyl_array_acquire(fxy_c);
  gkyl_proj_on_basis_advance(projf, 0.0, &local_ext, fxy_c_ho);
  gkyl_array_copy(fxy_c, fxy_c_ho);

  gkyl_proj_on_basis_release(projf);

  gkyl_proj_on_basis *projw =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalWeight_2x, NULL);

  struct gkyl_array *wxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *wxy_c_ho = use_gpu ? mkarr(wxy_c->ncomp, wxy_c->size, false) :
                                          gkyl_array_acquire(wxy_c);
  gkyl_proj_on_basis_advance(projw, 0.0, &local_ext, wxy_c_ho);
  gkyl_array_copy(wxy_c, wxy_c_ho);

  gkyl_proj_on_basis_release(projw);

  // create and run the array average updater to average on x only
  int avg_dim_x[] = {1, 0, 0};
  struct gkyl_array_average_inp inp_avg_x = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = basis_y,
    .local = &local,
    .local_avg = &local_y,
    .local_avg_ext = &local_y_ext,
    .weight = wxy_c,
    .avg_dim = avg_dim_x,
    .use_gpu = use_gpu,
  };
  struct gkyl_array *fy_c = mkarr(basis_y.num_basis, local_y_ext.volume, use_gpu);

  struct gkyl_array_average *avg_x = gkyl_array_average_inew(&inp_avg_x);
  gkyl_array_average_advance(
    avg_x, fxy_c,
    fy_c
  ); // fy_c is DG coeff of int[w(x,y) f(x,y)]dx / int[w(x,y)]dx

  gkyl_array_average_release(avg_x);

  // obtain x integral of the weight too
  struct gkyl_array_average_inp inp_int_x = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = basis_y,
    .local = &local,
    .local_avg = &local_y,
    .local_avg_ext = &local_y_ext,
    .weight = NULL,
    .avg_dim = avg_dim_x,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *int_x = gkyl_array_average_inew(&inp_int_x);

  struct gkyl_array *wy_c = mkarr(basis_y.num_basis, local_y_ext.volume, use_gpu);
  gkyl_array_average_advance(int_x, wxy_c, wy_c); // wy_c is DG coeff of int[w(x,y)]dx

  gkyl_array_average_release(int_x);

  // we now remove manually the denominator
  gkyl_dg_mul_op_range(
    &basis_y, 0, fy_c, 0, fy_c, 0, wy_c,
    &local_y
  ); // fy_c is DG coeff of int[w(x,y) f(x,y)]dx

  // average over y now
  int avg_dim_y[] = {1, 0, 0};
  struct gkyl_array_average_inp inp_int_y = {
    .grid = &grid_y,
    .basis = basis_y,
    .basis_avg = red_basis,
    .local = &local_y,
    .local_avg = &red_local,
    .local_avg_ext = &red_local_ext,
    .weight = NULL,
    .avg_dim = avg_dim_y,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *int_y = gkyl_array_average_inew(&inp_int_y);

  struct gkyl_array *intf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(
    int_y, fy_c,
    intf_c
  ); // intf_c is DG coeff of int[int[w(x,y) f(x,y)]dx]dy

  struct gkyl_array *intw_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(int_y, wy_c, intw_c); // intw_c is DG coeff of int[int[w(x,y)]dx]dy

  gkyl_array_average_release(int_y);

  // retrieve the computed average from the device (if applicable)
  double *intf_c0 = gkyl_array_fetch(intf_c, 0);
  double intf_c0_ho[1];
  if (use_gpu) {
    gkyl_cu_memcpy(intf_c0_ho, intf_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(intf_c0_ho, intf_c0, sizeof(double));
  }

  double *intw_c0 = gkyl_array_fetch(intw_c, 0);
  double intw_c0_ho[1];
  if (use_gpu) {
    gkyl_cu_memcpy(intw_c0_ho, intw_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(intw_c0_ho, intw_c0, sizeof(double));
  }

  double integral_wf = intf_c0_ho[0];
  double integral_w = intw_c0_ho[0];
  double result = integral_wf / integral_w;

  double intwf_ref = solution_array_integrate(grid, basis, local_ext, local, wxy_c, fxy_c, use_gpu);
  double intw_ref = solution_array_integrate(grid, basis, local_ext, local, NULL, wxy_c, use_gpu);
  double solution = intwf_ref / intw_ref;

  double rel_err = fabs(result - solution) / fabs(solution);

  // check results two step avg
  TEST_CHECK(gkyl_compare(rel_err, 0, 1e-12));

  // clean up
  gkyl_array_release(fxy_c);
  gkyl_array_release(wxy_c);
  gkyl_array_release(fy_c);
  gkyl_array_release(wy_c);
  gkyl_array_release(intf_c);
  gkyl_array_release(intw_c);
  gkyl_array_release(fxy_c_ho);
  gkyl_array_release(wxy_c_ho);
}

// two steps averaging x,y -> x -> avg
void
test_2x_avgy_avgx(int poly_order, bool use_gpu)
{
  // define grid and basis
  double lower[] = {-4.0, -3.0}, upper[] = {6.0, 5.0};
  int cells[] = {16, 8};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_rect_grid grid_x;
  gkyl_rect_grid_init(&grid_x, 1, &lower[0], &upper[0], &cells[0]);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // reduced grid and basis for two stage reduction
  struct gkyl_range local_x, local_x_ext;
  int ghost_x[] = {ghost[0]};
  gkyl_create_grid_ranges(&grid_x, ghost_x, &local_x_ext, &local_x);

  struct gkyl_basis basis_x;
  gkyl_cart_modal_serendip(&basis_x, 1, poly_order);

  // define the reduced range and basis for scalar result
  struct gkyl_range red_local, red_local_ext;
  gkyl_range_init(&red_local, 1, &local.lower[1], &local.lower[1]);
  gkyl_range_init(&red_local_ext, 1, &local_ext.lower[1], &local_ext.lower[1]);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  // project the target function and weight
  gkyl_proj_on_basis *projf =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalFunc_2x, NULL);

  struct gkyl_array *fxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *fxy_c_ho = use_gpu ? mkarr(fxy_c->ncomp, fxy_c->size, false) :
                                          gkyl_array_acquire(fxy_c);
  gkyl_proj_on_basis_advance(projf, 0.0, &local_ext, fxy_c_ho);
  gkyl_array_copy(fxy_c, fxy_c_ho);

  gkyl_proj_on_basis_release(projf);

  gkyl_proj_on_basis *projw =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalWeight_2x, NULL);

  struct gkyl_array *wxy_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);

  struct gkyl_array *wxy_c_ho = use_gpu ? mkarr(wxy_c->ncomp, wxy_c->size, false) :
                                          gkyl_array_acquire(wxy_c);
  gkyl_proj_on_basis_advance(projw, 0.0, &local_ext, wxy_c_ho);
  gkyl_array_copy(wxy_c, wxy_c_ho);

  gkyl_proj_on_basis_release(projw);

  // create and run the array average updater to average on y only
  int avg_dim_y[] = {0, 1, 0};
  struct gkyl_array_average_inp inp_avg_x = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = basis_x,
    .local = &local,
    .local_avg = &local_x,
    .local_avg_ext = &local_x_ext,
    .weight = wxy_c,
    .avg_dim = avg_dim_y,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *avg_x = gkyl_array_average_inew(&inp_avg_x);

  struct gkyl_array *fx_c = mkarr(basis_x.num_basis, local_x_ext.volume, use_gpu);

  gkyl_array_average_advance(
    avg_x, fxy_c,
    fx_c
  ); // fx_c is DG coeff of int[w(x,y) f(x,y)]dx / int[w(x,y)]dx

  gkyl_array_average_release(avg_x);

  // obtain x integral of the weight too
  struct gkyl_array_average_inp inp_int_x = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = basis_x,
    .local = &local,
    .local_avg = &local_x,
    .local_avg_ext = &local_x_ext,
    .weight = NULL,
    .avg_dim = avg_dim_y,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *int_x = gkyl_array_average_inew(&inp_int_x);

  struct gkyl_array *wx_c = mkarr(basis_x.num_basis, local_x_ext.volume, use_gpu);
  gkyl_array_average_advance(int_x, wxy_c, wx_c); // wx_c is DG coeff of int[w(x,y)]dx

  gkyl_array_average_release(int_x);

  // we now remove manually the denominator
  gkyl_dg_mul_op_range(
    &basis_x, 0, fx_c, 0, fx_c, 0, wx_c,
    &local_x
  ); // fx_c is DG coeff of int[w(x,y) f(x,y)]dx

  // create and run the array average updater to integrate on y
  int avg_dim_x[] = {1, 0, 0};
  struct gkyl_array_average_inp inp_int_y = {
    .grid = &grid_x,
    .basis = basis_x,
    .basis_avg = red_basis,
    .local = &local_x,
    .local_avg = &red_local,
    .local_avg_ext = &red_local_ext,
    .weight = NULL,
    .avg_dim = avg_dim_x,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *int_y = gkyl_array_average_inew(&inp_int_y);

  struct gkyl_array *intf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(
    int_y, fx_c,
    intf_c
  ); // intf_c is DG coeff of int[int[w(x,y) f(x,y)]dx]dy

  // obtain full integral of weight too
  struct gkyl_array *intw_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(int_y, wx_c, intw_c); // intw_c is DG coeff of int[int[w(x,y)]dx]dy

  gkyl_array_average_release(int_y);

  // check results two step avg
  double intf_c0_ho[1];
  double *intf_c0 = gkyl_array_fetch(intf_c, 0);
  if (use_gpu) {
    gkyl_cu_memcpy(intf_c0_ho, intf_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(intf_c0_ho, intf_c0, sizeof(double));
  }

  double intw_c0_ho[1];
  double *intw_c0 = gkyl_array_fetch(intw_c, 0);
  if (use_gpu) {
    gkyl_cu_memcpy(intw_c0_ho, intw_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(intw_c0_ho, intw_c0, sizeof(double));
  }

  double integral_wf = intf_c0_ho[0];
  double integral_w = intw_c0_ho[0];
  double result = integral_wf / integral_w;

  double intwf_ref = solution_array_integrate(grid, basis, local_ext, local, wxy_c, fxy_c, use_gpu);
  double intw_ref = solution_array_integrate(grid, basis, local_ext, local, NULL, wxy_c, use_gpu);
  double solution = intwf_ref / intw_ref;

  double rel_err = fabs(result - solution) / fabs(solution);

  TEST_CHECK(gkyl_compare(rel_err, 0.0, 1e-12));

  // clean up
  gkyl_array_release(fxy_c);
  gkyl_array_release(wxy_c);
  gkyl_array_release(fx_c);
  gkyl_array_release(wx_c);
  gkyl_array_release(intf_c);
  gkyl_array_release(intw_c);
  gkyl_array_release(wxy_c_ho);
  gkyl_array_release(fxy_c_ho);
}

// test 3x
void
evalFunc_3x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  double lower[] = {-4., -3., -2.}, upper[] = {6., 5., 4.};
  double Lx = upper[0] - lower[0];
  double Ly = upper[1] - lower[1];
  double Lz = upper[2] - lower[2];
  double k_x = 2. * M_PI / Lx;
  double k_y = 2. * M_PI / Ly;
  double k_z = 2. * M_PI / Lz;
  double phi = 0.5;
  fout[0] = x * y * z * sin(1.5 * k_x * x + 0.75 * k_y * y + 0.5 * k_z * z + phi) *
            cos(1.42 * k_y * y) * cos(4.20 * k_z * z);
}
void
evalWeight_3x(double t, const double *xn, double *restrict fout, void *ctx)
{
  double x = xn[0];
  double y = xn[1];
  double z = xn[2];
  fout[0] = 1 + x * x + y * y + z * z;
}
// two steps average x,y,z -> y,z -> avg
void
test_3x_avgx_avgyz(int poly_order, bool use_gpu)
{
  // define grids and basis
  double lower[] = {-4.0, -3.0, -2.0}, upper[] = {6.0, 5.0, 4.0};
  // int cells[] = {64, 48, 32};
  int cells[] = {16, 8, 4};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1, 1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // 2dim subgrid
  struct gkyl_rect_grid grid_yz;
  double yz_grid_lower[] = {grid.lower[1], grid.lower[2]};
  double yz_grid_upper[] = {grid.upper[1], grid.upper[2]};
  int yz_grid_cells[] = {grid.cells[1], grid.cells[2]};
  gkyl_rect_grid_init(&grid_yz, 2, yz_grid_lower, yz_grid_upper, yz_grid_cells);

  struct gkyl_range local_yz, local_yz_ext;
  int ghost_yz[] = {ghost[1], ghost[2]};
  gkyl_create_grid_ranges(&grid_yz, ghost_yz, &local_yz_ext, &local_yz);

  struct gkyl_basis basis_yz;
  gkyl_cart_modal_serendip(&basis_yz, 2, poly_order);

  // define the reduced range and basis for full reduction
  struct gkyl_range red_local, red_local_ext;
  gkyl_range_init(&red_local, 1, &local.lower[0], &local.lower[0]);
  gkyl_range_init(&red_local_ext, 1, &local_ext.lower[0], &local_ext.lower[0]);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  // project the target function and weight
  gkyl_proj_on_basis *projf =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalFunc_3x, NULL);

  struct gkyl_array *fxyz_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *fxyz_c_ho = use_gpu ? mkarr(fxyz_c->ncomp, fxyz_c->size, false) :
                                           gkyl_array_acquire(fxyz_c);
  gkyl_proj_on_basis_advance(projf, 0.0, &local_ext, fxyz_c_ho);
  gkyl_array_copy(fxyz_c, fxyz_c_ho);

  gkyl_proj_on_basis_release(projf);

  gkyl_proj_on_basis *projw =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalWeight_3x, NULL);

  struct gkyl_array *wxyz_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *wxyz_c_ho = use_gpu ? mkarr(wxyz_c->ncomp, wxyz_c->size, false) :
                                           gkyl_array_acquire(wxyz_c);
  gkyl_proj_on_basis_advance(projw, 0.0, &local_ext, wxyz_c_ho);
  gkyl_array_copy(wxyz_c, wxyz_c_ho);

  gkyl_proj_on_basis_release(projw);

  // create and run the array average updater to average on x only
  int avg_dim_x[] = {1, 0, 0};
  struct gkyl_array_average_inp inp_avg_xyz_to_yz = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = basis_yz,
    .local = &local,
    .local_avg = &local_yz,
    .local_avg_ext = &local_yz_ext,
    .weight = wxyz_c,
    .avg_dim = avg_dim_x,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *avg_xyz_to_yz = gkyl_array_average_inew(&inp_avg_xyz_to_yz);

  struct gkyl_array *fyz_c = mkarr(basis_yz.num_basis, local_yz_ext.volume, use_gpu);
  gkyl_array_average_advance(
    avg_xyz_to_yz, fxyz_c, fyz_c
  ); // fy_c is DG coeff of int[w(x,y,z) f(x,y,z)]dx / int[w(x,y,z)]dx
  gkyl_array_average_release(avg_xyz_to_yz);

  // obtain x integral of the weight too
  struct gkyl_array_average_inp inp_int_xyz_to_yz = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = basis_yz,
    .local = &local,
    .local_avg = &local_yz,
    .local_avg_ext = &local_yz_ext,
    .weight = NULL,
    .avg_dim = avg_dim_x,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *int_xyz_to_yz = gkyl_array_average_inew(&inp_int_xyz_to_yz);

  struct gkyl_array *wyz_c = mkarr(basis_yz.num_basis, local_yz_ext.volume, use_gpu);
  gkyl_array_average_advance(int_xyz_to_yz, wxyz_c, wyz_c); // wy_c is DG coeff of int[w(x,y)]dy
  gkyl_array_average_release(int_xyz_to_yz);

  // we now remove manually the denominator
  gkyl_dg_mul_op_range(
    &basis_yz, 0, fyz_c, 0, fyz_c, 0, wyz_c,
    &local_yz
  ); // fy_c is DG coeff of int[w(x,y) f(x,y)]dy

  // create and run the array average updater to average on y and z (first second dim)
  int avg_dim_yz[] = {1, 1, 0};
  struct gkyl_array_average_inp inp_int_yz = {
    .grid = &grid_yz,
    .basis = basis_yz,
    .basis_avg = red_basis,
    .local = &local_yz,
    .local_avg = &red_local,
    .local_avg_ext = &red_local_ext,
    .weight = NULL,
    .avg_dim = avg_dim_yz,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *int_yz = gkyl_array_average_inew(&inp_int_yz);

  struct gkyl_array *intf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(
    int_yz, fyz_c,
    intf_c
  ); // intf_c is DG coeff of int[int[w(x,y) f(x,y)]dy]dx

  // obtain full integral of weight too
  struct gkyl_array *intw_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(int_yz, wyz_c, intw_c); // intw_c is DG coeff of int[int[w(x,y)]dy]dx

  gkyl_array_average_release(int_yz);

  // check results
  double *intf_c0 = gkyl_array_fetch(intf_c, 0);
  double *intf_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu) {
    gkyl_cu_memcpy(intf_c0_ho, intf_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(intf_c0_ho, intf_c0, sizeof(double));
  }

  double *intw_c0 = gkyl_array_fetch(intw_c, 0);
  double *intw_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu) {
    gkyl_cu_memcpy(intw_c0_ho, intw_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(intw_c0_ho, intw_c0, sizeof(double));
  }

  double integral_wf = intf_c0_ho[0];
  double integral_w = intw_c0_ho[0];
  double result = integral_wf / integral_w;

  double intwf_ref =
    solution_array_integrate(grid, basis, local_ext, local, wxyz_c, fxyz_c, use_gpu);
  double intw_ref = solution_array_integrate(grid, basis, local_ext, local, NULL, wxyz_c, use_gpu);
  double solution = intwf_ref / intw_ref;

  double rel_err = fabs(result - solution) / fabs(solution);

  TEST_CHECK(gkyl_compare(rel_err, 0, 1e-12));

  // clean up
  gkyl_array_release(fxyz_c);
  gkyl_array_release(wxyz_c);
  gkyl_array_release(fyz_c);
  gkyl_array_release(wyz_c);
  gkyl_array_release(intf_c);
  gkyl_array_release(intw_c);
  gkyl_free(intf_c0_ho);
  gkyl_free(intw_c0_ho);
  gkyl_array_release(fxyz_c_ho);
  gkyl_array_release(wxyz_c_ho);
}
// two steps average x,y,z -> x -> avg
void
test_3x_avgyz_avgx(int poly_order, bool use_gpu)
{
  // define grids and basis
  double lower[] = {-4.0, -3.0, -2.0}, upper[] = {6.0, 5.0, 4.0};
  int cells[] = {16, 8, 4};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int ghost[] = {1, 1, 1};
  struct gkyl_range local, local_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // 2dim subgrid
  struct gkyl_rect_grid grid_x;
  gkyl_rect_grid_init(&grid_x, 1, &grid.lower[0], &grid.upper[0], &grid.cells[0]);

  // define the reduced range and basis for averaging along x only
  struct gkyl_range local_x, local_x_ext;
  int ghost_x[] = {ghost[0]};
  gkyl_create_grid_ranges(&grid_x, ghost_x, &local_x_ext, &local_x);

  struct gkyl_basis basis_x;
  gkyl_cart_modal_serendip(&basis_x, 1, poly_order);

  // define the reduced range and basis for full reduction
  struct gkyl_range red_local, red_local_ext;
  gkyl_range_init(&red_local, 1, &local.lower[0], &local.lower[0]);
  gkyl_range_init(&red_local_ext, 1, &local_ext.lower[0], &local_ext.lower[0]);

  struct gkyl_basis red_basis;
  gkyl_cart_modal_serendip(&red_basis, 1, poly_order);

  // project the target function and weight
  gkyl_proj_on_basis *projf =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalFunc_3x, NULL);

  struct gkyl_array *fxyz_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *fxyz_c_ho = use_gpu ? mkarr(fxyz_c->ncomp, fxyz_c->size, false) :
                                           gkyl_array_acquire(fxyz_c);
  gkyl_proj_on_basis_advance(projf, 0.0, &local_ext, fxyz_c_ho);
  gkyl_array_copy(fxyz_c, fxyz_c_ho);

  gkyl_proj_on_basis_release(projf);

  gkyl_proj_on_basis *projw =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalWeight_3x, NULL);

  struct gkyl_array *wxyz_c = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *wxyz_c_ho = use_gpu ? mkarr(wxyz_c->ncomp, wxyz_c->size, false) :
                                           gkyl_array_acquire(wxyz_c);
  gkyl_proj_on_basis_advance(projw, 0.0, &local_ext, wxyz_c_ho);
  gkyl_array_copy(wxyz_c, wxyz_c_ho);

  gkyl_proj_on_basis_release(projw);

  // create and run the array average updater to average y and z
  int avg_dim_yz[] = {0, 1, 1};
  struct gkyl_array_average_inp inp_avg_xyz_to_x = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = basis_x,
    .local = &local,
    .local_avg = &local_x,
    .local_avg_ext = &local_x_ext,
    .weight = wxyz_c,
    .avg_dim = avg_dim_yz,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *avg_xyz_to_x = gkyl_array_average_inew(&inp_avg_xyz_to_x);

  struct gkyl_array *fx_c = mkarr(basis_x.num_basis, local_x_ext.volume, use_gpu);
  gkyl_array_average_advance(avg_xyz_to_x, fxyz_c, fx_c); //

  gkyl_array_average_release(avg_xyz_to_x);

  // obtain x integral of the weight too
  struct gkyl_array_average_inp inp_int_xyz_to_x = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = basis_x,
    .local = &local,
    .local_avg = &local_x,
    .local_avg_ext = &local_x_ext,
    .weight = NULL,
    .avg_dim = avg_dim_yz,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *int_xyz_to_x = gkyl_array_average_inew(&inp_int_xyz_to_x);

  struct gkyl_array *wx_c = mkarr(basis_x.num_basis, local_x_ext.volume, use_gpu);
  gkyl_array_average_advance(int_xyz_to_x, wxyz_c, wx_c); // wy_c is DG coeff of int[w(x,y)]dydz

  gkyl_array_average_release(int_xyz_to_x);

  // remove manually the denominator
  gkyl_dg_mul_op_range(
    &basis_x, 0, fx_c, 0, fx_c, 0, wx_c,
    &local_x
  ); // fy_c is DG coeff of int[w(x,y) f(x,y)]dy

  // create and run the array average updater to average on x
  int avg_dim_x[] = {1, 0, 0};
  struct gkyl_array_average_inp inp_int_x = {
    .grid = &grid_x,
    .basis = basis_x,
    .basis_avg = red_basis,
    .local = &local_x,
    .local_avg = &red_local,
    .local_avg_ext = &red_local_ext,
    .weight = NULL,
    .avg_dim = avg_dim_x,
    .use_gpu = use_gpu,
  };
  struct gkyl_array_average *int_x = gkyl_array_average_inew(&inp_int_x);

  struct gkyl_array *intf_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(
    int_x, fx_c,
    intf_c
  ); // intf_c is DG coeff of int[int[int[w(x,y) f(x,y)]dy]dz]dx

  // obtain full integral of weight too
  struct gkyl_array *intw_c = mkarr(red_basis.num_basis, red_local.volume, use_gpu);
  gkyl_array_average_advance(
    int_x, wx_c,
    intw_c
  ); // intw_c is DG coeff of int[int[int[w(x,y)]dy]dz]dx

  gkyl_array_average_release(int_x);

  // check results two step avg
  double *intf_c0 = gkyl_array_fetch(intf_c, 0);
  double *intf_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu) {
    gkyl_cu_memcpy(intf_c0_ho, intf_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(intf_c0_ho, intf_c0, sizeof(double));
  }

  double *intw_c0 = gkyl_array_fetch(intw_c, 0);
  double *intw_c0_ho = gkyl_malloc(sizeof(double));
  if (use_gpu) {
    gkyl_cu_memcpy(intw_c0_ho, intw_c0, sizeof(double), GKYL_CU_MEMCPY_D2H);
  } else {
    memcpy(intw_c0_ho, intw_c0, sizeof(double));
  }

  double integral_wf = intf_c0_ho[0];
  double integral_w = intw_c0_ho[0];
  double result = integral_wf / integral_w;

  double intwf_ref =
    solution_array_integrate(grid, basis, local_ext, local, wxyz_c, fxyz_c, use_gpu);
  double intw_ref = solution_array_integrate(grid, basis, local_ext, local, NULL, wxyz_c, use_gpu);
  double solution = intwf_ref / intw_ref;

  double rel_err = fabs(result - solution) / fabs(solution);

  TEST_CHECK(gkyl_compare(rel_err, 0, 1e-12));

  // clean up
  gkyl_array_release(fxyz_c);
  gkyl_array_release(wxyz_c);
  gkyl_array_release(fx_c);
  gkyl_array_release(wx_c);
  gkyl_array_release(intf_c);
  gkyl_array_release(intw_c);
  gkyl_free(intf_c0_ho);
  gkyl_free(intw_c0_ho);
  gkyl_array_release(fxyz_c_ho);
  gkyl_array_release(wxyz_c_ho);
}

// Average over y with the output in the 2x basis, compared with the reduced-dim output.
void
test_2x_fulldim(int poly_order, bool use_gpu)
{
  double lower[] = {-4.0, -3.0}, upper[] = {6.0, 5.0};
  int cells[] = {16, 8};
  int ndim = sizeof(lower) / sizeof(lower[0]);

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  struct gkyl_rect_grid grid_x;
  gkyl_rect_grid_init(&grid_x, 1, &lower[0], &upper[0], &cells[0]);

  struct gkyl_basis basis, basis_x;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);
  gkyl_cart_modal_serendip(&basis_x, 1, poly_order);

  int ghost[] = {1, 1};
  struct gkyl_range local, local_ext, local_x, local_x_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  gkyl_create_grid_ranges(&grid_x, &ghost[0], &local_x_ext, &local_x);

  gkyl_proj_on_basis *projf =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalFunc_2x, NULL);
  gkyl_proj_on_basis *projw =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, evalWeight_2x, NULL);
  struct gkyl_array *f_ho = mkarr(basis.num_basis, local_ext.volume, false);
  struct gkyl_array *w_ho = mkarr(basis.num_basis, local_ext.volume, false);
  gkyl_proj_on_basis_advance(projf, 0.0, &local_ext, f_ho);
  gkyl_proj_on_basis_advance(projw, 0.0, &local_ext, w_ho);
  gkyl_proj_on_basis_release(projf);
  gkyl_proj_on_basis_release(projw);
  struct gkyl_array *f = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  struct gkyl_array *w = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_array_copy(f, f_ho);
  gkyl_array_copy(w, w_ho);

  int avg_dim[] = {0, 1};
  struct gkyl_array *avg_ho = mkarr(basis_x.num_basis, local_x_ext.volume, false);
  struct gkyl_array *favg_ho = mkarr(basis.num_basis, local_ext.volume, false);
  struct gkyl_array *avg = mkarr(basis_x.num_basis, local_x_ext.volume, use_gpu);
  struct gkyl_array *favg = mkarr(basis.num_basis, local_ext.volume, use_gpu);

  for (int weighted = 0; weighted < 2; weighted++) {
    struct gkyl_array_average_inp inp = {
      .grid = &grid,
      .basis = basis,
      .basis_avg = basis_x,
      .local = &local,
      .local_avg = &local_x,
      .local_avg_ext = &local_x_ext,
      .weight = weighted ? w : NULL,
      .avg_dim = avg_dim,
      .use_gpu = use_gpu,
    };
    struct gkyl_array_average *up_red = gkyl_array_average_inew(&inp);
    gkyl_array_average_advance(up_red, f, avg);
    gkyl_array_average_release(up_red);

    inp.full_dim_out = true;
    inp.local_avg_ext = &local_ext;
    struct gkyl_array_average *up_full = gkyl_array_average_inew(&inp);
    gkyl_array_average_advance(up_full, f, favg);
    gkyl_array_average_release(up_full);

    gkyl_array_copy(avg_ho, avg);
    gkyl_array_copy(favg_ho, favg);

    // The full-dim output evaluated anywhere in y must match the reduced output.
    double zs[][2] = {{-0.7, 0.3}, {0.2, -0.9}, {0.95, 0.85}};
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &local);
    while (gkyl_range_iter_next(&iter)) {
      const double *favg_c = gkyl_array_cfetch(favg_ho, gkyl_range_idx(&local, iter.idx));
      const double *avg_c = gkyl_array_cfetch(avg_ho, gkyl_range_idx(&local_x, &iter.idx[0]));
      for (int i = 0; i < 3; i++) {
        double ref = basis_x.eval_expand(&zs[i][0], avg_c);
        TEST_CHECK(gkyl_compare(basis.eval_expand(zs[i], favg_c), ref, 1e-12 * fabs(ref)));
      }
    }
  }

  gkyl_array_release(f);
  gkyl_array_release(w);
  gkyl_array_release(f_ho);
  gkyl_array_release(w_ho);
  gkyl_array_release(avg);
  gkyl_array_release(favg);
  gkyl_array_release(avg_ho);
  gkyl_array_release(favg_ho);
}

struct ctx_3x2v {
  int ndim; // Dimensionality of the function.
  int avg_dir; // Direction averaged over.
  bool averaged; // Whether to evaluate the average (on the reduced coordinates).
};

// a(s)*g(other coordinates), with a quadratic and g a product of Gaussians.
void
evalFunc_3x2v(double t, const double *xn, double *restrict fout, void *ctx)
{
  struct ctx_3x2v *pctx = ctx;
  if (pctx->averaged) {
    fout[0] = 1.0 + 0.1 / 3.0;
    for (int r = 0; r < pctx->ndim; r++) {
      int d = r < pctx->avg_dir ? r : r + 1;
      fout[0] *= exp(-(0.5 + 0.1 * d) * xn[r] * xn[r]);
    }
  } else {
    double s = xn[pctx->avg_dir];
    fout[0] = 1.0 + 0.3 * s + 0.1 * s * s;
    for (int d = 0; d < pctx->ndim; d++) {
      if (d != pctx->avg_dir) {
        fout[0] *= exp(-(0.5 + 0.1 * d) * xn[d] * xn[d]);
      }
    }
  }
}

// Basis of a cdim x vdim phase-space field, and of the field left after
// averaging over avg_dir.
static void
phase_bases(
  int cdim, int vdim, enum gkyl_basis_type b_type, int poly_order, int avg_dir,
  struct gkyl_basis *basis, struct gkyl_basis *basis_red
)
{
  int ndim = cdim + vdim;
  int cdim_red = avg_dir < cdim ? cdim - 1 : cdim;
  int vdim_red = avg_dir < cdim ? vdim : vdim - 1;
  if (b_type == GKYL_BASIS_MODAL_GKHYBRID) {
    gkyl_cart_modal_gkhybrid(basis, cdim, vdim);
    if (vdim_red == 0) {
      gkyl_cart_modal_serendip(basis_red, cdim_red, 1);
    } else if (cdim_red == 0) {
      if (avg_dir == cdim) {
        gkyl_cart_modal_serendip(basis_red, 1, 1); // mu only.
      } else {
        gkyl_cart_modal_gkhybrid_vel(basis_red, vdim_red);
      }
    } else if (avg_dir == cdim) {
      gkyl_cart_modal_serendip(basis_red, cdim_red + 1, 1); // vpar averaged.
    } else if (avg_dir == cdim + 1) {
      gkyl_cart_modal_hybrid(basis_red, cdim_red, 1); // mu averaged.
    } else {
      gkyl_cart_modal_gkhybrid(basis_red, cdim_red, vdim_red);
    }
  } else if (b_type == GKYL_BASIS_MODAL_HYBRID) {
    gkyl_cart_modal_hybrid(basis, cdim, vdim);
    if (vdim_red == 0) {
      gkyl_cart_modal_serendip(basis_red, cdim_red, 1);
    } else if (cdim_red == 0) {
      gkyl_cart_modal_serendip(basis_red, vdim_red, 2);
    } else {
      gkyl_cart_modal_hybrid(basis_red, cdim_red, vdim_red);
    }
  } else {
    gkyl_cart_modal_serendip(basis, ndim, poly_order);
    gkyl_cart_modal_serendip(basis_red, ndim - 1, poly_order);
  }
}

// Average a phase-space array over the dimensions flagged in avg_dim. The
// full-dim output is compared with the mean of the cells along those
// dimensions and, for a single averaged dimension, the reduced output with
// the projection of the analytic average.
void
test_3x2v(
  int cdim, int vdim, enum gkyl_basis_type b_type, int poly_order, const int *avg_dim, bool use_gpu
)
{
  int ndim = cdim + vdim, ndim_red = ndim - 1, avg_dir = -1, num_avg_dim = 0;
  for (int d = 0; d < ndim; d++) {
    if (avg_dim[d]) {
      avg_dir = d;
      num_avg_dim++;
    }
  }
  bool single_dir = num_avg_dim == 1;

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM], lower_red[GKYL_MAX_DIM], upper_red[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM], cells_red[GKYL_MAX_DIM], ghost_red[GKYL_MAX_DIM];
  for (int d = 0, r = 0; d < ndim; d++) {
    lower[d] = -1.0;
    upper[d] = 1.0;
    cells[d] = avg_dim[d] ? 6 : 4;
    ghost[d] = 0;
    if (d != avg_dir) {
      lower_red[r] = lower[d];
      upper_red[r] = upper[d];
      cells_red[r] = cells[d];
      ghost_red[r] = 0;
      r++;
    }
  }

  struct gkyl_rect_grid grid, grid_red;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);
  gkyl_rect_grid_init(&grid_red, ndim_red, lower_red, upper_red, cells_red);

  struct gkyl_basis basis, basis_red;
  phase_bases(cdim, vdim, b_type, poly_order, avg_dir, &basis, &basis_red);

  struct gkyl_range local, local_ext, local_red, local_red_ext;
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);
  gkyl_create_grid_ranges(&grid_red, ghost_red, &local_red_ext, &local_red);

  struct ctx_3x2v ctx = {.ndim = ndim, .avg_dir = avg_dir, .averaged = false};
  struct ctx_3x2v ctx_red = {.ndim = ndim_red, .avg_dir = avg_dir, .averaged = true};
  gkyl_proj_on_basis *projf = gkyl_proj_on_basis_new(&grid, &basis, 3, 1, evalFunc_3x2v, &ctx);
  gkyl_proj_on_basis *projf_red =
    gkyl_proj_on_basis_new(&grid_red, &basis_red, 3, 1, evalFunc_3x2v, &ctx_red);
  struct gkyl_array *f_ho = mkarr(basis.num_basis, local_ext.volume, false);
  struct gkyl_array *ref_ho = mkarr(basis_red.num_basis, local_red_ext.volume, false);
  gkyl_proj_on_basis_advance(projf, 0.0, &local, f_ho);
  gkyl_proj_on_basis_advance(projf_red, 0.0, &local_red, ref_ho);
  gkyl_proj_on_basis_release(projf);
  gkyl_proj_on_basis_release(projf_red);
  struct gkyl_array *f = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_array_copy(f, f_ho);

  struct gkyl_array_average_inp inp = {
    .grid = &grid,
    .basis = basis,
    .basis_avg = basis_red,
    .vdim = vdim,
    .local = &local,
    .local_avg = &local_red,
    .local_avg_ext = &local_red_ext,
    .avg_dim = avg_dim,
    .use_gpu = use_gpu,
  };
  struct gkyl_array *avg = mkarr(basis_red.num_basis, local_red_ext.volume, use_gpu);
  struct gkyl_array *avg_ho = mkarr(basis_red.num_basis, local_red_ext.volume, false);
  if (single_dir) {
    struct gkyl_array_average *up_red = gkyl_array_average_inew(&inp);
    gkyl_array_average_advance(up_red, f, avg);
    gkyl_array_average_release(up_red);
    gkyl_array_copy(avg_ho, avg);
  }

  inp.full_dim_out = true;
  inp.local_avg_ext = &local_ext;
  struct gkyl_array_average *up_full = gkyl_array_average_inew(&inp);
  struct gkyl_array *favg = mkarr(basis.num_basis, local_ext.volume, use_gpu);
  gkyl_array_average_advance(up_full, f, favg);
  gkyl_array_average_release(up_full);
  struct gkyl_array *favg_ho = mkarr(basis.num_basis, local_ext.volume, false);
  gkyl_array_copy(favg_ho, favg);

  double f_max = 0.0;
  for (long i = 0; i < f_ho->size * f_ho->ncomp; i++) {
    f_max = fmax(f_max, fabs(((const double *)f_ho->data)[i]));
  }

  struct gkyl_range_iter iter;
  if (single_dir) {
    gkyl_range_iter_init(&iter, &local_red);
    while (gkyl_range_iter_next(&iter)) {
      long lidx = gkyl_range_idx(&local_red, iter.idx);
      const double *avg_c = gkyl_array_cfetch(avg_ho, lidx);
      const double *ref_c = gkyl_array_cfetch(ref_ho, lidx);
      for (int k = 0; k < basis_red.num_basis; k++) {
        TEST_CHECK(fabs(avg_c[k] - ref_c[k]) < 1e-11 * f_max);
        TEST_MSG(
          "cdim=%d vdim=%d b_type=%d avg_dir=%d k=%d: %g vs %g", cdim, vdim, b_type, avg_dir, k,
          avg_c[k], ref_c[k]
        );
      }
    }
  }

  // Basis functions independent of the averaged dimensions.
  int indep[160], dim_remains[GKYL_MAX_DIM];
  double z0[GKYL_MAX_DIM], z1[GKYL_MAX_DIM], b0[160], b1[160];
  for (int d = 0; d < ndim; d++) {
    dim_remains[d] = 1 - avg_dim[d];
    z0[d] = z1[d] = 0.17 + 0.11 * d;
    if (avg_dim[d]) {
      z0[d] = -0.6 + 0.07 * d;
      z1[d] = 0.35 - 0.05 * d;
    }
  }
  basis.eval(z0, b0);
  basis.eval(z1, b1);
  for (int k = 0; k < basis.num_basis; k++) {
    indep[k] = fabs(b0[k] - b1[k]) < 1e-12;
  }

  int num_cells_avg = 1;
  for (int d = 0; d < ndim; d++) {
    num_cells_avg *= avg_dim[d] ? cells[d] : 1;
  }
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    const double *favg_c = gkyl_array_cfetch(favg_ho, gkyl_range_idx(&local, iter.idx));
    struct gkyl_range rng_avg;
    struct gkyl_range_iter iter_avg;
    gkyl_range_deflate(&rng_avg, &local, dim_remains, iter.idx);
    for (int k = 0; k < basis.num_basis; k++) {
      double mean = 0.0;
      gkyl_range_iter_no_split_init(&iter_avg, &rng_avg);
      while (gkyl_range_iter_next(&iter_avg)) {
        const double *f_c = gkyl_array_cfetch(f_ho, gkyl_range_idx(&rng_avg, iter_avg.idx));
        mean += indep[k] ? f_c[k] / num_cells_avg : 0.0;
      }
      TEST_CHECK(fabs(favg_c[k] - mean) < 1e-12 * f_max);
      TEST_MSG("cdim=%d vdim=%d b_type=%d k=%d: %g vs %g", cdim, vdim, b_type, k, favg_c[k], mean);
    }
  }

  gkyl_array_release(f);
  gkyl_array_release(f_ho);
  gkyl_array_release(ref_ho);
  gkyl_array_release(avg);
  gkyl_array_release(favg);
  gkyl_array_release(avg_ho);
  gkyl_array_release(favg_ho);
}

void
test_3x2v_all(bool use_gpu)
{
  enum gkyl_basis_type gkhyb = GKYL_BASIS_MODAL_GKHYBRID, hyb = GKYL_BASIS_MODAL_HYBRID,
                       ser = GKYL_BASIS_MODAL_SERENDIPITY;
  // Single averaged dimension (conf or velocity).
  test_3x2v(1, 1, gkhyb, 1, (int[]){1, 0}, use_gpu);
  test_3x2v(1, 1, gkhyb, 1, (int[]){0, 1}, use_gpu);
  test_3x2v(1, 2, gkhyb, 1, (int[]){1, 0, 0}, use_gpu);
  test_3x2v(1, 2, gkhyb, 1, (int[]){0, 1, 0}, use_gpu);
  test_3x2v(1, 2, gkhyb, 1, (int[]){0, 0, 1}, use_gpu);
  test_3x2v(2, 2, gkhyb, 1, (int[]){0, 1, 0, 0}, use_gpu);
  test_3x2v(2, 2, gkhyb, 1, (int[]){0, 0, 0, 1}, use_gpu);
  test_3x2v(3, 2, gkhyb, 1, (int[]){0, 1, 0, 0, 0}, use_gpu);
  test_3x2v(3, 2, gkhyb, 1, (int[]){0, 0, 0, 1, 0}, use_gpu);
  test_3x2v(1, 2, hyb, 1, (int[]){1, 0, 0}, use_gpu);
  test_3x2v(1, 2, hyb, 1, (int[]){0, 1, 0}, use_gpu);
  test_3x2v(1, 3, hyb, 1, (int[]){0, 0, 0, 1}, use_gpu);
  test_3x2v(2, 2, hyb, 1, (int[]){1, 0, 0, 0}, use_gpu);
  test_3x2v(2, 2, hyb, 1, (int[]){0, 0, 0, 1}, use_gpu);
  test_3x2v(2, 3, hyb, 1, (int[]){0, 1, 0, 0, 0}, use_gpu);
  test_3x2v(2, 3, hyb, 1, (int[]){0, 0, 1, 0, 0}, use_gpu);
  test_3x2v(1, 2, ser, 2, (int[]){1, 0, 0}, use_gpu);
  test_3x2v(1, 3, ser, 2, (int[]){0, 0, 0, 1}, use_gpu);
  test_3x2v(2, 2, ser, 2, (int[]){0, 1, 0, 0}, use_gpu);
  test_3x2v(2, 2, ser, 2, (int[]){0, 0, 1, 0}, use_gpu);
  // Several averaged dimensions (full-dim output only).
  test_3x2v(1, 2, gkhyb, 1, (int[]){1, 1, 0}, use_gpu);
  test_3x2v(2, 2, gkhyb, 1, (int[]){1, 0, 1, 1}, use_gpu);
  test_3x2v(3, 2, gkhyb, 1, (int[]){0, 1, 0, 1, 1}, use_gpu);
  test_3x2v(2, 2, hyb, 1, (int[]){0, 1, 1, 0}, use_gpu);
  test_3x2v(2, 3, hyb, 1, (int[]){1, 0, 0, 1, 1}, use_gpu);
  test_3x2v(2, 2, ser, 2, (int[]){1, 1, 0, 1}, use_gpu);
}

void
test_array_average_1x_ho()
{
  for (int p = 1; p <= 2; p++) {
    test_1x(p, false);
  }
}

void
test_array_average_2x_ho()
{
  for (int p = 1; p <= 2; p++) {
    test_2x_1step(p, false);
    test_2x_intx_inty(p, false);
    test_2x_avgx_avgy(p, false);
    test_2x_avgy_avgx(p, false);
    test_2x_fulldim(p, false);
  }
}

void
test_array_average_3x2v_ho()
{
  test_3x2v_all(false);
}

void
test_array_average_3x_ho()
{
  for (int p = 1; p <= 2; p++) {
    test_3x_avgx_avgyz(p, false);
    test_3x_avgyz_avgx(p, false);
  }
}

#ifdef GKYL_HAVE_CUDA
void
test_array_average_1x_dev()
{
  for (int p = 1; p <= 2; p++) {
    test_1x(p, true);
  }
}

void
test_array_average_2x_dev()
{
  for (int p = 1; p <= 2; p++) {
    test_2x_1step(p, true);
    test_2x_intx_inty(p, true);
    test_2x_avgx_avgy(p, true);
    test_2x_avgy_avgx(p, true);
    test_2x_fulldim(p, true);
  }
}

void
test_array_average_3x2v_dev()
{
  test_3x2v_all(true);
}

void
test_array_average_3x_dev()
{
  for (int p = 1; p <= 2; p++) {
    test_3x_avgx_avgyz(p, true);
    test_3x_avgyz_avgx(p, true);
  }
}

#endif

TEST_LIST = {
  {"test_array_average_1x_ho", test_array_average_1x_ho},
  {"test_array_average_2x_ho", test_array_average_2x_ho},
  {"test_array_average_3x_ho", test_array_average_3x_ho},
  {"test_array_average_3x2v_ho", test_array_average_3x2v_ho},
#ifdef GKYL_HAVE_CUDA
  {"test_array_average_1x_dev", test_array_average_1x_dev},
  {"test_array_average_2x_dev", test_array_average_2x_dev},
  {"test_array_average_3x_dev", test_array_average_3x_dev},
  {"test_array_average_3x2v_dev", test_array_average_3x2v_dev},
#endif
  {NULL, NULL}
};
