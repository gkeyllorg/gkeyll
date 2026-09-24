#include <acutest.h>

#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array_rio.h>
#include <math.h>

#include <assert.h>

void evalFunc_1x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = x*x;
}

void evalFunc_1x_trig(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = sin((2.*M_PI/(4.))*x);
}

void test_1x(int poly_order, int test_func_op)
{
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_eval_on_nodes *evup;
  if (test_func_op==0)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_1x_quad, NULL);
  else if (test_func_op==1)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_1x_trig, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_eval_on_nodes_advance(evup, 0.0, &arr_range, distf);
  gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, "ctest_eval_on_nodes_distf.gkyl");

  double *dfl = gkyl_array_fetch(distf, 0);  // left cell
  double *dfr = gkyl_array_fetch(distf, 1);  // right cell
  if (poly_order == 1) {
    if (test_func_op==0) {
      TEST_CHECK( gkyl_compare( 2.82842712474619e+00, dfl[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.63299316185545e+00, dfl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.82842712474619e+00, dfr[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.63299316185545e+00, dfr[1], 1e-12) );
    } else if (test_func_op==1) {
      TEST_CHECK( gkyl_compare(-8.65956056235493e-17, dfl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.99959962173949e-17, dfl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 8.65956056235493e-17, dfr[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.99959962173949e-17, dfr[1], 1e-12) );
    }
  } else if (poly_order == 2) {
    if (test_func_op==0) {
      TEST_CHECK( gkyl_compare( 1.88561808316413e+00, dfl[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.63299316185545e+00, dfl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.21637021355784e-01, dfl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.88561808316413e+00, dfr[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.63299316185545e+00, dfr[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.21637021355784e-01, dfr[2], 1e-12) );
    } else if (test_func_op==1) {
      TEST_CHECK( gkyl_compare(-9.42809041582064e-01, dfl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.99959962173949e-17, dfl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.21637021355784e-01, dfl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 9.42809041582064e-01, dfr[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.99959962173949e-17, dfr[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-4.21637021355784e-01, dfr[2], 1e-12) );
    }
  }

  gkyl_eval_on_nodes_release(evup);
  gkyl_array_release(distf);
}

void evalFunc_2x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = x*x + x*y + y*y;
}

void evalFunc_2x_trig(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = sin((2.*M_PI/(4.))*x) + sin((2.*M_PI/(4.))*x)*cos((2.*M_PI/(4.))*y) + cos((2.*M_PI/(4.))*y);
}

void test_2x(int poly_order, int test_func_op)
{
  double lower[] = {-2.0, -2.0}, upper[] = {2.0, 2.0};
  int cells[] = {2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_eval_on_nodes *evup;
  if (test_func_op==0)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_2x_quad, NULL);
  else if (test_func_op==1)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_2x_trig, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_eval_on_nodes_advance(evup, 0.0, &arr_range, distf);
  gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, "ctest_eval_on_nodes_distf.gkyl");

  double *dfll = gkyl_array_fetch(distf, 0);  // left, low cell
  double *dflu = gkyl_array_fetch(distf, 1);  // left, up cell
  double *dfrl = gkyl_array_fetch(distf, 2);  // right, low cell
  double *dfru = gkyl_array_fetch(distf, 3);  // right, up cell
  if (poly_order == 1) {
    if (test_func_op==0) {
      TEST_CHECK( gkyl_compare( 1.00000000000000e+01, dfll[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.46410161513776e+00, dfll[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.46410161513776e+00, dfll[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfll[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.00000000000000e+00, dflu[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00, dflu[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00, dflu[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dflu[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.00000000000000e+00, dfrl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00, dfrl[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00, dfrl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfrl[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.00000000000000e+01, dfru[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.46410161513776e+00, dfru[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.46410161513776e+00, dfru[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfru[3], 1e-12) );
    } else if (test_func_op==1) {
      TEST_CHECK( gkyl_compare(-1.11022302462516e-16 , dfll[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17 , dfll[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00 , dfll[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.77555756156289e-17 , dfll[3], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.11022302462516e-16 , dflu[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17 , dflu[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00 , dflu[2], 1e-12) );
      TEST_CHECK( gkyl_compare(-2.77555756156289e-17 , dflu[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.11022302462516e-16 , dfrl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17 , dfrl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00 , dfrl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.77555756156289e-17 , dfrl[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00 , dfru[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17 , dfru[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00 , dfru[2], 1e-12) );
      TEST_CHECK( gkyl_compare(-2.77555756156289e-17 , dfru[3], 1e-12) );
    }
  } else if (poly_order == 2) {
    if (test_func_op==0) {
      TEST_CHECK( gkyl_compare( 7.33333333333333e+00, dfll[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.46410161513775e+00, dfll[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.46410161513775e+00, dfll[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfll[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999945e-01, dfll[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999945e-01, dfll[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.11022302462516e-16, dfll[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.11022302462516e-16, dfll[7], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.33333333333333e+00, dflu[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00, dflu[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00, dflu[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dflu[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dflu[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dflu[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17, dflu[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17, dflu[7], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.33333333333333e+00, dfrl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00, dfrl[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00, dfrl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfrl[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dfrl[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dfrl[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17, dfrl[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17, dfrl[7], 1e-12) );
      TEST_CHECK( gkyl_compare( 7.33333333333333e+00, dfru[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.46410161513775e+00, dfru[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.46410161513775e+00, dfru[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01, dfru[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dfru[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dfru[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dfru[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dfru[7], 1e-12) );
    } else if (test_func_op==1) {
      TEST_CHECK( gkyl_compare(-1.33333333333333e+00, dfll[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.93889390390723e-17, dfll[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.84900179459750e-01, dfll[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.77555756156289e-17, dfll[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dfll[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dfll[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.44265186329548e-01, dfll[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dfll[7], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.33333333333333e+00, dflu[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.93889390390723e-17, dflu[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.84900179459750e-01, dflu[2], 1e-12) );
      TEST_CHECK( gkyl_compare(-2.77555756156289e-17, dflu[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.96284793999944e-01, dflu[4], 1e-12) );
      TEST_CHECK( gkyl_compare(-2.77555756156289e-17, dflu[5], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.44265186329548e-01, dflu[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dflu[7], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.33333333333333e+00, dfrl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.55111512312578e-17, dfrl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.92450089729875e+00, dfrl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.77555756156289e-17, dfrl[3], 1e-12) );
      TEST_CHECK( gkyl_compare(-5.96284793999944e-01, dfrl[4], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.38777878078145e-16, dfrl[5], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.44265186329548e-01, dfrl[6], 1e-12) );
      TEST_CHECK( gkyl_compare( 0.00000000000000e+00, dfrl[7], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.33333333333333e+00, dfru[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.93889390390723e-17, dfru[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.92450089729875e+00, dfru[2], 1e-12) );
      TEST_CHECK( gkyl_compare(-2.77555756156289e-17, dfru[3], 1e-12) );
      TEST_CHECK( gkyl_compare(-5.96284793999944e-01, dfru[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 8.32667268468867e-17, dfru[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.44265186329548e-01, dfru[6], 1e-12) );
      TEST_CHECK( gkyl_compare(-6.93889390390723e-17, dfru[7], 1e-12) );
    }
  }

  gkyl_eval_on_nodes_release(evup);
  gkyl_array_release(distf);
}

void test_1x1v_hyb(int poly_order, int test_func_op)
{
  double lower[] = {-2.0, -2.0}, upper[] = {2.0, 2.0};
  int cells[] = {2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  int cdim = 1; int vdim = 1;
  gkyl_cart_modal_hybrid(&basis, cdim, vdim);

  // projection updater for dist-function
  gkyl_eval_on_nodes *evup;
  if (test_func_op==0)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_2x_quad, NULL);
  else if (test_func_op==1)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_2x_trig, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_eval_on_nodes_advance(evup, 0.0, &arr_range, distf);
  gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, "ctest_eval_on_nodes_distf_hyb.gkyl");

  double *dfll = gkyl_array_fetch(distf, 0);  // left, low cell
  double *dflu = gkyl_array_fetch(distf, 1);  // left, up cell
  double *dfrl = gkyl_array_fetch(distf, 2);  // right, low cell
  double *dfru = gkyl_array_fetch(distf, 3);  // right, up cell
  if (poly_order == 1) {
    if (test_func_op==0) {
      TEST_CHECK( gkyl_compare( 8.6666666666666679e+00, dfll[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.46410161513776e+00  , dfll[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.46410161513776e+00  , dfll[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01  , dfll[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.9628479399994472e-01, dfll[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.7755575615628914e-16, dfll[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.6666666666666670e+00, dflu[0], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00  , dflu[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00  , dflu[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01  , dflu[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.9628479399994450e-01, dflu[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.1102230246251565e-16, dflu[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 4.6666666666666670e+00, dfrl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.15470053837925e+00  , dfrl[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.15470053837925e+00  , dfrl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01  , dfrl[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.9628479399994450e-01, dfrl[4], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.1102230246251565e-16, dfrl[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 8.6666666666666679e+00, dfru[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.46410161513776e+00  , dfru[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.46410161513776e+00  , dfru[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.66666666666667e-01  , dfru[3], 1e-12) );
      TEST_CHECK( gkyl_compare( 5.9628479399994472e-01, dfru[4], 1e-12) );
      TEST_CHECK( gkyl_compare(-2.7755575615628914e-16, dfru[5], 1e-12) );
    } else if (test_func_op==1) {
      TEST_CHECK( gkyl_compare(-3.7007434154171901e-17 , dfll[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.8502929347591034e-17 , dfll[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.1547005383792517e+00 , dfll[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 3.7007434154171895e-17 , dfll[3], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.3100455376630308e-17 , dfll[4], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.9696061028207229e-18 , dfll[5], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.7007434154171901e-17 , dflu[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 6.8502929347591034e-17 , dflu[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.1547005383792517e+00 , dflu[2], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.7007434154171895e-17 , dflu[3], 1e-12) );
      TEST_CHECK( gkyl_compare(-3.3100455376630308e-17 , dflu[4], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.9696061028207229e-18 , dflu[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.3730110819465753e-16 , dfrl[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 8.9869181418519517e-17 , dfrl[1], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.1547005383792517e+00 , dfrl[2], 1e-12) );
      TEST_CHECK( gkyl_compare( 7.4014868308343790e-17 , dfrl[3], 1e-12) );
      TEST_CHECK( gkyl_compare(-6.8229156819663812e-18 , dfrl[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.7140950719175985e-17 , dfrl[5], 1e-12) );
      TEST_CHECK( gkyl_compare( 2.3730110819465753e-16 , dfru[0], 1e-12) );
      TEST_CHECK( gkyl_compare( 8.9869181418519517e-17 , dfru[1], 1e-12) );
      TEST_CHECK( gkyl_compare(-1.1547005383792517e+00 , dfru[2], 1e-12) );
      TEST_CHECK( gkyl_compare(-7.4014868308343790e-17 , dfru[3], 1e-12) );
      TEST_CHECK( gkyl_compare(-6.8229156819663812e-18 , dfru[4], 1e-12) );
      TEST_CHECK( gkyl_compare( 1.7140950719175985e-17 , dfru[5], 1e-12) );
    }
  } else {
    assert(true);
  }
  gkyl_eval_on_nodes_release(evup);
  gkyl_array_release(distf);
}

void evalFunc_3x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = x*x + x*y + y*y + x*z + y*z + z*z;
}

void evalFunc_3x_trig(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = sin((2.*M_PI/(4.))*x) + sin((2.*M_PI/(4.))*x)*cos((2.*M_PI/(4.))*y) + cos((2.*M_PI/(4.))*y)
           + sin((2.*M_PI/(4.))*x)*cos((2.*M_PI/(4.))*z)
           + cos((2.*M_PI/(4.))*y)*cos((2.*M_PI/(4.))*z) + cos((2.*M_PI/(4.))*z);
}

void test_3x(int poly_order, int test_func_op)
{
  double lower[] = {-2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0};
  int cells[] = {2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_eval_on_nodes *evup;
  if (test_func_op==0)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_3x_quad, NULL);
  else if (test_func_op==1)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_3x_trig, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_eval_on_nodes_advance(evup, 0.0, &arr_range, distf);

  gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, "ctest_eval_on_nodes_distf.gkyl");

  double *dfllf = gkyl_array_fetch(distf, 0);  // left, low, front cell
  double *dfllb = gkyl_array_fetch(distf, 1);  // left, low, back cell
  double *dfluf = gkyl_array_fetch(distf, 2);  // left, up, front cell
  double *dflub = gkyl_array_fetch(distf, 3);  // left, up, back cell
  double *dfrlf = gkyl_array_fetch(distf, 4);  // right, low, front cell
  double *dfrlb = gkyl_array_fetch(distf, 5);  // right, low, back cell
  double *dfruf = gkyl_array_fetch(distf, 6);  // right, up, front cell
  double *dfrub = gkyl_array_fetch(distf, 7);  // right, up, back cell
  if (poly_order == 1) {
    if (test_func_op==0) {
      double fref[] = {
         2.54558441227157e+01, -6.53197264742181e+00, -6.53197264742181e+00, -6.53197264742181e+00, 9.42809041582064e-01, 9.42809041582062e-01,
         9.42809041582062e-01, -5.55111512312578e-17, 1.41421356237310e+01, -3.26598632371091e+00, -3.26598632371090e+00, 7.77156117237610e-16,
         9.42809041582064e-01, 9.42809041582063e-01, 9.42809041582063e-01, -5.55111512312578e-17, 1.41421356237309e+01, -3.26598632371091e+00,
         3.33066907387547e-16, -3.26598632371090e+00, 9.42809041582064e-01, 9.42809041582063e-01, 9.42809041582064e-01, 5.55111512312578e-17,
         1.41421356237310e+01, 0.00000000000000e+00, 3.26598632371091e+00, 3.26598632371091e+00, 9.42809041582064e-01, 9.42809041582064e-01,
         9.42809041582064e-01, 0.00000000000000e+00, 1.41421356237309e+01, 2.22044604925031e-16, -3.26598632371090e+00, -3.26598632371090e+00,
         9.42809041582064e-01, 9.42809041582064e-01, 9.42809041582063e-01, 0.00000000000000e+00, 1.41421356237310e+01, 3.26598632371091e+00,
         8.88178419700125e-16, 3.26598632371091e+00, 9.42809041582064e-01, 9.42809041582064e-01, 9.42809041582064e-01, 0.00000000000000e+00,
         1.41421356237310e+01, 3.26598632371091e+00, 3.26598632371091e+00, 8.88178419700125e-16, 9.42809041582064e-01, 9.42809041582064e-01,
         9.42809041582064e-01, 0.00000000000000e+00, 2.54558441227157e+01, 6.53197264742181e+00, 6.53197264742181e+00, 6.53197264742181e+00,
         9.42809041582064e-01, 9.42809041582065e-01, 9.42809041582065e-01, 0.00000000000000e+00,
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    } else if (test_func_op==1) {
      double fref[] = {
         8.88178419700125e-16, 2.22044604925031e-16, 1.63299316185545e+00, 1.63299316185545e+00, 5.55111512312578e-17, 5.55111512312578e-17,
         9.42809041582064e-01, 2.77555756156289e-17, -1.22124532708767e-15, 8.32667268468867e-17, 1.63299316185545e+00, -1.63299316185545e+00,
         6.93889390390723e-17, -6.93889390390723e-17, -9.42809041582063e-01, 2.77555756156289e-17, 7.77156117237610e-16, -2.22044604925031e-16,
         -1.63299316185545e+00, 1.63299316185545e+00, -5.55111512312578e-17, 5.55111512312578e-17, -9.42809041582064e-01, -5.55111512312578e-17,
         -1.16573417585641e-15, 8.32667268468867e-17, -1.63299316185545e+00, -1.63299316185545e+00, -6.93889390390723e-17, -6.93889390390723e-17,
         9.42809041582063e-01, 1.38777878078145e-17, 6.66133814775094e-16, 2.22044604925031e-16, 1.63299316185545e+00, 1.63299316185545e+00,
         1.11022302462516e-16, 1.11022302462516e-16, 9.42809041582064e-01, 0.00000000000000e+00, -7.77156117237610e-16, 1.66533453693773e-16,
         1.63299316185545e+00, -1.63299316185545e+00, 8.32667268468867e-17, -8.32667268468867e-17, -9.42809041582063e-01, 2.77555756156289e-17,
         9.99200722162641e-16, -2.22044604925031e-16, -1.63299316185545e+00, 1.63299316185545e+00, -1.11022302462516e-16, 1.11022302462516e-16,
         -9.42809041582064e-01, -8.32667268468867e-17, -1.27675647831893e-15, 1.66533453693773e-16, -1.63299316185545e+00, -1.63299316185545e+00,
         -8.32667268468867e-17, -8.32667268468867e-17, 9.42809041582063e-01, 1.38777878078145e-17,
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    }
  } else if (poly_order == 2) {
    if (test_func_op==0) {
      double fref[] = {
         1.97989898732233e+01, -6.53197264742181e+00, -6.53197264742181e+00, -6.53197264742181e+00, 9.42809041582064e-01, 9.42809041582063e-01,
         9.42809041582064e-01, 8.43274042711568e-01, 8.43274042711567e-01, 8.43274042711568e-01, 0.00000000000000e+00, 1.52655665885959e-16,
         5.55111512312578e-17, 1.52655665885959e-16, 2.77555756156289e-16, 2.22044604925031e-16, 2.22044604925031e-16, 0.00000000000000e+00,
         0.00000000000000e+00, 0.00000000000000e+00, 8.48528137423857e+00, -3.26598632371090e+00, -3.26598632371090e+00, 2.22044604925031e-16,
         9.42809041582063e-01, 9.42809041582063e-01, 9.42809041582063e-01, 8.43274042711568e-01, 8.43274042711568e-01, 8.43274042711568e-01,
         0.00000000000000e+00, 5.55111512312578e-17, 5.55111512312578e-17, -5.55111512312578e-17, -5.55111512312578e-17, 1.66533453693773e-16,
         1.66533453693773e-16, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 8.48528137423858e+00, -3.26598632371090e+00,
         2.22044604925031e-16, -3.26598632371090e+00, 9.42809041582063e-01, 9.42809041582063e-01, 9.42809041582063e-01, 8.43274042711568e-01,
         8.43274042711568e-01, 8.43274042711568e-01, 0.00000000000000e+00, -5.55111512312578e-17, 1.66533453693773e-16, 5.55111512312578e-17,
         1.66533453693773e-16, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00,
         8.48528137423857e+00, 0.00000000000000e+00, 3.26598632371090e+00, 3.26598632371090e+00, 9.42809041582063e-01, 9.42809041582063e-01,
         9.42809041582063e-01, 8.43274042711568e-01, 8.43274042711568e-01, 8.43274042711568e-01, 0.00000000000000e+00, -1.11022302462516e-16,
         1.11022302462516e-16, -1.11022302462516e-16, -1.11022302462516e-16, 0.00000000000000e+00, 1.11022302462516e-16, 0.00000000000000e+00,
         0.00000000000000e+00, 0.00000000000000e+00, 8.48528137423857e+00, 0.00000000000000e+00, -3.26598632371090e+00, -3.26598632371090e+00,
         9.42809041582063e-01, 9.42809041582063e-01, 9.42809041582063e-01, 8.43274042711568e-01, 8.43274042711568e-01, 8.43274042711568e-01,
         0.00000000000000e+00, 1.66533453693773e-16, -5.55111512312578e-17, 1.66533453693773e-16, 1.66533453693773e-16, -5.55111512312578e-17,
         1.11022302462516e-16, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 8.48528137423857e+00, 3.26598632371090e+00,
         0.00000000000000e+00, 3.26598632371090e+00, 9.42809041582063e-01, 9.42809041582063e-01, 9.42809041582063e-01, 8.43274042711568e-01,
         8.43274042711568e-01, 8.43274042711568e-01, 0.00000000000000e+00, 1.11022302462516e-16, -2.22044604925031e-16, 0.00000000000000e+00,
         -2.22044604925031e-16, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00,
         8.48528137423857e+00, 3.26598632371090e+00, 3.26598632371090e+00, 4.44089209850063e-16, 9.42809041582063e-01, 9.42809041582063e-01,
         9.42809041582063e-01, 8.43274042711568e-01, 8.43274042711568e-01, 8.43274042711568e-01, 0.00000000000000e+00, 0.00000000000000e+00,
         -1.11022302462516e-16, 1.11022302462516e-16, 0.00000000000000e+00, -1.11022302462516e-16, -1.11022302462516e-16, 0.00000000000000e+00,
         0.00000000000000e+00, 0.00000000000000e+00, 1.97989898732233e+01, 6.53197264742181e+00, 6.53197264742181e+00, 6.53197264742181e+00,
         9.42809041582063e-01, 9.42809041582063e-01, 9.42809041582063e-01, 8.43274042711567e-01, 8.43274042711567e-01, 8.43274042711568e-01,
         0.00000000000000e+00, -2.22044604925031e-16, -2.22044604925031e-16, -2.22044604925031e-16, -2.22044604925031e-16, -2.22044604925031e-16,
         -2.22044604925031e-16, 0.00000000000000e+00, 0.00000000000000e+00, 0.00000000000000e+00,
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    } else if (test_func_op==1) {
      double fref[] = {
        -1.88561808316413e+00, 8.32667268468867e-17, 5.44331053951817e-01, 5.44331053951817e-01, 1.38777878078145e-17, 6.93889390390723e-17,
        9.42809041582063e-01, 8.43274042711568e-01, -5.55111512312578e-17, -5.55111512312578e-17, 0.00000000000000e+00, 4.86864495560148e-01,
        -2.77555756156289e-17, 4.86864495560148e-01, -8.32667268468867e-17, 2.77555756156289e-17, -8.32667268468867e-17, -1.38777878078145e-17,
        -1.38777878078145e-17, 2.77555756156289e-17, -1.88561808316413e+00, 5.55111512312578e-17, 5.44331053951817e-01, -5.44331053951817e-01,
        3.46944695195361e-17, -8.32667268468867e-17, -9.42809041582063e-01, 8.43274042711568e-01, -2.77555756156289e-17, -1.38777878078145e-16,
        -1.38777878078145e-17, 4.86864495560148e-01, -4.16333634234434e-17, -4.86864495560147e-01, 1.11022302462516e-16, 2.08166817117217e-17,
        -2.77555756156289e-17, 0.00000000000000e+00, 0.00000000000000e+00, 2.08166817117217e-17, -1.88561808316413e+00, 5.55111512312578e-17,
        -5.44331053951817e-01, 5.44331053951817e-01, -1.38777878078145e-17, 8.32667268468867e-17, -9.42809041582063e-01, 8.43274042711568e-01,
        -5.55111512312578e-17, -8.32667268468867e-17, 0.00000000000000e+00, -4.86864495560148e-01, -1.38777878078145e-17, 4.86864495560148e-01,
        -9.71445146547012e-17, 2.77555756156289e-17, 8.32667268468867e-17, 1.38777878078145e-17, 0.00000000000000e+00, -1.38777878078145e-17,
        -1.88561808316413e+00, 6.93889390390723e-17, -5.44331053951817e-01, -5.44331053951817e-01, -6.93889390390723e-18, -9.02056207507940e-17,
        9.42809041582063e-01, 8.43274042711568e-01, -4.16333634234434e-17, -6.93889390390723e-17, 1.38777878078145e-17, -4.86864495560147e-01,
        -2.08166817117217e-17, -4.86864495560147e-01, 1.31838984174237e-16, 2.08166817117217e-17, 1.17961196366423e-16, 6.93889390390723e-18,
        6.93889390390723e-18, -2.08166817117217e-17, 1.88561808316413e+00, 8.32667268468867e-17, 2.72165526975909e+00, 2.72165526975909e+00,
        4.16333634234434e-17, 6.93889390390723e-17, 9.42809041582063e-01, -8.43274042711568e-01, -1.66533453693773e-16, -5.55111512312578e-17,
        2.77555756156289e-17, -4.86864495560148e-01, 0.00000000000000e+00, -4.86864495560148e-01, -8.32667268468867e-17, 0.00000000000000e+00,
        -2.77555756156289e-17, 1.38777878078145e-17, 0.00000000000000e+00, 1.38777878078145e-17, 1.88561808316413e+00, 1.11022302462516e-16,
        2.72165526975909e+00, -2.72165526975909e+00, 6.24500451351651e-17, -8.32667268468867e-17, -9.42809041582063e-01, -8.43274042711568e-01,
        -1.94289029309402e-16, -1.94289029309402e-16, -1.38777878078145e-17, -4.86864495560148e-01, -1.38777878078145e-17, 4.86864495560148e-01,
        9.71445146547012e-17, -1.38777878078145e-17, -4.16333634234434e-17, 0.00000000000000e+00, 0.00000000000000e+00, 6.93889390390723e-18,
        1.88561808316413e+00, 1.11022302462516e-16, -2.72165526975909e+00, 2.72165526975909e+00, -4.16333634234434e-17, 8.32667268468867e-17,
        -9.42809041582063e-01, -8.43274042711568e-01, -1.38777878078145e-16, 2.77555756156289e-17, -2.77555756156289e-17, 4.86864495560148e-01,
        -2.77555756156289e-17, -4.86864495560148e-01, -1.24900090270330e-16, 0.00000000000000e+00, 5.55111512312578e-17, -1.38777878078145e-17,
        0.00000000000000e+00, 0.00000000000000e+00, 1.88561808316413e+00, 6.93889390390723e-17, -2.72165526975909e+00, -2.72165526975909e+00,
        -3.46944695195361e-17, -9.02056207507940e-17, 9.42809041582063e-01, -8.43274042711568e-01, -1.66533453693773e-16, -1.66533453693773e-16,
        1.38777878078145e-17, 4.86864495560148e-01, -1.38777878078145e-17, 4.86864495560148e-01, 8.32667268468867e-17, -1.38777878078145e-17,
        4.16333634234434e-17, -6.93889390390723e-18, -6.93889390390723e-18, -6.93889390390723e-18,
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    }
  }

  gkyl_eval_on_nodes_release(evup);
  gkyl_array_release(distf);
}

void evalFunc_4x_quad(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2], vx = xn[3];
  fout[0] = x*x + x*y + x*z + x*vx + y*y + y*z + y*vx + z*z + z*vx + vx*vx;
}

void evalFunc_4x_trig(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2], vx = xn[3];
  fout[0] =  sin((2.*M_PI/(4.))*x) 
           + sin((2.*M_PI/(4.))*x)*cos((2.*M_PI/(4.))*y) 
           + sin((2.*M_PI/(4.))*x)*cos((2.*M_PI/(4.))*z)
           + sin((2.*M_PI/(4.))*x)*cos((2.*M_PI/(4.))*vx)
           + cos((2.*M_PI/(4.))*y)
           + cos((2.*M_PI/(4.))*y)*cos((2.*M_PI/(4.))*z) 
           + cos((2.*M_PI/(4.))*y)*cos((2.*M_PI/(4.))*vx) 
           + cos((2.*M_PI/(4.))*z)
           + cos((2.*M_PI/(4.))*z)*cos((2.*M_PI/(4.))*vx) 
           + cos((2.*M_PI/(4.))*vx);
}


void test_4x(int poly_order, int test_func_op)
{
  double lower[] = {-2.0, -2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0, 2.0};
  int cells[] = {2, 2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // projection updater for dist-function
  gkyl_eval_on_nodes *evup;
  if (test_func_op==0)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_4x_quad, NULL);
  else if (test_func_op==1)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_4x_trig, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_eval_on_nodes_advance(evup, 0.0, &arr_range, distf);
  gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, "ctest_eval_on_nodes_distf.gkyl");

  if (poly_order == 1) {
    if (test_func_op==0) {
      double fref[] = {
         5.6000000000000000e+01, -1.1547005383792520e+01, -1.1547005383792520e+01, -1.1547005383792520e+01, 
        -1.1547005383792520e+01, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        3.2000000000000000e+01, -6.9282032302755114e+00, -6.9282032302755114e+00, -6.9282032302755114e+00, 
        -2.3094010767585038e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        3.2000000000000000e+01, -6.9282032302755114e+00, -6.9282032302755114e+00, -2.3094010767585038e+00, 
        -6.9282032302755114e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.4000000000000000e+01, -2.3094010767585038e+00, -2.3094010767585038e+00, 2.3094010767585038e+00, 
        2.3094010767585038e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        3.2000000000000000e+01, -6.9282032302755114e+00, -2.3094010767585038e+00, -6.9282032302755114e+00, 
        -6.9282032302755114e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.4000000000000000e+01, -2.3094010767585038e+00, 2.3094010767585038e+00, -2.3094010767585038e+00, 
        2.3094010767585038e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.4000000000000000e+01, -2.3094010767585038e+00, 2.3094010767585038e+00, 2.3094010767585038e+00, 
        -2.3094010767585038e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        3.2000000000000000e+01, 2.3094010767585038e+00, 6.9282032302755114e+00, 6.9282032302755114e+00, 
        6.9282032302755114e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        3.2000000000000000e+01, -2.3094010767585038e+00, -6.9282032302755114e+00, -6.9282032302755114e+00, 
        -6.9282032302755114e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.4000000000000000e+01, 2.3094010767585038e+00, -2.3094010767585038e+00, -2.3094010767585038e+00, 
        2.3094010767585038e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.4000000000000000e+01, 2.3094010767585038e+00, -2.3094010767585038e+00, 2.3094010767585038e+00, 
        -2.3094010767585038e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        3.2000000000000000e+01, 6.9282032302755114e+00, 2.3094010767585038e+00, 6.9282032302755114e+00, 
        6.9282032302755114e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.4000000000000000e+01, 2.3094010767585038e+00, 2.3094010767585038e+00, -2.3094010767585038e+00, 
        -2.3094010767585038e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        3.2000000000000000e+01, 6.9282032302755114e+00, 6.9282032302755114e+00, 2.3094010767585038e+00, 
        6.9282032302755114e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        3.2000000000000000e+01, 6.9282032302755114e+00, 6.9282032302755114e+00, 6.9282032302755114e+00, 
        2.3094010767585038e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        5.6000000000000000e+01, 1.1547005383792520e+01, 1.1547005383792520e+01, 1.1547005383792520e+01, 
        1.1547005383792520e+01, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00,
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    } else if (test_func_op==1) {
      double fref[] = {
          0.0000000000000000e+00, -3.2049378106392743e-17, 2.3094010767585038e+00, 2.3094010767585038e+00, 
          2.3094010767585038e+00, 1.8503717077085941e-17, 1.8503717077085941e-17, 1.3333333333333333e+00, 
          1.8503717077085941e-17, 1.3333333333333333e+00, 1.3333333333333333e+00, -1.0683126035464246e-17, 
          -3.2049378106392736e-17, -3.2049378106392736e-17, -1.0683126035464246e-17, 6.1679056923619812e-18, 
          -2.2204460492503131e-16, -3.2049378106392743e-17, 2.3094010767585038e+00, 2.3094010767585038e+00, 
          -2.3094010767585038e+00, -1.8503717077085941e-17, 1.8503717077085941e-17, 1.3333333333333333e+00, 
          1.8503717077085941e-17, -1.3333333333333333e+00, -1.3333333333333333e+00, -1.0683126035464246e-17, 
          1.0683126035464246e-17, -2.1366252070928492e-17, -4.2732504141856984e-17, 3.0839528461809905e-17, 
          0.0000000000000000e+00, -3.2049378106392743e-17, 2.3094010767585038e+00, -2.3094010767585038e+00, 
          2.3094010767585038e+00, 1.8503717077085941e-17, -5.5511151231257827e-17, -1.3333333333333333e+00, 
          1.8503717077085941e-17, 1.3333333333333333e+00, -1.3333333333333333e+00, -1.0683126035464246e-17, 
          -5.3415630177321232e-17, 1.0683126035464246e-17, 0.0000000000000000e+00, 1.8503717077085944e-17, 
          0.0000000000000000e+00, 6.4098756212785485e-17, 2.3094010767585038e+00, -2.3094010767585038e+00, 
          -2.3094010767585038e+00, 3.7007434154171883e-17, -7.4014868308343765e-17, -1.3333333333333333e+00, 
          -3.7007434154171883e-17, -1.3333333333333333e+00, 1.3333333333333333e+00, 2.1366252070928492e-17, 
          6.4098756212785473e-17, 0.0000000000000000e+00, 2.1366252070928492e-17, -1.2335811384723962e-17, 
          0.0000000000000000e+00, 0.0000000000000000e+00, -2.3094010767585038e+00, 2.3094010767585038e+00, 
          2.3094010767585038e+00, -7.4014868308343765e-17, 0.0000000000000000e+00, -1.3333333333333333e+00, 
          0.0000000000000000e+00, -1.3333333333333333e+00, 1.3333333333333333e+00, 1.0683126035464246e-17, 
          4.2732504141856984e-17, -3.2049378106392736e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
          -1.6653345369377348e-16, 3.2049378106392743e-17, -2.3094010767585038e+00, 2.3094010767585038e+00, 
          -2.3094010767585038e+00, -1.1102230246251565e-16, -1.8503717077085941e-17, -1.3333333333333333e+00, 
          -1.8503717077085941e-17, 1.3333333333333333e+00, -1.3333333333333333e+00, 0.0000000000000000e+00, 
          -2.1366252070928492e-17, 0.0000000000000000e+00, 1.0683126035464246e-17, -1.8503717077085944e-17, 
          0.0000000000000000e+00, 0.0000000000000000e+00, -2.3094010767585038e+00, -2.3094010767585038e+00, 
          2.3094010767585038e+00, -5.5511151231257827e-17, -3.7007434154171883e-17, 1.3333333333333333e+00, 
          0.0000000000000000e+00, -1.3333333333333333e+00, -1.3333333333333333e+00, 3.2049378106392736e-17, 
          3.2049378106392736e-17, 0.0000000000000000e+00, 2.1366252070928492e-17, -1.2335811384723962e-17, 
          0.0000000000000000e+00, 0.0000000000000000e+00, -2.3094010767585038e+00, -2.3094010767585038e+00, 
          -2.3094010767585038e+00, -1.4802973661668753e-16, -7.4014868308343765e-17, 1.3333333333333333e+00, 
          0.0000000000000000e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
          0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 4.9343245538895850e-17, 
          2.2204460492503131e-16, 3.8459253727671291e-16, 2.3094010767585038e+00, 2.3094010767585038e+00, 
          2.3094010767585038e+00, 1.8503717077085941e-16, 1.4802973661668753e-16, 1.3333333333333333e+00, 
          7.4014868308343765e-17, 1.3333333333333333e+00, 1.3333333333333333e+00, 8.5465008283713968e-17, 
          2.1366252070928492e-17, 0.0000000000000000e+00, -4.2732504141856984e-17, 0.0000000000000000e+00, 
          3.3306690738754696e-16, 3.8459253727671291e-16, 2.3094010767585038e+00, 2.3094010767585038e+00, 
          -2.3094010767585038e+00, 1.8503717077085941e-16, 1.4802973661668753e-16, 1.3333333333333333e+00, 
          -7.4014868308343765e-17, -1.3333333333333333e+00, -1.3333333333333333e+00, 4.2732504141856984e-17, 
          -2.1366252070928492e-17, 2.1366252070928492e-17, 4.2732504141856984e-17, 0.0000000000000000e+00, 
          2.2204460492503131e-16, 3.8459253727671291e-16, 2.3094010767585038e+00, -2.3094010767585038e+00, 
          2.3094010767585038e+00, 2.2204460492503131e-16, -1.4802973661668753e-16, -1.3333333333333333e+00, 
          7.4014868308343765e-17, 1.3333333333333333e+00, -1.3333333333333333e+00, -4.2732504141856984e-17, 
          4.2732504141856984e-17, 0.0000000000000000e+00, 4.2732504141856984e-17, 2.4671622769447925e-17, 
          4.4408920985006262e-16, 4.1664191538310567e-16, 2.3094010767585038e+00, -2.3094010767585038e+00, 
          -2.3094010767585038e+00, 1.4802973661668753e-16, -1.1102230246251565e-16, -1.3333333333333333e+00, 
          -9.2518585385429703e-17, -1.3333333333333333e+00, 1.3333333333333330e+00, -4.2732504141856984e-17, 
          0.0000000000000000e+00, -4.2732504141856984e-17, -1.0683126035464246e-17, -1.2335811384723962e-17, 
          2.2204460492503131e-16, 3.8459253727671291e-16, -2.3094010767585038e+00, 2.3094010767585038e+00, 
          2.3094010767585038e+00, -1.8503717077085941e-16, 1.4802973661668753e-16, -1.3333333333333333e+00, 
          7.4014868308343765e-17, -1.3333333333333333e+00, 1.3333333333333333e+00, -6.4098756212785473e-17, 
          -2.1366252070928492e-17, 0.0000000000000000e+00, 4.2732504141856984e-17, 1.2335811384723962e-17, 
          5.5511151231257827e-16, 3.2049378106392745e-16, -2.3094010767585038e+00, 2.3094010767585038e+00, 
          -2.3094010767585038e+00, -1.8503717077085941e-16, 1.4802973661668753e-16, -1.3333333333333333e+00, 
          -3.7007434154171883e-17, 1.3333333333333333e+00, -1.3333333333333333e+00, -6.4098756212785473e-17, 
          2.1366252070928492e-17, 0.0000000000000000e+00, -2.1366252070928492e-17, 0.0000000000000000e+00, 
          6.6613381477509392e-16, 3.8459253727671291e-16, -2.3094010767585038e+00, -2.3094010767585038e+00, 
          2.3094010767585038e+00, -1.4802973661668753e-16, -7.4014868308343765e-17, 1.3333333333333333e+00, 
          7.4014868308343765e-17, -1.3333333333333333e+00, -1.3333333333333333e+00, 4.2732504141856984e-17, 
          0.0000000000000000e+00, 0.0000000000000000e+00, -8.5465008283713968e-17, -2.4671622769447925e-17, 
          4.4408920985006262e-16, 3.8459253727671291e-16, -2.3094010767585038e+00, -2.3094010767585038e+00, 
          -2.3094010767585038e+00, -1.4802973661668753e-16, -1.4802973661668753e-16, 1.3333333333333333e+00, 
          -7.4014868308343765e-17, 1.3333333333333333e+00, 1.3333333333333330e+00, 4.2732504141856984e-17, 
          0.0000000000000000e+00, -4.2732504141856984e-17, 0.0000000000000000e+00, 2.4671622769447925e-17, 
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    }
  } else if (poly_order == 2) {
    if (test_func_op==0) {
      double fref[] = {
        4.5333333333333321e+01, -1.1547005383792515e+01, -1.1547005383792516e+01, -1.1547005383792516e+01, 
        -1.1547005383792516e+01, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998819e+00, 
        1.1925695879998819e+00, 1.1925695879998819e+00, 1.1925695879998823e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 5.5511151231257827e-16, 
        5.5511151231257827e-16, 5.5511151231257827e-16, 5.5511151231257827e-16, 4.4408920985006262e-16, 
        5.5511151231257827e-16, 5.5511151231257827e-16, 5.5511151231257827e-16, 5.5511151231257827e-16, 
        5.5511151231257827e-16, 5.5511151231257827e-16, 5.5511151231257827e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.1333333333333321e+01, -6.9282032302755097e+00, -6.9282032302755097e+00, -6.9282032302755097e+00, 
        -2.3094010767585029e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998846e+00, 
        1.1925695879998846e+00, 1.1925695879998846e+00, 1.1925695879998841e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.3306690738754696e-16, 
        3.3306690738754696e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 4.4408920985006262e-16, 
        3.3306690738754696e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        3.3306690738754696e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.1333333333333321e+01, -6.9282032302755097e+00, -6.9282032302755097e+00, -2.3094010767585029e+00, 
        -6.9282032302755097e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998846e+00, 
        1.1925695879998846e+00, 1.1925695879998846e+00, 1.1925695879998841e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.3306690738754696e-16, 
        3.3306690738754696e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 4.4408920985006262e-16, 
        3.3306690738754696e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 
        3.3306690738754696e-16, 3.3306690738754696e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.3333333333333327e+01, -2.3094010767585029e+00, -2.3094010767585029e+00, 2.3094010767585029e+00, 
        2.3094010767585029e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998854e+00, 
        1.1925695879998854e+00, 1.1925695879998854e+00, 1.1925695879998859e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1102230246251565e-16, 
        1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, 1.1102230246251565e-16, 
        1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        1.1102230246251565e-16, 1.1102230246251565e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.1333333333333321e+01, -6.9282032302755097e+00, -2.3094010767585029e+00, -6.9282032302755097e+00, 
        -6.9282032302755097e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998846e+00, 
        1.1925695879998846e+00, 1.1925695879998846e+00, 1.1925695879998841e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1102230246251565e-16, 
        3.3306690738754696e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 4.4408920985006262e-16, 
        1.1102230246251565e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 
        3.3306690738754696e-16, 1.1102230246251565e-16, 3.3306690738754696e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.3333333333333327e+01, -2.3094010767585029e+00, 2.3094010767585029e+00, -2.3094010767585029e+00, 
        2.3094010767585029e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998854e+00, 
        1.1925695879998854e+00, 1.1925695879998854e+00, 1.1925695879998859e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.1102230246251565e-16, 
        1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        -1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        1.1102230246251565e-16, -1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.3333333333333327e+01, -2.3094010767585029e+00, 2.3094010767585029e+00, 2.3094010767585029e+00, 
        -2.3094010767585029e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998854e+00, 
        1.1925695879998854e+00, 1.1925695879998854e+00, 1.1925695879998859e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.1102230246251565e-16, 
        1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, 1.1102230246251565e-16, 
        -1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.1333333333333321e+01, 2.3094010767585029e+00, 6.9282032302755097e+00, 6.9282032302755097e+00, 
        6.9282032302755097e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998846e+00, 
        1.1925695879998846e+00, 1.1925695879998846e+00, 1.1925695879998841e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -3.3306690738754696e-16, 
        -1.1102230246251565e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, -1.1102230246251565e-16, 
        -3.3306690738754696e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, 
        -1.1102230246251565e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.1333333333333321e+01, -2.3094010767585029e+00, -6.9282032302755097e+00, -6.9282032302755097e+00, 
        -6.9282032302755097e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998846e+00, 
        1.1925695879998846e+00, 1.1925695879998846e+00, 1.1925695879998841e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.3306690738754696e-16, 
        1.1102230246251565e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 1.1102230246251565e-16, 
        3.3306690738754696e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 
        1.1102230246251565e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.3333333333333327e+01, 2.3094010767585029e+00, -2.3094010767585029e+00, -2.3094010767585029e+00, 
        2.3094010767585029e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998854e+00, 
        1.1925695879998854e+00, 1.1925695879998854e+00, 1.1925695879998859e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1102230246251565e-16, 
        -1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, -1.1102230246251565e-16, 
        1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        -1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.3333333333333327e+01, 2.3094010767585029e+00, -2.3094010767585029e+00, 2.3094010767585029e+00, 
        -2.3094010767585029e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998854e+00, 
        1.1925695879998854e+00, 1.1925695879998854e+00, 1.1925695879998859e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1102230246251565e-16, 
        -1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        -1.1102230246251565e-16, 1.1102230246251565e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.1333333333333321e+01, 6.9282032302755097e+00, 2.3094010767585029e+00, 6.9282032302755097e+00, 
        6.9282032302755097e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998846e+00, 
        1.1925695879998846e+00, 1.1925695879998846e+00, 1.1925695879998841e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.1102230246251565e-16, 
        -3.3306690738754696e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, -4.4408920985006262e-16, 
        -1.1102230246251565e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, 
        -3.3306690738754696e-16, -1.1102230246251565e-16, -3.3306690738754696e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.3333333333333327e+01, 2.3094010767585029e+00, 2.3094010767585029e+00, -2.3094010767585029e+00, 
        -2.3094010767585029e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998854e+00, 
        1.1925695879998854e+00, 1.1925695879998854e+00, 1.1925695879998859e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.1102230246251565e-16, 
        -1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, -1.1102230246251565e-16, 
        -1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        -1.1102230246251565e-16, -1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.1333333333333321e+01, 6.9282032302755097e+00, 6.9282032302755097e+00, 2.3094010767585029e+00, 
        6.9282032302755097e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998846e+00, 
        1.1925695879998846e+00, 1.1925695879998846e+00, 1.1925695879998841e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -3.3306690738754696e-16, 
        -3.3306690738754696e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, -4.4408920985006262e-16, 
        -3.3306690738754696e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, 
        -3.3306690738754696e-16, -3.3306690738754696e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.1333333333333321e+01, 6.9282032302755097e+00, 6.9282032302755097e+00, 6.9282032302755097e+00, 
        2.3094010767585029e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998846e+00, 
        1.1925695879998846e+00, 1.1925695879998846e+00, 1.1925695879998841e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -3.3306690738754696e-16, 
        -3.3306690738754696e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, -4.4408920985006262e-16, 
        -3.3306690738754696e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        -3.3306690738754696e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        4.5333333333333321e+01, 1.1547005383792515e+01, 1.1547005383792516e+01, 1.1547005383792516e+01, 
        1.1547005383792516e+01, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998819e+00, 
        1.1925695879998819e+00, 1.1925695879998819e+00, 1.1925695879998823e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -5.5511151231257827e-16, 
        -5.5511151231257827e-16, -5.5511151231257827e-16, -5.5511151231257827e-16, -4.4408920985006262e-16, 
        -5.5511151231257827e-16, -5.5511151231257827e-16, -5.5511151231257827e-16, -5.5511151231257827e-16, 
        -5.5511151231257827e-16, -5.5511151231257827e-16, -5.5511151231257827e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    } else if (test_func_op==1) {
      double fref[] = {
        -2.6666666666666665e+00, 2.0297939467382067e-16, 7.6980035891950083e-01, 7.6980035891950083e-01, 
        7.6980035891950083e-01, 4.3175339846533860e-17, 6.7846962615981782e-17, 1.3333333333333333e+00, 
        2.4671622769447922e-17, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998881e+00, 
        0.0000000000000000e+00, -3.3100455376630314e-17, -1.6550227688315157e-17, -3.5610420118214158e-18, 
        -4.6293546153678406e-17, -1.0683126035464248e-17, -3.5610420118214158e-18, 6.8853037265909611e-01, 
        3.8221113643993394e-17, 6.8853037265909611e-01, -1.1102230246251565e-16, -2.8665835232995051e-17, 
        -1.1102230246251565e-16, 6.8853037265909611e-01, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        -5.7331670465990103e-17, -1.1102230246251565e-16, -1.1102230246251565e-16, 6.1679056923619812e-18, 
        0.0000000000000000e+00, 2.7583712813858585e-17, -5.5167425627717169e-18, 0.0000000000000000e+00, 
        2.7583712813858585e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, -5.5167425627717169e-18, 
        0.0000000000000000e+00, -4.4133940502173736e-17, -3.8617197939402019e-17, 0.0000000000000000e+00, 
        -3.1850928036661163e-18, -3.1850928036661163e-18, 9.5552784109983484e-18, -1.2740371214664465e-17, 
        -2.6666666666666661e+00, 3.3117690709939164e-16, 7.6980035891950083e-01, 7.6980035891950083e-01, 
        -7.6980035891950083e-01, 3.0839528461809905e-17, 7.4014868308343765e-17, 1.3333333333333333e+00, 
        -2.4671622769447922e-17, -1.3333333333333333e+00, -1.3333333333333333e+00, 1.1925695879998881e+00, 
        -6.6200910753260603e-17, -9.9301366129890917e-17, -8.2751138441575754e-17, -3.5610420118214158e-18, 
        3.9171462130035574e-17, -7.1220840236428317e-18, -1.4244168047285663e-17, 6.8853037265909611e-01, 
        1.9110556821996697e-17, 6.8853037265909611e-01, -1.1102230246251565e-16, -1.9110556821996700e-17, 
        -1.1102230246251565e-16, -6.8853037265909611e-01, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        -5.7331670465990103e-17, -1.1102230246251565e-16, -1.1102230246251565e-16, 3.0839528461809905e-17, 
        0.0000000000000000e+00, 5.5167425627717169e-18, -5.5167425627717169e-18, 0.0000000000000000e+00, 
        5.5167425627717169e-18, 0.0000000000000000e+00, 0.0000000000000000e+00, -5.5167425627717169e-18, 
        0.0000000000000000e+00, -4.4133940502173736e-17, -3.8617197939402019e-17, 0.0000000000000000e+00, 
        -1.2740371214664465e-17, 0.0000000000000000e+00, -1.5925464018330582e-17, -1.2740371214664465e-17, 
        -2.6666666666666665e+00, 2.0297939467382067e-16, 7.6980035891950083e-01, -7.6980035891950083e-01, 
        7.6980035891950083e-01, 6.7846962615981782e-17, -4.3175339846533860e-17, -1.3333333333333333e+00, 
        -4.9343245538895844e-17, 1.3333333333333333e+00, -1.3333333333333333e+00, 1.1925695879998881e+00, 
        0.0000000000000000e+00, -3.3100455376630314e-17, 3.3100455376630314e-17, -3.5610420118214158e-18, 
        -5.3415630177321238e-17, 3.5610420118214158e-18, 0.0000000000000000e+00, 6.8853037265909611e-01, 
        1.9110556821996697e-17, -6.8853037265909611e-01, 1.1102230246251565e-16, -1.9110556821996700e-17, 
        -1.1102230246251565e-16, 6.8853037265909611e-01, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        -5.7331670465990103e-17, -1.1102230246251565e-16, 3.4580075174528854e-17, 1.8503717077085944e-17, 
        0.0000000000000000e+00, -2.2066970251086868e-17, -5.5167425627717169e-18, 0.0000000000000000e+00, 
        5.5167425627717169e-18, 0.0000000000000000e+00, 0.0000000000000000e+00, -5.5167425627717169e-18, 
        0.0000000000000000e+00, -2.2066970251086868e-17, 1.6550227688315151e-17, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 3.1850928036661159e-18, 0.0000000000000000e+00, 
        -2.6666666666666665e+00, 1.0683126035464245e-16, 7.6980035891950083e-01, -7.6980035891950083e-01, 
        -7.6980035891950083e-01, 1.1102230246251565e-16, -4.9343245538895844e-17, -1.3333333333333333e+00, 
        0.0000000000000000e+00, -1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998881e+00, 
        0.0000000000000000e+00, -3.3100455376630314e-17, 3.3100455376630314e-17, 7.1220840236428317e-18, 
        4.9854588165499822e-17, 0.0000000000000000e+00, 7.1220840236428317e-18, 6.8853037265909611e-01, 
        1.9110556821996697e-17, -6.8853037265909611e-01, 1.1102230246251565e-16, -1.9110556821996700e-17, 
        -1.1102230246251565e-16, -6.8853037265909611e-01, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        -5.7331670465990103e-17, -1.1102230246251565e-16, 3.4580075174528854e-17, -1.2335811384723962e-17, 
        0.0000000000000000e+00, -1.1033485125543434e-17, 2.2066970251086868e-17, 0.0000000000000000e+00, 
        -1.1033485125543434e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, -3.3100455376630302e-17, 1.1033485125543434e-17, 0.0000000000000000e+00, 
        6.3701856073322325e-18, 0.0000000000000000e+00, -7.7037197775489434e-34, -6.3701856073322325e-18, 
        -2.6666666666666665e+00, 1.7093001656742794e-16, -7.6980035891950083e-01, 7.6980035891950083e-01, 
        7.6980035891950083e-01, -1.7270135938613544e-16, 5.5511151231257821e-17, -1.3333333333333333e+00, 
        0.0000000000000000e+00, -1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998881e+00, 
        0.0000000000000000e+00, -6.6200910753260628e-17, 0.0000000000000000e+00, 3.5610420118214158e-18, 
        4.6293546153678406e-17, -1.0683126035464248e-17, 0.0000000000000000e+00, -6.8853037265909611e-01, 
        0.0000000000000000e+00, 6.8853037265909611e-01, -1.1102230246251565e-16, -2.8665835232995051e-17, 
        1.1102230246251565e-16, 6.8853037265909611e-01, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        -5.7331670465990103e-17, 1.1102230246251565e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 1.1033485125543434e-17, -2.2066970251086868e-17, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.1033485125543434e-17, 
        0.0000000000000000e+00, 0.0000000000000000e+00, -4.4133940502173736e-17, -5.5511151231257827e-17, 
        0.0000000000000000e+00, -9.5552784109983484e-18, -6.3701856073322325e-18, 3.1850928036661163e-18, 
        -2.6666666666666661e+00, 2.6707815088660613e-16, -7.6980035891950083e-01, 7.6980035891950083e-01, 
        -7.6980035891950083e-01, -6.1679056923619798e-17, 9.8686491077791687e-17, -1.3333333333333333e+00, 
        -2.4671622769447922e-17, 1.3333333333333333e+00, -1.3333333333333333e+00, 1.1925695879998881e+00, 
        -4.9650683064945452e-17, -1.1585159381820608e-16, -4.9650683064945452e-17, 0.0000000000000000e+00, 
        -6.4098756212785485e-17, 0.0000000000000000e+00, 3.5610420118214158e-18, -6.8853037265909611e-01, 
        0.0000000000000000e+00, 6.8853037265909611e-01, -1.1102230246251565e-16, 0.0000000000000000e+00, 
        1.1102230246251565e-16, -6.8853037265909611e-01, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        -5.7331670465990103e-17, 1.1102230246251565e-16, -1.1102230246251565e-16, -1.8503717077085944e-17, 
        0.0000000000000000e+00, 1.6550227688315151e-17, -2.7583712813858585e-17, 0.0000000000000000e+00, 
        -1.1033485125543434e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, -5.5167425627717169e-18, 
        0.0000000000000000e+00, 0.0000000000000000e+00, -4.9650683064945452e-17, -5.5511151231257827e-17, 
        3.1850928036661163e-18, 0.0000000000000000e+00, 1.9110556821996697e-17, -9.5552784109983484e-18, 
        -2.6666666666666665e+00, 1.7093001656742794e-16, -7.6980035891950083e-01, -7.6980035891950083e-01, 
        7.6980035891950083e-01, -1.0485439677015368e-16, -4.9343245538895844e-17, 1.3333333333333333e+00, 
        -4.9343245538895844e-17, -1.3333333333333333e+00, -1.3333333333333333e+00, 1.1925695879998881e+00, 
        0.0000000000000000e+00, -6.6200910753260628e-17, 6.6200910753260628e-17, 1.0683126035464248e-17, 
        5.3415630177321238e-17, 0.0000000000000000e+00, 7.1220840236428317e-18, -6.8853037265909611e-01, 
        0.0000000000000000e+00, -6.8853037265909611e-01, 1.1102230246251565e-16, -1.9110556821996700e-17, 
        1.1102230246251565e-16, 6.8853037265909611e-01, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        -5.7331670465990103e-17, 1.1102230246251565e-16, 3.4580075174528854e-17, -1.2335811384723962e-17, 
        0.0000000000000000e+00, -3.3100455376630302e-17, -1.1033485125543434e-17, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -2.2066970251086868e-17, 
        0.0000000000000000e+00, 3.8617197939402019e-17, 2.2066970251086868e-17, 0.0000000000000000e+00, 
        6.3701856073322325e-18, -1.2740371214664465e-17, -9.5552784109983484e-18, 0.0000000000000000e+00, 
        -2.6666666666666665e+00, 1.7093001656742794e-16, -7.6980035891950083e-01, -7.6980035891950083e-01, 
        -7.6980035891950083e-01, -4.9343245538895837e-17, -4.9343245538895844e-17, 1.3333333333333333e+00, 
        0.0000000000000000e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.1925695879998881e+00, 
        0.0000000000000000e+00, -6.6200910753260628e-17, 6.6200910753260628e-17, 0.0000000000000000e+00, 
        -5.6976672189142654e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, -6.8853037265909611e-01, 
        3.8221113643993394e-17, -6.8853037265909611e-01, 1.1102230246251565e-16, 0.0000000000000000e+00, 
        1.1102230246251565e-16, -6.8853037265909611e-01, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        -5.7331670465990103e-17, 1.1102230246251565e-16, 3.4580075174528854e-17, 4.9343245538895850e-17, 
        0.0000000000000000e+00, -2.2066970251086868e-17, -2.2066970251086868e-17, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.1033485125543434e-17, 
        0.0000000000000000e+00, 2.2066970251086868e-17, 1.1033485125543434e-17, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 2.5480742429328930e-17, 0.0000000000000000e+00, 
        2.6666666666666661e+00, 4.2732504141856929e-17, 3.8490017945975055e+00, 3.8490017945975055e+00, 
        3.8490017945975055e+00, 1.7270135938613546e-16, 6.1679056923619798e-17, 1.3333333333333333e+00, 
        1.7270135938613544e-16, 1.3333333333333333e+00, 1.3333333333333333e+00, -1.1925695879998881e+00, 
        6.6200910753260603e-17, -1.3240182150652128e-16, -3.3100455376630339e-17, 8.5465008283713980e-17, 
        3.5610420118214158e-17, 0.0000000000000000e+00, -1.4244168047285663e-17, -6.8853037265909633e-01, 
        1.1466334093198018e-16, -6.8853037265909633e-01, -1.1102230246251565e-16, 1.9110556821996684e-17, 
        -2.6390675703848923e-16, -6.8853037265909633e-01, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        3.8221113643993381e-17, -1.1102230246251565e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 4.4133940502173736e-17, -1.1033485125543434e-17, 0.0000000000000000e+00, 
        2.2066970251086868e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, -3.3100455376630302e-17, 
        0.0000000000000000e+00, 1.1033485125543434e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        -1.2740371214664465e-17, 0.0000000000000000e+00, -6.3701856073322325e-18, 0.0000000000000000e+00, 
        2.6666666666666661e+00, -5.5466782398352393e-32, 3.8490017945975055e+00, 3.8490017945975055e+00, 
        -3.8490017945975055e+00, 1.8503717077085943e-16, 7.4014868308343765e-17, 1.3333333333333333e+00, 
        -1.1102230246251564e-16, -1.3333333333333333e+00, -1.3333333333333333e+00, -1.1925695879998881e+00, 
        8.2751138441575754e-17, -1.6550227688315161e-16, -3.6977854932234928e-32, 7.1220840236428317e-17, 
        -4.2732504141856990e-17, 0.0000000000000000e+00, 1.4244168047285663e-17, -6.8853037265909633e-01, 
        1.1466334093198018e-16, -6.8853037265909633e-01, -1.1102230246251565e-16, 1.9110556821996684e-17, 
        -2.6390675703848923e-16, 6.8853037265909633e-01, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        3.8221113643993381e-17, -1.1102230246251565e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 3.8617197939402019e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        -2.2066970251086868e-17, 0.0000000000000000e+00, 5.5511151231257827e-17, 3.3100455376630302e-17, 
        0.0000000000000000e+00, 1.1033485125543434e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.2740371214664465e-17, 9.5552784109983484e-18, 6.3701856073322325e-18, 0.0000000000000000e+00, 
        2.6666666666666661e+00, -4.2732504141857039e-17, 3.8490017945975055e+00, -3.8490017945975055e+00, 
        3.8490017945975055e+00, 2.2204460492503128e-16, -3.7007434154171889e-17, -1.3333333333333333e+00, 
        3.7007434154171876e-17, 1.3333333333333333e+00, -1.3333333333333333e+00, -1.1925695879998881e+00, 
        4.9650683064945446e-17, -1.3240182150652128e-16, -6.6200910753260653e-17, -7.1220840236428317e-17, 
        4.2732504141856990e-17, 7.1220840236428317e-18, 1.4244168047285663e-17, -6.8853037265909633e-01, 
        1.3377389775397689e-16, 6.8853037265909633e-01, 1.1102230246251565e-16, 3.8221113643993381e-17, 
        -2.6390675703848923e-16, -6.8853037265909633e-01, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        3.8221113643993381e-17, -2.6390675703848923e-16, 1.1102230246251565e-16, 2.4671622769447925e-17, 
        0.0000000000000000e+00, -3.8617197939402019e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        2.2066970251086868e-17, 0.0000000000000000e+00, 5.5511151231257827e-17, -3.3100455376630302e-17, 
        0.0000000000000000e+00, 5.5167425627717169e-18, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.2740371214664465e-17, -3.1850928036661163e-18, 0.0000000000000000e+00, 1.5925464018330582e-17, 
        2.6666666666666656e+00, -7.4781882248249794e-17, 3.8490017945975055e+00, -3.8490017945975055e+00, 
        -3.8490017945975055e+00, 2.4671622769447919e-16, -4.9343245538895837e-17, -1.3333333333333333e+00, 
        -1.9737298215558337e-16, -1.3333333333333333e+00, 1.3333333333333333e+00, -1.1925695879998879e+00, 
        6.6200910753260579e-17, -1.3240182150652131e-16, -4.9303806576313238e-32, -7.1220840236428317e-17, 
        -2.8488336094571327e-17, -1.4244168047285663e-17, -3.5610420118214158e-18, -6.8853037265909633e-01, 
        7.6442227287986787e-17, 6.8853037265909633e-01, 1.1102230246251565e-16, -1.2325951644078309e-32, 
        -2.6390675703848923e-16, 6.8853037265909633e-01, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        3.8221113643993381e-17, -1.5288445457597360e-16, 1.1102230246251565e-16, -1.2335811384723962e-17, 
        0.0000000000000000e+00, -3.3100455376630302e-17, -2.2066970251086868e-17, 0.0000000000000000e+00, 
        0.0000000000000000e+00, -5.5511151231257827e-17, 0.0000000000000000e+00, 3.3100455376630302e-17, 
        0.0000000000000000e+00, -1.6550227688315151e-17, 1.1033485125543434e-17, 0.0000000000000000e+00, 
        -3.1850928036661163e-18, 0.0000000000000000e+00, 1.2740371214664465e-17, 3.1850928036661163e-18, 
        2.6666666666666661e+00, 4.2732504141856929e-17, -3.8490017945975055e+00, 3.8490017945975055e+00, 
        3.8490017945975055e+00, -2.0970879354030732e-16, 4.9343245538895837e-17, -1.3333333333333333e+00, 
        1.4802973661668753e-16, -1.3333333333333333e+00, 1.3333333333333333e+00, -1.1925695879998881e+00, 
        6.6200910753260603e-17, -1.6550227688315161e-16, -8.2751138441575803e-17, -7.8342924260071149e-17, 
        -3.5610420118214158e-17, 0.0000000000000000e+00, 1.4244168047285663e-17, 6.8853037265909633e-01, 
        9.5552784109983484e-17, -6.8853037265909633e-01, -1.1102230246251565e-16, 1.9110556821996684e-17, 
        1.1102230246251565e-16, -6.8853037265909633e-01, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        3.8221113643993381e-17, 2.6390675703848923e-16, -1.1102230246251565e-16, 1.2335811384723962e-17, 
        0.0000000000000000e+00, 4.4133940502173736e-17, 1.1033485125543434e-17, 0.0000000000000000e+00, 
        1.1033485125543434e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, -4.4133940502173736e-17, 
        0.0000000000000000e+00, 1.6550227688315151e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.2740371214664465e-17, -6.3701856073322325e-18, 6.3701856073322325e-18, 1.2740371214664465e-17, 
        2.6666666666666656e+00, 6.4098756212785411e-17, -3.8490017945975055e+00, 3.8490017945975055e+00, 
        -3.8490017945975055e+00, -1.3569392523196356e-16, 6.1679056923619811e-17, -1.3333333333333333e+00, 
        -1.3569392523196356e-16, 1.3333333333333333e+00, -1.3333333333333333e+00, -1.1925695879998879e+00, 
        1.4895204919483636e-16, -9.9301366129891004e-17, 1.6550227688315102e-17, -7.8342924260071149e-17, 
        3.5610420118214158e-17, -7.1220840236428317e-18, -7.1220840236428317e-18, 6.8853037265909633e-01, 
        1.2421861934297854e-16, -6.8853037265909633e-01, -1.1102230246251565e-16, 1.9110556821996684e-17, 
        2.6390675703848923e-16, 6.8853037265909633e-01, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        3.8221113643993381e-17, 2.6390675703848923e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 3.8617197939402019e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        -1.1033485125543434e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 3.3100455376630302e-17, 
        0.0000000000000000e+00, 1.6550227688315151e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        -6.3701856073322325e-18, 3.1850928036661163e-18, -6.3701856073322325e-18, 0.0000000000000000e+00, 
        2.6666666666666652e+00, -1.2819751242557102e-16, -3.8490017945975055e+00, -3.8490017945975055e+00, 
        3.8490017945975055e+00, -1.9737298215558340e-16, -8.6350679693067732e-17, 1.3333333333333333e+00, 
        3.7007434154171876e-17, -1.3333333333333333e+00, -1.3333333333333333e+00, -1.1925695879998879e+00, 
        1.8205250457146666e-16, -3.3100455376630388e-17, 6.6200910753260554e-17, 7.1220840236428317e-17, 
        -2.8488336094571327e-17, 7.1220840236428317e-18, -2.8488336094571327e-17, 6.8853037265909633e-01, 
        1.2421861934297854e-16, 6.8853037265909633e-01, 1.1102230246251565e-16, 3.8221113643993381e-17, 
        1.1102230246251565e-16, -6.8853037265909633e-01, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        3.8221113643993381e-17, 2.6390675703848923e-16, 1.1102230246251565e-16, -2.4671622769447925e-17, 
        0.0000000000000000e+00, -1.6550227688315151e-17, 2.2066970251086868e-17, 0.0000000000000000e+00, 
        2.2066970251086868e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, -3.3100455376630302e-17, 
        0.0000000000000000e+00, 1.1033485125543434e-17, 2.2066970251086868e-17, 0.0000000000000000e+00, 
        -2.5480742429328930e-17, -3.1850928036661163e-18, 1.2740371214664465e-17, 0.0000000000000000e+00, 
        2.6666666666666670e+00, -4.2732504141857039e-17, -3.8490017945975055e+00, -3.8490017945975055e+00, 
        -3.8490017945975055e+00, -1.9737298215558340e-16, -4.9343245538895837e-17, 1.3333333333333333e+00, 
        -1.7270135938613544e-16, 1.3333333333333333e+00, 1.3333333333333333e+00, -1.1925695879998879e+00, 
        6.6200910753260579e-17, -1.3240182150652131e-16, -4.9303806576313238e-32, 7.1220840236428317e-17, 
        2.8488336094571327e-17, -1.4244168047285663e-17, 0.0000000000000000e+00, 6.8853037265909633e-01, 
        1.1466334093198018e-16, 6.8853037265909633e-01, 1.1102230246251565e-16, -1.2325951644078309e-32, 
        2.6390675703848923e-16, 6.8853037265909633e-01, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        3.8221113643993381e-17, 2.6390675703848923e-16, 1.1102230246251565e-16, 2.4671622769447925e-17, 
        0.0000000000000000e+00, -4.4133940502173736e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, -5.5511151231257827e-17, -5.5511151231257827e-17, 4.4133940502173736e-17, 
        0.0000000000000000e+00, -1.1033485125543434e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, -1.2740371214664465e-17, -1.2740371214664465e-17, 1.2740371214664465e-17, 
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    }
  }

  gkyl_eval_on_nodes_release(evup);
  gkyl_array_release(distf);
}

void test_2x2v_hyb(int poly_order, int test_func_op)
{
  double lower[] = {-2.0, -2.0, -2.0, -2.0}, upper[] = {2.0, 2.0, 2.0, 2.0};
  int cells[] = {2, 2, 2, 2};
  int ndim = sizeof(lower)/sizeof(lower[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  int cdim = 2; int vdim = 2;
  gkyl_cart_modal_hybrid(&basis, cdim, vdim);

  // projection updater for dist-function
  gkyl_eval_on_nodes *evup;
  if (test_func_op==0)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_4x_quad, NULL);
  else if (test_func_op==1)
    evup = gkyl_eval_on_nodes_new(&grid, &basis, 1, evalFunc_4x_trig, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_eval_on_nodes_advance(evup, 0.0, &arr_range, distf);
  gkyl_grid_sub_array_write(&grid, &arr_range, 0, distf, "ctest_eval_on_nodes_distf_hyb.gkyl");

  if (poly_order == 1) {
    if (test_func_op==0) {
      double fref[] = {
        5.0666666666666664e+01, -1.1547005383792518e+01, -1.1547005383792518e+01, -1.1547005383792518e+01, -1.1547005383792518e+01, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998819e+00, 5.5511151231257827e-16, 
        5.5511151231257827e-16, 5.5511151231257827e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.1925695879998819e+00, 5.5511151231257827e-16, 5.5511151231257827e-16, 5.5511151231257827e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 2.6666666666666664e+01, -6.9282032302755097e+00, -6.9282032302755097e+00, -6.9282032302755105e+00, 
        -2.3094010767585034e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.1925695879998846e+00, 3.3306690738754696e-16, 3.3306690738754696e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, 3.3306690738754696e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 2.6666666666666664e+01, -6.9282032302755097e+00, 
        -6.9282032302755097e+00, -2.3094010767585034e+00, -6.9282032302755105e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, 3.3306690738754696e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, 3.3306690738754696e-16, 
        3.3306690738754696e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.8666666666666664e+01, -2.3094010767585034e+00, -2.3094010767585034e+00, 2.3094010767585034e+00, 2.3094010767585034e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998854e+00, 1.1102230246251565e-16, 
        1.1102230246251565e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.1925695879998854e+00, 1.1102230246251565e-16, 1.1102230246251565e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 2.6666666666666664e+01, -6.9282032302755097e+00, -2.3094010767585034e+00, -6.9282032302755105e+00, 
        -6.9282032302755105e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.1925695879998846e+00, 3.3306690738754696e-16, 1.1102230246251565e-16, 3.3306690738754696e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, 3.3306690738754696e-16, 1.1102230246251565e-16, 3.3306690738754696e-16, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.8666666666666664e+01, -2.3094010767585034e+00, 
        2.3094010767585034e+00, -2.3094010767585034e+00, 2.3094010767585034e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998854e+00, 1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998854e+00, 1.1102230246251565e-16, 
        -1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.8666666666666664e+01, -2.3094010767585034e+00, 2.3094010767585034e+00, 2.3094010767585034e+00, -2.3094010767585034e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998854e+00, 1.1102230246251565e-16, 
        -1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.1925695879998854e+00, 1.1102230246251565e-16, -1.1102230246251565e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 2.6666666666666664e+01, 2.3094010767585034e+00, 6.9282032302755097e+00, 6.9282032302755105e+00, 
        6.9282032302755105e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.1925695879998846e+00, -1.1102230246251565e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, -1.1102230246251565e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 2.6666666666666664e+01, -2.3094010767585034e+00, 
        -6.9282032302755097e+00, -6.9282032302755105e+00, -6.9282032302755105e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, 1.1102230246251565e-16, 3.3306690738754696e-16, 3.3306690738754696e-16, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, 1.1102230246251565e-16, 
        3.3306690738754696e-16, 3.3306690738754696e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.8666666666666664e+01, 2.3094010767585034e+00, -2.3094010767585034e+00, -2.3094010767585034e+00, 2.3094010767585034e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998854e+00, -1.1102230246251565e-16, 
        1.1102230246251565e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.1925695879998854e+00, -1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.8666666666666664e+01, 2.3094010767585034e+00, -2.3094010767585034e+00, 2.3094010767585034e+00, 
        -2.3094010767585034e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.1925695879998854e+00, -1.1102230246251565e-16, 1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998854e+00, -1.1102230246251565e-16, 1.1102230246251565e-16, -1.1102230246251565e-16, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 2.6666666666666664e+01, 6.9282032302755097e+00, 
        2.3094010767585034e+00, 6.9282032302755105e+00, 6.9282032302755105e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, -3.3306690738754696e-16, -1.1102230246251565e-16, -3.3306690738754696e-16, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, -3.3306690738754696e-16, 
        -1.1102230246251565e-16, -3.3306690738754696e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.8666666666666664e+01, 2.3094010767585034e+00, 2.3094010767585034e+00, -2.3094010767585034e+00, -2.3094010767585034e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998854e+00, -1.1102230246251565e-16, 
        -1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.1925695879998854e+00, -1.1102230246251565e-16, -1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 2.6666666666666664e+01, 6.9282032302755097e+00, 6.9282032302755097e+00, 2.3094010767585034e+00, 
        6.9282032302755105e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.1925695879998846e+00, -3.3306690738754696e-16, -3.3306690738754696e-16, -3.3306690738754696e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, -3.3306690738754696e-16, -3.3306690738754696e-16, -1.1102230246251565e-16, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 2.6666666666666664e+01, 6.9282032302755097e+00, 
        6.9282032302755097e+00, 6.9282032302755105e+00, 2.3094010767585034e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, -3.3306690738754696e-16, -3.3306690738754696e-16, -1.1102230246251565e-16, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998846e+00, -3.3306690738754696e-16, 
        -3.3306690738754696e-16, -3.3306690738754696e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        5.0666666666666664e+01, 1.1547005383792518e+01, 1.1547005383792518e+01, 1.1547005383792518e+01, 1.1547005383792518e+01, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.1925695879998819e+00, -5.5511151231257827e-16, 
        -5.5511151231257827e-16, -5.5511151231257827e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        1.1925695879998819e+00, -5.5511151231257827e-16, -5.5511151231257827e-16, -5.5511151231257827e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00,
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    } else if (test_func_op==1) {
      double fref[] = {
        -2.9605947323337506e-16, -6.4098756212785473e-17, 2.3094010767585034e+00, 2.3094010767585034e+00, 2.3094010767585034e+00, -3.7007434154171889e-17, 
        9.8686491077791687e-17, 1.3333333333333333e+00, 2.4671622769447922e-17, 1.3333333333333333e+00, 1.3333333333333333e+00, 0.0000000000000000e+00, 
        -4.2732504141856990e-17, 1.0683126035464246e-17, 5.3415630177321232e-17, -1.8503717077085944e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        -1.1102230246251565e-16, -1.1102230246251565e-16, 1.1033485125543434e-17, 1.6550227688315151e-17, 0.0000000000000000e+00, 6.3701856073322325e-18, 
        8.2751138441575779e-17, -4.7776392054991748e-17, -1.1102230246251565e-16, -3.4580075174528854e-17, -4.4133940502173736e-17, -3.3100455376630302e-17, 
        0.0000000000000000e+00, -9.5552784109983484e-18, -2.5905203907920316e-16, 2.1366252070928492e-17, 2.3094010767585034e+00, 2.3094010767585034e+00, 
        -2.3094010767585034e+00, -1.2335811384723962e-17, 9.8686491077791687e-17, 1.3333333333333333e+00, -2.4671622769447922e-17, -1.3333333333333333e+00, 
        -1.3333333333333333e+00, 0.0000000000000000e+00, 2.4927294082749911e-17, -1.0683126035464246e-17, -5.3415630177321232e-17, 3.0839528461809905e-17, 
        -3.3100455376630302e-17, -1.9110556821996700e-17, -1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 1.6550227688315151e-17, 
        0.0000000000000000e+00, -9.5552784109983484e-18, 4.9650683064945483e-17, -4.7776392054991748e-17, -1.1102230246251565e-16, -3.4580075174528854e-17, 
        -3.3100455376630302e-17, -3.3100455376630302e-17, 0.0000000000000000e+00, -9.5552784109983484e-18, -2.2204460492503131e-16, 7.4781882248249720e-17, 
        2.3094010767585034e+00, -2.3094010767585034e+00, 2.3094010767585034e+00, -1.2335811384723962e-17, -8.6350679693067732e-17, -1.3333333333333333e+00, 
        4.3175339846533866e-17, 1.3333333333333333e+00, -1.3333333333333333e+00, -1.4244168047285663e-17, -4.2732504141856990e-17, -1.0683126035464246e-17, 
        -5.3415630177321232e-17, 3.0839528461809905e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.1102230246251565e-16, -1.1102230246251565e-16, 
        0.0000000000000000e+00, 5.5167425627717169e-18, 0.0000000000000000e+00, 0.0000000000000000e+00, 8.2751138441575779e-17, -6.6886948876988445e-17, 
        -1.1102230246251565e-16, 1.1102230246251565e-16, -2.2066970251086868e-17, 2.2066970251086868e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        -1.4802973661668756e-16, -4.2732504141856984e-17, 2.3094010767585034e+00, -2.3094010767585034e+00, -2.3094010767585034e+00, -3.7007434154171889e-17, 
        -8.6350679693067732e-17, -1.3333333333333333e+00, -4.3175339846533866e-17, -1.3333333333333333e+00, 1.3333333333333333e+00, -1.4244168047285663e-17, 
        1.4244168047285663e-17, 4.2732504141856984e-17, 4.2732504141856984e-17, -2.4671622769447925e-17, -6.6200910753260603e-17, -1.9110556821996700e-17, 
        -1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.2740371214664465e-17, 
        1.6550227688315182e-17, -7.6442227287986800e-17, -1.1102230246251565e-16, 1.1102230246251565e-16, -1.1033485125543434e-17, 2.2066970251086868e-17, 
        0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 2.1366252070928492e-17, -2.3094010767585034e+00, 2.3094010767585034e+00, 
        2.3094010767585034e+00, 2.4671622769447925e-17, 7.4014868308343765e-17, -1.3333333333333333e+00, 2.4671622769447922e-17, -1.3333333333333333e+00, 
        1.3333333333333333e+00, -7.1220840236428317e-18, 3.5610420118214158e-17, -4.2732504141856984e-17, -4.2732504141856984e-17, 3.7007434154171889e-17, 
        0.0000000000000000e+00, -1.9110556821996700e-17, 1.1102230246251565e-16, -1.1102230246251565e-16, -2.2066970251086868e-17, 1.1033485125543434e-17, 
        0.0000000000000000e+00, -6.3701856073322325e-18, 6.6200910753260628e-17, -7.6442227287986800e-17, 1.1102230246251565e-16, -3.4580075174528854e-17, 
        2.2066970251086868e-17, -2.2066970251086868e-17, 0.0000000000000000e+00, -6.3701856073322325e-18, -1.4802973661668753e-16, 9.6148134319178215e-17, 
        -2.3094010767585034e+00, 2.3094010767585034e+00, -2.3094010767585034e+00, 0.0000000000000000e+00, 7.4014868308343765e-17, -1.3333333333333333e+00, 
        -6.1679056923619812e-18, 1.3333333333333333e+00, -1.3333333333333333e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 1.0683126035464246e-17, 
        5.3415630177321232e-17, -3.0839528461809905e-17, 0.0000000000000000e+00, 1.9110556821996700e-17, 1.1102230246251565e-16, 3.4580075174528854e-17, 
        0.0000000000000000e+00, -5.5167425627717169e-18, 0.0000000000000000e+00, 1.2740371214664465e-17, 6.6200910753260628e-17, -6.6886948876988445e-17, 
        1.1102230246251565e-16, -3.4580075174528854e-17, 2.2066970251086868e-17, -2.2066970251086868e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 
        -1.1102230246251565e-16, 8.5465008283713968e-17, -2.3094010767585034e+00, -2.3094010767585034e+00, 2.3094010767585034e+00, -4.9343245538895844e-17, 
        -4.9343245538895844e-17, 1.3333333333333333e+00, 2.4671622769447922e-17, -1.3333333333333333e+00, -1.3333333333333333e+00, 2.8488336094571327e-17, 
        3.5610420118214158e-17, 2.1366252070928492e-17, 4.2732504141856984e-17, -3.0839528461809905e-17, -3.3100455376630302e-17, -1.9110556821996700e-17, 
        1.1102230246251565e-16, -1.1102230246251565e-16, -1.1033485125543434e-17, -1.1033485125543434e-17, 0.0000000000000000e+00, -6.3701856073322325e-18, 
        3.3100455376630326e-17, -2.8665835232995051e-17, 1.1102230246251565e-16, 1.1102230246251565e-16, 0.0000000000000000e+00, 3.3100455376630302e-17, 
        0.0000000000000000e+00, -1.2740371214664465e-17, 1.1102230246251565e-16, 8.5465008283713968e-17, -2.3094010767585034e+00, -2.3094010767585034e+00, 
        -2.3094010767585034e+00, -4.9343245538895844e-17, -4.9343245538895844e-17, 1.3333333333333333e+00, 0.0000000000000000e+00, 1.3333333333333333e+00, 
        1.3333333333333333e+00, 2.8488336094571327e-17, -1.4244168047285663e-17, -4.2732504141856984e-17, -4.2732504141856984e-17, 2.4671622769447925e-17, 
        -9.9301366129890905e-17, 1.9110556821996700e-17, 1.1102230246251565e-16, 3.4580075174528854e-17, -2.2066970251086868e-17, 0.0000000000000000e+00, 
        0.0000000000000000e+00, 0.0000000000000000e+00, -3.3100455376630277e-17, -2.8665835232995051e-17, 1.1102230246251565e-16, 1.1102230246251565e-16, 
        0.0000000000000000e+00, 3.3100455376630302e-17, 0.0000000000000000e+00, -1.2740371214664465e-17, 7.4014868308343773e-16, 3.4186003313485587e-16, 
        2.3094010767585043e+00, 2.3094010767585034e+00, 2.3094010767585034e+00, 2.4671622769447919e-16, 1.2335811384723962e-16, 1.3333333333333333e+00, 
        1.2335811384723962e-16, 1.3333333333333333e+00, 1.3333333333333333e+00, 4.2732504141856990e-17, 5.6976672189142654e-17, 0.0000000000000000e+00, 
        -4.2732504141856984e-17, -2.4671622769447925e-17, -6.6200910753260677e-17, -3.8221113643993412e-17, -2.6390675703848923e-16, -1.1102230246251565e-16, 
        -2.2066970251086868e-17, -5.5167425627717169e-17, 0.0000000000000000e+00, -1.9110556821996697e-17, -3.3100455376630363e-17, -3.8221113643993412e-17, 
        -2.6390675703848923e-16, -2.6390675703848923e-16, -2.2066970251086868e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.9110556821996697e-17, 
        7.4014868308343773e-16, 3.4186003313485587e-16, 2.3094010767585043e+00, 2.3094010767585034e+00, -2.3094010767585034e+00, 2.4671622769447919e-16, 
        1.2335811384723962e-16, 1.3333333333333333e+00, -1.2335811384723962e-16, -1.3333333333333333e+00, -1.3333333333333333e+00, 4.2732504141856990e-17, 
        -7.1220840236428317e-17, 0.0000000000000000e+00, 3.2049378106392736e-17, 2.4671622769447925e-17, -6.6200910753260677e-17, -1.9110556821996712e-17, 
        -2.6390675703848923e-16, 1.1102230246251565e-16, -2.2066970251086868e-17, 6.6200910753260603e-17, 0.0000000000000000e+00, 1.2740371214664465e-17, 
        -3.3100455376630363e-17, -3.8221113643993412e-17, -2.6390675703848923e-16, -2.6390675703848923e-16, -2.2066970251086868e-17, 0.0000000000000000e+00, 
        0.0000000000000000e+00, -1.9110556821996697e-17, 8.1416355139178143e-16, 4.0595878934764128e-16, 2.3094010767585043e+00, -2.3094010767585034e+00, 
        2.3094010767585034e+00, 2.4671622769447919e-16, -1.4802973661668753e-16, -1.3333333333333333e+00, 1.2335811384723962e-16, 1.3333333333333333e+00, 
        -1.3333333333333333e+00, -7.1220840236428317e-17, 5.6976672189142654e-17, 4.2732504141856984e-17, 3.2049378106392736e-17, 2.4671622769447925e-17, 
        -6.6200910753260677e-17, -3.8221113643993412e-17, -2.6390675703848923e-16, -1.1102230246251565e-16, -2.2066970251086868e-17, -3.3100455376630302e-17, 
        0.0000000000000000e+00, -6.3701856073322325e-18, -3.3100455376630363e-17, -5.7331670465990115e-17, -2.6390675703848923e-16, 1.1102230246251565e-16, 
        -4.4133940502173736e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, -6.3701856073322325e-18, 6.4763009769800795e-16, 3.2049378106392740e-16, 
        2.3094010767585043e+00, -2.3094010767585034e+00, -2.3094010767585034e+00, 2.4671622769447919e-16, -1.4802973661668753e-16, -1.3333333333333333e+00, 
        -1.2335811384723962e-16, -1.3333333333333333e+00, 1.3333333333333333e+00, -7.1220840236428317e-17, -7.1220840236428317e-17, 0.0000000000000000e+00, 
        -4.2732504141856984e-17, -1.8503717077085944e-17, -4.9650683064945526e-17, -9.5552784109983577e-18, -2.6390675703848923e-16, 1.1102230246251565e-16, 
        -2.2066970251086868e-17, 6.6200910753260603e-17, 0.0000000000000000e+00, 1.9110556821996697e-17, -1.6550227688315212e-17, -4.7776392054991761e-17, 
        -2.6390675703848923e-16, 1.1102230246251565e-16, -4.4133940502173736e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, -6.3701856073322325e-18, 
        5.9211894646675012e-16, 2.7776127692207036e-16, -2.3094010767585043e+00, 2.3094010767585034e+00, 2.3094010767585034e+00, -2.4671622769447919e-16, 
        1.3569392523196359e-16, -1.3333333333333333e+00, 1.2335811384723962e-16, -1.3333333333333333e+00, 1.3333333333333333e+00, -7.8342924260071149e-17, 
        -4.2732504141856990e-17, -8.5465008283713968e-17, 2.1366252070928492e-17, 1.2335811384723962e-17, -9.9301366129890991e-17, -1.9110556821996712e-17, 
        2.6390675703848923e-16, -1.1102230246251565e-16, 1.1033485125543434e-17, -3.3100455376630302e-17, 0.0000000000000000e+00, 1.9110556821996697e-17, 
        -8.2751138441575828e-17, -4.7776392054991761e-17, 2.6390675703848923e-16, -1.1102230246251565e-16, 0.0000000000000000e+00, 2.2066970251086868e-17, 
        0.0000000000000000e+00, -3.1850928036661163e-18, 5.5511151231257817e-16, 2.7776127692207036e-16, -2.3094010767585043e+00, 2.3094010767585034e+00, 
        -2.3094010767585034e+00, -2.2204460492503128e-16, 1.3569392523196359e-16, -1.3333333333333333e+00, -1.2335811384723962e-16, 1.3333333333333333e+00, 
        -1.3333333333333333e+00, -7.8342924260071149e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, -1.2335811384723962e-17, 
        -3.3100455376630376e-17, -1.5407439555097887e-32, 2.6390675703848923e-16, 1.1102230246251565e-16, 2.2066970251086868e-17, 3.3100455376630302e-17, 
        0.0000000000000000e+00, -2.5480742429328930e-17, -4.9650683064945533e-17, -2.8665835232995064e-17, 2.6390675703848923e-16, -1.1102230246251565e-16, 
        0.0000000000000000e+00, 2.2066970251086868e-17, 0.0000000000000000e+00, -3.1850928036661163e-18, 6.6613381477509392e-16, 3.8459253727671281e-16, 
        -2.3094010767585043e+00, -2.3094010767585034e+00, 2.3094010767585034e+00, -2.4671622769447919e-16, -1.4802973661668753e-16, 1.3333333333333333e+00, 
        1.2335811384723962e-16, -1.3333333333333333e+00, -1.3333333333333333e+00, 5.6976672189142654e-17, -5.6976672189142654e-17, 0.0000000000000000e+00, 
        -4.2732504141856984e-17, -2.4671622769447925e-17, -3.3100455376630388e-17, 1.9110556821996681e-17, 2.6390675703848923e-16, -1.1102230246251565e-16, 
        2.2066970251086868e-17, -1.1033485125543434e-17, 0.0000000000000000e+00, 1.9110556821996697e-17, -7.3955709864469857e-32, -1.9110556821996718e-17, 
        2.6390675703848923e-16, 1.1102230246251565e-16, 2.2066970251086868e-17, 0.0000000000000000e+00, 0.0000000000000000e+00, -2.5480742429328930e-17, 
        6.6613381477509392e-16, 3.8459253727671281e-16, -2.3094010767585043e+00, -2.3094010767585034e+00, -2.3094010767585034e+00, -2.4671622769447919e-16, 
        -1.4802973661668753e-16, 1.3333333333333333e+00, -1.2335811384723962e-16, 1.3333333333333333e+00, 1.3333333333333333e+00, 5.6976672189142654e-17, 
        0.0000000000000000e+00, -8.5465008283713968e-17, 4.2732504141856984e-17, 2.4671622769447925e-17, -7.3955709864469857e-32, 1.9110556821996681e-17, 
        2.6390675703848923e-16, 1.1102230246251565e-16, 2.2066970251086868e-17, 2.2066970251086868e-17, 0.0000000000000000e+00, -2.5480742429328930e-17, 
        -7.3955709864469857e-32, -1.9110556821996718e-17, 2.6390675703848923e-16, 1.1102230246251565e-16, 2.2066970251086868e-17, 0.0000000000000000e+00, 
        0.0000000000000000e+00, -2.5480742429328930e-17, 
      };
      for (int i=0; i<arr_range.volume; i++) {
        double *f_p = gkyl_array_fetch(distf,i);
        for (int k=0; k<basis.num_basis; k++)
          TEST_CHECK( gkyl_compare(fref[i*basis.num_basis+k], f_p[k], 1e-12) );
      }
    }
  }
  gkyl_eval_on_nodes_release(evup);
  gkyl_array_release(distf);
}


void test_1x_p1_quad() { test_1x(1, 0); };
void test_1x_p1_trig() { test_1x(1, 1); };
void test_1x1v_hyb_quad() { test_1x1v_hyb(1, 0); };
void test_1x1v_hyb_trig() { test_1x1v_hyb(1, 1); };
void test_2x_p1_quad() { test_2x(1, 0); };
void test_2x_p1_trig() { test_2x(1, 1); };
void test_3x_p1_quad() { test_3x(1, 0); };
void test_3x_p1_trig() { test_3x(1, 1); };
void test_4x_p1_quad() { test_4x(1, 0); };
void test_4x_p1_trig() { test_4x(1, 1); };
void test_2x2v_hyb_quad() { test_2x2v_hyb(1, 0); };
void test_2x2v_hyb_trig() { test_2x2v_hyb(1, 1); };

void test_1x_p2_quad() { test_1x(2, 0); };
void test_1x_p2_trig() { test_1x(2, 1); };
void test_2x_p2_quad() { test_2x(2, 0); };
void test_2x_p2_trig() { test_2x(2, 1); };
void test_3x_p2_quad() { test_3x(2, 0); };
void test_3x_p2_trig() { test_3x(2, 1); };
void test_4x_p2_quad() { test_4x(2, 0); };
void test_4x_p2_trig() { test_4x(2, 1); };

TEST_LIST = {
  { "test_1x_p1_quad", test_1x_p1_quad },
  { "test_1x_p1_trig", test_1x_p1_trig },
  { "test_1x1v_hyb_quad", test_1x1v_hyb_quad },
  { "test_1x1v_hyb_trig", test_1x1v_hyb_trig },
  { "test_2x_p1_quad", test_2x_p1_quad },
  { "test_2x_p1_trig", test_2x_p1_trig },
  { "test_3x_p1_quad", test_3x_p1_quad },
  { "test_3x_p1_trig", test_3x_p1_trig },
  { "test_4x_p1_quad", test_4x_p1_quad },
  { "test_4x_p1_trig", test_4x_p1_trig },
  { "test_2x2v_hyb_quad", test_2x2v_hyb_quad },
  { "test_2x2v_hyb_trig", test_2x2v_hyb_trig },
//
  { "test_1x_p2_quad", test_1x_p2_quad },
  { "test_1x_p2_trig", test_1x_p2_trig },
  { "test_2x_p2_quad", test_2x_p2_quad },
  { "test_2x_p2_trig", test_2x_p2_trig },
  { "test_3x_p2_quad", test_3x_p2_quad },
  { "test_3x_p2_trig", test_3x_p2_trig },
  { "test_4x_p2_quad", test_4x_p2_quad },
  { "test_4x_p2_trig", test_4x_p2_trig },
  { NULL, NULL },
};
