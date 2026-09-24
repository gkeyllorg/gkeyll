#include <acutest.h>

#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <math.h>

void evalFunc(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  fout[0] = x*x;
}

void
test_1()
{
  int poly_order = 1;
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  // projection updater for dist-function
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc, NULL);

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  // left cell
  double *dfl = gkyl_array_fetch(distf, 0);
  TEST_CHECK( gkyl_compare(1.885618083164127, dfl[0], 1e-12) );
  TEST_CHECK( gkyl_compare(-1.632993161855453, dfl[1], 1e-12) );

  // right cell
  double *dfr = gkyl_array_fetch(distf, 1);
  TEST_CHECK( gkyl_compare(1.885618083164127, dfr[0], 1e-12) );
  TEST_CHECK( gkyl_compare(1.632993161855453, dfr[1], 1e-12) );

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_2()
{
  int poly_order = 1;
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = 3,
      .num_ret_vals = 1,
      .eval = evalFunc,
    }
  );

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  // left cell
  double *dfl = gkyl_array_fetch(distf, 0);
  TEST_CHECK( gkyl_compare(1.885618083164127, dfl[0], 1e-12) );
  TEST_CHECK( gkyl_compare(-1.632993161855453, dfl[1], 1e-12) );

  // right cell
  double *dfr = gkyl_array_fetch(distf, 1);
  TEST_CHECK( gkyl_compare(1.885618083164127, dfr[0], 1e-12) );
  TEST_CHECK( gkyl_compare(1.632993161855453, dfr[1], 1e-12) );

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_2_2d()
{
  int poly_order = 1;
  double lower[] = {-2.0,-2.0}, upper[] = {2.0,2.0};
  int cells[] = {2, 2};
  int ndim = sizeof(cells)/sizeof(cells[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = poly_order+1,
      .num_ret_vals = 1,
      .eval = evalFunc,
    }
  );

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  double xval, xc, dx, xlog, basisval, fval, dgval;
  dx = 2.;

  // left cell
  double *dfl = gkyl_array_fetch(distf, 0);
  xval = -2.;
  basisval = 1./pow(sqrt(2.),ndim);
  fval = pow(xval, 2);
  dgval = 2.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfl[0], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfl[0]);

  xc = -1.;
  xlog = (xval-xc)/(dx/2.);
  basisval = (sqrt(3.)/pow(sqrt(2.),ndim))*xlog;
  fval = pow(xval, 2);
  dgval = 2.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfl[1], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfl[1]);

  dgval = 0.;
  TEST_CHECK( gkyl_compare(dgval, dfl[2], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfl[2]);
  TEST_CHECK( gkyl_compare(dgval, dfl[3], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfl[3]);

  // right cell
  double *dfr = gkyl_array_fetch(distf, 2);
  xval = 2.;
  basisval = 1./pow(sqrt(2.),ndim);
  fval = pow(xval, 2);
  dgval = 2.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfr[0], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfr[0]);

  xc = 1.;
  xlog = (xval-xc)/(dx/2.);
  basisval = (sqrt(3.)/pow(sqrt(2.),ndim))*xlog;
  fval = pow(xval, 2);
  dgval = 2.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfr[1], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfr[1]);

  dgval = 0.;
  TEST_CHECK( gkyl_compare(dgval, dfr[2], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfr[2]);
  TEST_CHECK( gkyl_compare(dgval, dfr[3], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfr[3]);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void
test_2_3d()
{
  int poly_order = 1;
  double lower[] = {-2.0,-2.0,-2.0}, upper[] = {2.0,2.0,2.0};
  int cells[] = {2, 2, 2};
  int ndim = sizeof(cells)/sizeof(cells[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = poly_order+1,
      .num_ret_vals = 1,
      .eval = evalFunc,
    }
  );

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  double xval, xc, dx, xlog, basisval, fval, dgval;
  dx = 2.;

  // left cell
  double *dfl = gkyl_array_fetch(distf, 0);
  xval = -2.;
  basisval = 1./pow(sqrt(2.),ndim);
  fval = pow(xval, 2);
  dgval = 4.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfl[0], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfl[0]);

  xc = -1.;
  xlog = (xval-xc)/(dx/2.);
  basisval = (sqrt(3.)/pow(sqrt(2.),ndim))*xlog;
  fval = pow(xval, 2);
  dgval = 4.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfl[1], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfl[1]);

  dgval = 0.;
  TEST_CHECK( gkyl_compare(dgval, dfl[2], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfl[3], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfl[4], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfl[5], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfl[6], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfl[7], 1e-12) );

  // right cell
  double *dfr = gkyl_array_fetch(distf, 4);
  xval = 2.;
  basisval = 1./pow(sqrt(2.),ndim);
  fval = pow(xval, 2);
  dgval = 4.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfr[0], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfr[0]);

  xc = 1.;
  xlog = (xval-xc)/(dx/2.);
  basisval = (sqrt(3.)/pow(sqrt(2.),ndim))*xlog;
  fval = pow(xval, 2);
  dgval = 4.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfr[1], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfr[1]);

  dgval = 0.;
  TEST_CHECK( gkyl_compare(dgval, dfr[2], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfr[3], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfr[4], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfr[5], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfr[6], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfr[7], 1e-12) );

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

void evalFuncP(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = z*z;
}

void
test_3_3d()
{
  int poly_order = 1;
  double lower[] = {-2.0,-2.0,-2.0}, upper[] = {2.0,2.0,2.0};
  int cells[] = {2, 2, 2};
  int ndim = sizeof(cells)/sizeof(cells[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // basis functions
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = poly_order+1,
      .num_ret_vals = 1,
      .eval = evalFuncP,
    }
  );

  // create array range: no ghost-cells
  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // create distribution function
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);

  // project distribution function on basis
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  double xval, xc, dx, xlog, basisval, fval, dgval;
  dx = 2.;

  // left cell
  double *dfl = gkyl_array_fetch(distf, 0);
  xval = -2.;
  basisval = 1./pow(sqrt(2.),ndim);
  fval = pow(xval, 2);
  dgval = 4.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfl[0], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfl[0]);

  xc = -1.;
  xlog = (xval-xc)/(dx/2.);
  basisval = (sqrt(3.)/pow(sqrt(2.),ndim))*xlog;
  fval = pow(xval, 2);
  dgval = 4.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfl[3], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfl[3]);

  dgval = 0.;
  TEST_CHECK( gkyl_compare(dgval, dfl[1], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfl[2], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfl[4], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfl[5], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfl[6], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfl[7], 1e-12) );

  // right cell
  double *dfr = gkyl_array_fetch(distf, 1);
  xval = 2.;
  basisval = 1./pow(sqrt(2.),ndim);
  fval = pow(xval, 2);
  dgval = 4.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfr[0], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfr[0]);

  xc = 1.;
  xlog = (xval-xc)/(dx/2.);
  basisval = (sqrt(3.)/pow(sqrt(2.),ndim))*xlog;
  fval = pow(xval, 2);
  dgval = 4.*fval*basisval;
  TEST_CHECK( gkyl_compare(dgval, dfr[3], 1e-12) );
  TEST_MSG("Expected: %.13e | Produced: %.13e", dgval, dfr[3]);

  dgval = 0.;
  TEST_CHECK( gkyl_compare(dgval, dfr[1], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfr[2], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfr[4], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfr[5], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfr[6], 1e-12) );
  TEST_CHECK( gkyl_compare(dgval, dfr[7], 1e-12) );

  gkyl_proj_on_basis_release(projDistf);
  gkyl_array_release(distf);
}

TEST_LIST = {
  { "test_1", test_1 },
  { "test_2", test_2 },  
  { "test_2_2d", test_2_2d },  
  { "test_2_3d", test_2_3d },  
  { "test_3_3d", test_3_3d },  
  { NULL, NULL },
};
