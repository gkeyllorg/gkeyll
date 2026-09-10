#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
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

// Cuda specific tests
#ifdef GKYL_HAVE_CUDA

// Getters for device function addresses, defined in ctest_proj_on_basis_cu.cu.
// The device functions must match their host counterparts in this file.
evalf_t ctest_proj_on_basis_f_1d_cu_dev_ptr(void);
evalf_t ctest_proj_on_basis_f_2d_2c_cu_dev_ptr(void);
proj_on_basis_c2p_t ctest_proj_on_basis_c2p_1d_cu_dev_ptr(void);

// Context for the 2d two-component function; must match the definition in
// ctest_proj_on_basis_cu.cu.
struct ctest_proj_on_basis_2d_ctx {
  double c0, c1;
};

// Host counterpart of the device function ctest_pob_f_2d_2c.
void evalFunc_2d_2c(double t, const double *xn, double* restrict fout, void *ctx)
{
  struct ctest_proj_on_basis_2d_ctx *tctx = ctx;
  double x = xn[0], y = xn[1];
  fout[0] = tctx->c0 + x*y;
  fout[1] = tctx->c1*x*x + y;
}

// Host counterpart of the device function ctest_pob_c2p_1d.
void c2pFunc_1d(const double *xcomp, double *xphys, void *ctx)
{
  xphys[0] = 0.5*xcomp[0] + 0.1;
}

// Compare a device array against a host reference, coefficient by coefficient.
static void
check_same_dev_ho(const struct gkyl_array *ref_ho, struct gkyl_array *arr_cu)
{
  struct gkyl_array *arr_ho = gkyl_array_new(GKYL_DOUBLE, arr_cu->ncomp, arr_cu->size);
  gkyl_array_copy(arr_ho, arr_cu);

  for (size_t i=0; i<ref_ho->size; ++i) {
    const double *ref_c = gkyl_array_cfetch(ref_ho, i);
    const double *arr_c = gkyl_array_cfetch(arr_ho, i);
    for (size_t k=0; k<ref_ho->ncomp; ++k) {
      TEST_CHECK( gkyl_compare(ref_c[k], arr_c[k], 1e-14) );
      TEST_MSG("cell %zu coeff %zu | Expected: %.13e | Produced: %.13e", i, k, ref_c[k], arr_c[k]);
    }
  }

  gkyl_array_release(arr_ho);
}

void
test_1_cu()
{
  int poly_order = 1;
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {2};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // Reference: project on the host.
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, 1, evalFunc, NULL);
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  // Project the same function on the device.
  gkyl_proj_on_basis *projDistf_cu = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = poly_order+1,
      .num_ret_vals = 1,
      .eval = ctest_proj_on_basis_f_1d_cu_dev_ptr(),
      .ctx = NULL,
      .use_gpu = true,
    }
  );
  struct gkyl_array *distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_proj_on_basis_advance(projDistf_cu, 0.0, &arr_range, distf_cu);

  check_same_dev_ho(distf, distf_cu);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistf_cu);
  gkyl_array_release(distf);
  gkyl_array_release(distf_cu);
}

void
test_2d_2c_cu()
{
  int poly_order = 2;
  double lower[] = {-2.0,-1.0}, upper[] = {2.0,3.0};
  int cells[] = {4, 3};
  int ndim = sizeof(cells)/sizeof(cells[0]);
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  int num_ret_vals = 2;
  struct ctest_proj_on_basis_2d_ctx tctx = { .c0 = 0.5, .c1 = 2.0 };

  // Reference: project on the host.
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_new(&grid, &basis,
    poly_order+1, num_ret_vals, evalFunc_2d_2c, &tctx);
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE,
    num_ret_vals*basis.num_basis, arr_range.volume);
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  // Project the same function on the device; the context must be placed
  // in device memory by the user.
  struct ctest_proj_on_basis_2d_ctx *tctx_cu = gkyl_cu_malloc(sizeof(struct ctest_proj_on_basis_2d_ctx));
  gkyl_cu_memcpy(tctx_cu, &tctx, sizeof(struct ctest_proj_on_basis_2d_ctx), GKYL_CU_MEMCPY_H2D);

  gkyl_proj_on_basis *projDistf_cu = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = poly_order+1,
      .num_ret_vals = num_ret_vals,
      .eval = ctest_proj_on_basis_f_2d_2c_cu_dev_ptr(),
      .ctx = tctx_cu,
      .use_gpu = true,
    }
  );
  struct gkyl_array *distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE,
    num_ret_vals*basis.num_basis, arr_range.volume);
  gkyl_proj_on_basis_advance(projDistf_cu, 0.0, &arr_range, distf_cu);

  check_same_dev_ho(distf, distf_cu);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistf_cu);
  gkyl_array_release(distf);
  gkyl_array_release(distf_cu);
  gkyl_cu_free(tctx_cu);
}

void
test_c2p_1d_cu()
{
  int poly_order = 1;
  double lower[] = {-2.0}, upper[] = {2.0};
  int cells[] = {8};
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 1, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, poly_order);

  int nghost[GKYL_MAX_DIM] = { 0 };
  struct gkyl_range arr_range, arr_ext_range;
  gkyl_create_grid_ranges(&grid, nghost, &arr_ext_range, &arr_range);

  // Reference: project on the host with the host c2p mapping.
  gkyl_proj_on_basis *projDistf = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = poly_order+1,
      .num_ret_vals = 1,
      .eval = evalFunc,
      .c2p_func = c2pFunc_1d,
    }
  );
  struct gkyl_array *distf = gkyl_array_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_proj_on_basis_advance(projDistf, 0.0, &arr_range, distf);

  // Project on the device with the equivalent device c2p mapping.
  gkyl_proj_on_basis *projDistf_cu = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = &grid,
      .basis = &basis,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = poly_order+1,
      .num_ret_vals = 1,
      .eval = ctest_proj_on_basis_f_1d_cu_dev_ptr(),
      .c2p_func = ctest_proj_on_basis_c2p_1d_cu_dev_ptr(),
      .use_gpu = true,
    }
  );
  struct gkyl_array *distf_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, basis.num_basis, arr_range.volume);
  gkyl_proj_on_basis_advance(projDistf_cu, 0.0, &arr_range, distf_cu);

  check_same_dev_ho(distf, distf_cu);

  gkyl_proj_on_basis_release(projDistf);
  gkyl_proj_on_basis_release(projDistf_cu);
  gkyl_array_release(distf);
  gkyl_array_release(distf_cu);
}

#endif

TEST_LIST = {
  { "test_1", test_1 },
  { "test_2", test_2 },
  { "test_2_2d", test_2_2d },
  { "test_2_3d", test_2_3d },
  { "test_3_3d", test_3_3d },
#ifdef GKYL_HAVE_CUDA
  { "test_1_cu", test_1_cu },
  { "test_2d_2c_cu", test_2d_2c_cu },
  { "test_c2p_1d_cu", test_c2p_1d_cu },
#endif
  { NULL, NULL },
};
