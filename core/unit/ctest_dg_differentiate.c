#include <acutest.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_differentiate.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_util.h>

// Allocate array (filled with zeros).
static struct gkyl_array*
mkarr(bool on_gpu, long nc, long size)
{
  struct gkyl_array* a;
  if (on_gpu)
    a = gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size);
  else
    a = gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

struct fin_ctx {
  double exp_c[20]; // Expansion coefficients in test function.
};

void fin_1x_func(double t, const double *xn, double* restrict fout, void *ctx)
{

  struct fin_ctx *params = ctx;
  const double *a = params->exp_c;

  double x = xn[0];
  fout[0] = a[0] + a[1]*x;
}

void
test_dg_differentiate_1x(int poly_order, bool use_gpu)
{
  double lower[] = {-M_PI}, upper[] = {M_PI};
  int cells[] = {6};

  // Parameters for projected function.
  struct fin_ctx inp_params = {0};
  inp_params.exp_c[0] = 2.0;
  inp_params.exp_c[1] = 0.3;

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Basis.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Ranges
  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project the input function.
  struct gkyl_proj_on_basis *proj_fin = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, fin_1x_func, &inp_params);

  // Input field array.
  struct gkyl_array *fin = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fin_ho = use_gpu? mkarr(false, fin->ncomp, fin->size)
                                     : gkyl_array_acquire(fin);

  gkyl_proj_on_basis_advance(proj_fin, 0.0, &local, fin_ho);
  gkyl_array_copy(fin, fin_ho);

  struct gkyl_array *derf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *derf_ho = use_gpu? mkarr(false, derf->ncomp, derf->size)
                                     : gkyl_array_acquire(derf);

  // Differentiate input field.
  int diff_dir = 0;
  int diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);
    int m;
    double ref_val;

    m = 0;
    ref_val = inp_params.exp_c[1]*pow(sqrt(2.0),ndim);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, derf_c[m], ref_val);

    m = 1;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], m, derf_c[m], ref_val);
  }

  gkyl_proj_on_basis_release(proj_fin);
  gkyl_array_release(fin);
  gkyl_array_release(fin_ho);
}

void fin_2x_func(double t, const double *xn, double* restrict fout, void *ctx)
{

  struct fin_ctx *params = ctx;
  const double *a = params->exp_c;

  double x = xn[0], y = xn[1];
  fout[0] = a[0] + a[1]*x + a[2]*y + a[3]*x*y;
}

void
test_dg_differentiate_2x(int poly_order, bool use_gpu)
{
  double lower[] = {-M_PI, -2.0}, upper[] = {M_PI, 2.0};
  int cells[] = {6, 4};

  // Parameters for projected function.
  struct fin_ctx inp_params = {0};
  inp_params.exp_c[0] = 2.0;
  inp_params.exp_c[1] = 0.3;
  inp_params.exp_c[2] = 1.7;
  inp_params.exp_c[3] = 1.1;

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Basis.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Ranges
  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project the input function.
  struct gkyl_proj_on_basis *proj_fin = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, fin_2x_func, &inp_params);

  // Input field array.
  struct gkyl_array *fin = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fin_ho = use_gpu? mkarr(false, fin->ncomp, fin->size)
                                     : gkyl_array_acquire(fin);

  gkyl_proj_on_basis_advance(proj_fin, 0.0, &local, fin_ho);
  gkyl_array_copy(fin, fin_ho);

  struct gkyl_array *derf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *derf_ho = use_gpu? mkarr(false, derf->ncomp, derf->size)
                                     : gkyl_array_acquire(derf);
  
  // Differentiate input field along x.
  int diff_dir = 0;
  int diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);
    int m;
    double ref_val;

    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

    m = 0;
    ref_val = (inp_params.exp_c[1]+xc[1]*inp_params.exp_c[3])*pow(sqrt(2.0),ndim);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 1;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 2;
    ref_val = inp_params.exp_c[3]*grid.dx[1]/sqrt(3.0);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 3;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);
  }

  // Differentiate input field along y.
  diff_dir = 1;
  diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);
    int m;
    double ref_val;

    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

    m = 0;
    ref_val = (inp_params.exp_c[2]+xc[0]*inp_params.exp_c[3])*pow(sqrt(2.0),ndim);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 1;
    ref_val = inp_params.exp_c[3]*grid.dx[0]/sqrt(3.0);
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 2;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);

    m = 3;
    ref_val = 0.0;
    TEST_CHECK( gkyl_compare(derf_c[m], ref_val, 1e-10) );
    TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val);
  }

  gkyl_proj_on_basis_release(proj_fin);
  gkyl_array_release(fin);
  gkyl_array_release(fin_ho);
}

void fin_3x_func(double t, const double *xn, double* restrict fout, void *ctx)
{

  struct fin_ctx *params = ctx;
  const double *a = params->exp_c;

  double x = xn[0], y = xn[1], z = xn[2];
  fout[0] = a[0] + a[1]*x + a[2]*y + a[3]*x*y + a[4]*z + a[5]*x*z + a[6]*y*z + a[7]*x*y*z;
   
}

void
test_dg_differentiate_3x(int poly_order, bool use_gpu)
{
  double lower[] = {-M_PI, -2.0, 1.0}, upper[] = {M_PI, 2.0, 3.5};
  int cells[] = {6, 4, 4};

  // Parameters for projected function.
  struct fin_ctx inp_params = {0};
  inp_params.exp_c[0] = 2.0;
  inp_params.exp_c[1] = 0.3;
  inp_params.exp_c[2] = 1.7;
  inp_params.exp_c[3] = 1.1;
  inp_params.exp_c[4] = 6.2;
  inp_params.exp_c[5] = 3.0;
  inp_params.exp_c[6] = 4.1;
  inp_params.exp_c[7] = 0.5;

  int ndim = sizeof(lower)/sizeof(lower[0]);

  // Basis.
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, ndim, poly_order);

  // Grid.
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, ndim, lower, upper, cells);

  // Ranges
  int ghost[GKYL_MAX_DIM] = {0};
  for (int d=0; d<ndim; d++) ghost[d] = 1;
  struct gkyl_range local, local_ext; // local, local-ext phase-space ranges
  gkyl_create_grid_ranges(&grid, ghost, &local_ext, &local);

  // Project the input function.
  struct gkyl_proj_on_basis *proj_fin = gkyl_proj_on_basis_new(&grid, &basis, poly_order+1, 1, fin_3x_func, &inp_params);

  // Input field array.
  struct gkyl_array *fin = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *fin_ho = use_gpu? mkarr(false, fin->ncomp, fin->size)
                                     : gkyl_array_acquire(fin);

  gkyl_proj_on_basis_advance(proj_fin, 0.0, &local, fin_ho);
  gkyl_array_copy(fin, fin_ho);

  struct gkyl_array *derf = mkarr(use_gpu, basis.num_basis, local_ext.volume);
  struct gkyl_array *derf_ho = use_gpu? mkarr(false, derf->ncomp, derf->size)
                                     : gkyl_array_acquire(derf);
  
  int diff_dir;
  int diff_order;
  struct gkyl_range_iter iter;

  // Differentiate input field along x.
  diff_dir = 0;
  diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);

    const double *dx = grid.dx;
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

    const double *a = inp_params.exp_c;

    const double ref_val[] = {
      pow(sqrt(2),3)*xc[1]*xc[2]*a[7]+pow(sqrt(2),3)*xc[2]*a[5]+pow(sqrt(2),3)*xc[1]*a[3]+pow(sqrt(2),3)*a[1],
      0,
      (sqrt(2)*sqrt(3)*dx[1]*xc[2]*a[7]+sqrt(2)*sqrt(3)*dx[1]*a[3])/3,
      (sqrt(2)*sqrt(3)*xc[1]*dx[2]*a[7]+sqrt(2)*sqrt(3)*dx[2]*a[5])/3,
      0,
      0,
      (dx[1]*dx[2]*a[7])/(3*sqrt(2)),
      0
    };

    for (int m=0; m<basis.num_basis; m++) {
      TEST_CHECK( gkyl_compare(derf_c[m], ref_val[m], 1e-10) );
      TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val[m]);
    }
  }
  
  // Differentiate input field along y.
  diff_dir = 1;
  diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);

    const double *dx = grid.dx;
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

    const double *a = inp_params.exp_c;

    const double ref_val[] = {
      pow(sqrt(2),3)*xc[0]*xc[2]*a[7]+pow(sqrt(2),3)*xc[2]*a[6]+pow(sqrt(2),3)*xc[0]*a[3]+pow(sqrt(2),3)*a[2],
      (sqrt(2)*sqrt(3)*dx[0]*xc[2]*a[7]+sqrt(2)*sqrt(3)*dx[0]*a[3])/3,
      0,
      (sqrt(2)*sqrt(3)*xc[0]*dx[2]*a[7]+sqrt(2)*sqrt(3)*dx[2]*a[6])/3,
      0,
      (dx[0]*dx[2]*a[7])/(3*sqrt(2)),
      0,
      0
    };

    for (int m=0; m<basis.num_basis; m++) {
      TEST_CHECK( gkyl_compare(derf_c[m], ref_val[m], 1e-10) );
      TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val[m]);
    }
  }
  
  // Differentiate input field along x.
  diff_dir = 2;
  diff_order = 1;
  gkyl_dg_differentiate_op_local_range(&basis, diff_dir, diff_order, grid.dx[diff_dir], 0, derf, 0, fin, &local);

  // Check results.
  gkyl_array_copy(derf_ho, derf);
  gkyl_range_iter_init(&iter, &local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&local, iter.idx);
    const double *derf_c = gkyl_array_cfetch(derf_ho, linidx);

    const double *dx = grid.dx;
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc);

    const double *a = inp_params.exp_c;

    const double ref_val[] = {
      pow(sqrt(2),3)*xc[0]*xc[1]*a[7]+pow(sqrt(2),3)*xc[1]*a[6]+pow(sqrt(2),3)*xc[0]*a[5]+pow(sqrt(2),3)*a[4],
      (sqrt(2)*sqrt(3)*dx[0]*xc[1]*a[7]+sqrt(2)*sqrt(3)*dx[0]*a[5])/3,
      (sqrt(2)*sqrt(3)*xc[0]*dx[1]*a[7]+sqrt(2)*sqrt(3)*dx[1]*a[6])/3,
      0,
      (dx[0]*dx[1]*a[7])/(3*sqrt(2)),
      0,
      0,
      0
    };

    for (int m=0; m<basis.num_basis; m++) {
      TEST_CHECK( gkyl_compare(derf_c[m], ref_val[m], 1e-10) );
      TEST_MSG( "idx=%d,%d | m=%d | Got: %.13e | Expected: %.13e\n", iter.idx[0], iter.idx[1], m, derf_c[m], ref_val[m]);
    }
  }

  gkyl_proj_on_basis_release(proj_fin);
  gkyl_array_release(fin);
  gkyl_array_release(fin_ho);
}

void test_dg_differentiate_1x_p1_ho() {
  test_dg_differentiate_1x(1, false);
}

void test_dg_differentiate_2x_p1_ho() {
  test_dg_differentiate_2x(1, false);
}

void test_dg_differentiate_3x_p1_ho() {
  test_dg_differentiate_3x(1, false);
}

#ifdef GKYL_HAVE_CUDA
void test_dg_differentiate_1x_p1_cu() {
  test_dg_differentiate_1x(1, true);
}

void test_dg_differentiate_2x_p1_cu() {
  test_dg_differentiate_2x(1, true);
}

void test_dg_differentiate_3x_p1_cu() {
  test_dg_differentiate_3x(1, true);
}
#endif

TEST_LIST = {
  { "test_dg_differentiate_1x_p1_ho", test_dg_differentiate_1x_p1_ho },
  { "test_dg_differentiate_2x_p1_ho", test_dg_differentiate_2x_p1_ho },
  { "test_dg_differentiate_3x_p1_ho", test_dg_differentiate_3x_p1_ho },
#ifdef GKYL_HAVE_CUDA
  { "test_dg_differentiate_1x_p1_cu", test_dg_differentiate_1x_p1_cu },
  { "test_dg_differentiate_2x_p1_cu", test_dg_differentiate_2x_p1_cu },
  { "test_dg_differentiate_3x_p1_cu", test_dg_differentiate_3x_p1_cu },
#endif
  { NULL, NULL },
};
