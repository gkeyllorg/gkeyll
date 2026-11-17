#include <acutest.h>

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_gaussian_filter.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct affine_ctx {
  double c0;
  double c1;
};

// allocate array (filled with zeros)
static struct gkyl_array *
mkarr(bool use_gpu, long nc, long size)
{
  struct gkyl_array *a = use_gpu ? gkyl_array_cu_dev_new(GKYL_DOUBLE, nc, size) :
                                   gkyl_array_new(GKYL_DOUBLE, nc, size);
  return a;
}

static void
affine_eval(double t, const double *xn, double *restrict fout, void *ctx)
{
  const struct affine_ctx *p = ctx;
  fout[0] = p->c0 + p->c1 * xn[0];
}

static void
accumulate_moments(double x0, double sigma, double dx, double xc, double *i0, double *i1)
{
  const double x_lower = xc - 0.5 * dx;
  const double x_upper = xc + 0.5 * dx;

  const double t_lower = (x_lower - x0) / (sqrt(2.0) * sigma);
  const double t_upper = (x_upper - x0) / (sqrt(2.0) * sigma);

  const double erf_term = erf(t_upper) - erf(t_lower);
  const double exp_lower = exp(-t_lower * t_lower);
  const double exp_upper = exp(-t_upper * t_upper);

  // Integral of constant term at quadrature point x0
  const double cell_i0 = sigma * sqrt(M_PI / 2.0) * erf_term;
  // Integral of linear term at quadrature point x0
  const double cell_i1 = x0 * cell_i0 + sigma * sigma * (exp_lower - exp_upper);

  i0[0] += cell_i0;
  i1[0] += cell_i1;
}

static double
expected_affine_value(
  const struct gkyl_rect_grid *grid, int cell_idx, double eta, double sigma, double c0, double c1
)
{
  double xc[GKYL_MAX_DIM];
  int idx_arr[GKYL_MAX_CDIM] = {cell_idx};
  gkyl_rect_grid_cell_center(grid, idx_arr, xc);

  const double dx = grid->dx[0];
  const double x0 = 0.5 * dx * eta + xc[0];

  double i0 = 0.0, i1 = 0.0;
  for (int offset = -1; offset <= 1; ++offset) {
    int nb_idx[GKYL_MAX_CDIM] = {cell_idx + offset};
    double xc_nb[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(grid, nb_idx, xc_nb);
    accumulate_moments(x0, sigma, dx, xc_nb[0], &i0, &i1);
  }

  return c0 + c1 * (i1 / i0);
}

static void
run_affine_gaussian_filter_test(double c0, double c1, double tol, bool use_gpu)
{
  const int cdim = 1;
  const int poly_order = 3;

  double lower[] = {-1.5};
  double upper[] = {1.5};
  int cells[] = {6};

  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);

  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  int ghost[] = {1};
  struct gkyl_range conf_range, conf_range_ext;
  gkyl_create_grid_ranges(&grid, ghost, &conf_range_ext, &conf_range);

  struct gkyl_array *conf_arr = mkarr(use_gpu, basis.num_basis, conf_range_ext.volume);
  struct gkyl_array *conf_arr_ho;
  if (use_gpu) {
    conf_arr_ho = mkarr(false, basis.num_basis, conf_range_ext.volume);
  } else {
    conf_arr_ho = gkyl_array_acquire(conf_arr);
  }

  struct affine_ctx ctx = {c0, c1};
  gkyl_proj_on_basis *proj =
    gkyl_proj_on_basis_new(&grid, &basis, poly_order + 1, 1, affine_eval, &ctx);
  gkyl_proj_on_basis_advance(proj, 0.0, &conf_range_ext, conf_arr_ho);

  if (use_gpu) {
    gkyl_array_copy(conf_arr, conf_arr_ho);
  }

  struct gkyl_dg_gaussian_filter_inp inp = {
    .conf_grid = &grid,
    .conf_basis = &basis,
    .conf_range = &conf_range,
    .conf_range_ext = &conf_range_ext,
    .extend_filter = false,
    .use_gpu = use_gpu,
  };
  gkyl_dg_gaussian_filter *up = gkyl_dg_gaussian_filter_inew(&inp);

  gkyl_dg_gaussian_filter_advance(up, &conf_range, conf_arr);

  // If we are using GPUs, copy the filtered array back to CPU in the conf_arr_ho array.
  if (use_gpu) {
    gkyl_array_copy(conf_arr_ho, conf_arr);
  }

  const int num_quad = poly_order + 1;
  const double *ords = gkyl_gauss_ordinates[num_quad];
  const double sigma = grid.dx[0];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&conf_range, iter.idx);
    const double *modal = gkyl_array_cfetch(conf_arr_ho, loc);

    double xc_target[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, iter.idx, xc_target);

    for (int n = 0; n < num_quad; ++n) {
      double eta_vec[GKYL_MAX_DIM] = {0.0};
      eta_vec[0] = ords[n];

      double phi[basis.num_basis];
      basis.eval(eta_vec, phi);

      double expected_nodal = expected_affine_value(&grid, iter.idx[0], ords[n], sigma, c0, c1);

      double actual_nodal = 0.0;
      for (int m = 0; m < basis.num_basis; ++m) {
        actual_nodal += modal[m] * phi[m];
      }

      TEST_CHECK(fabs(actual_nodal - expected_nodal) < tol);
      TEST_MSG(
        "Expected %.16e, got %.16e (cell %d node %d)", expected_nodal, actual_nodal, iter.idx[0], n
      );
    }
  }

  gkyl_array_release(conf_arr);
  gkyl_array_release(conf_arr_ho);
  gkyl_proj_on_basis_release(proj);
  gkyl_dg_gaussian_filter_release(up);
}

static void
test_gaussian_filter_constant()
{
  run_affine_gaussian_filter_test(1.25, 0.0, 1e-12, false);
}

static void
test_gaussian_filter_linear()
{
  run_affine_gaussian_filter_test(-0.4, 0.8, 1e-7, false);
}

#ifdef GKYL_HAVE_CUDA
static void
test_gaussian_filter_constant_gpu()
{
  run_affine_gaussian_filter_test(1.25, 0.0, 1e-12, true);
}

static void
test_gaussian_filter_linear_gpu()
{
  run_affine_gaussian_filter_test(-0.4, 0.8, 1e-7, true);
}
#endif

TEST_LIST = {
  {"gaussian_filter_constant", test_gaussian_filter_constant},
  {"gaussian_filter_linear", test_gaussian_filter_linear},
#ifdef GKYL_HAVE_CUDA
  {"gaussian_filter_constant_gpu", test_gaussian_filter_constant_gpu},
  {"gaussian_filter_linear_gpu", test_gaussian_filter_linear_gpu},
#endif
  {NULL, NULL}
};
