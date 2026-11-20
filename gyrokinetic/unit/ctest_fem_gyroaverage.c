// Test gkyl_fem_gyroaverage, which gyroaverages a field in the plane
// perpendicular to the magnetic field with a weak (FEM) discretization of
//   <f>(X) = 1/(2*pi) int_0^{2*pi} f(X + rho(alpha)) dalpha.
// Checks that constants are preserved, the response to a plane wave against
// the exact symbol of the discrete operator and against J_0(k_perp*rho), the
// adjoint identity, and the polarization response built from the operator.
//
#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_fem_gyroaverage.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <math.h>

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Parameters of the fields and geometry used in the tests.
struct test_ctx {
  int ndim; // Number of dimensions.
  double wavenum[2]; // Wavenumber of the plane wave in each direction.
  double metric[3]; // Reference g^xx, g^xy, g^yy.
  double bmag; // Reference magnetic field.
  double vary; // Relative amplitude of the spatial variation of the geometry.
};

// Species mass and charge used in the tests.
static const double test_mass = 2.0, test_charge = -1.5;

static double
mu_from_rho(double rho, double bmag)
{
  // Magnetic moment for which the gyroradius sqrt(2*m*mu/B)/|q| is rho.
  return 0.5 * bmag * pow(test_charge * rho, 2) / test_mass;
}

static void
eval_plane_wave(double t, const double *xn, double *fout, void *ctx)
{
  struct test_ctx *tc = ctx;
  double arg = 0.0;
  for (int d = 0; d < tc->ndim; d++) {
    arg += tc->wavenum[d] * xn[d];
  }
  fout[0] = cos(arg);
}

static void
eval_one(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 1.0;
}

static void
eval_linear(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 0.5 + 0.7 * xn[0];
}

static double
geo_variation(const struct test_ctx *tc, const double *xn)
{
  double y = tc->ndim > 1 ? xn[1] : 0.0;
  return 1.0 + tc->vary * sin(0.7 * xn[0] + 0.3) * cos(y);
}

static void
eval_gxx(double t, const double *xn, double *fout, void *ctx)
{
  struct test_ctx *tc = ctx;
  fout[0] = tc->metric[0] * geo_variation(tc, xn);
}

static void
eval_gxy(double t, const double *xn, double *fout, void *ctx)
{
  struct test_ctx *tc = ctx;
  fout[0] = tc->metric[1] * geo_variation(tc, xn);
}

static void
eval_gyy(double t, const double *xn, double *fout, void *ctx)
{
  struct test_ctx *tc = ctx;
  fout[0] = tc->metric[2] * geo_variation(tc, xn);
}

static void
eval_bmag(double t, const double *xn, double *fout, void *ctx)
{
  struct test_ctx *tc = ctx;
  fout[0] = tc->bmag / geo_variation(tc, xn);
}

// Grid, basis, ranges, geometry and updater of a test.
struct test_setup {
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext;
  bool is_periodic[2];
  struct gkyl_array *gxx, *gxy, *gyy, *bmag;
  struct gkyl_fem_gyroaverage *gyroavg;
};

static void
eval_on_nodes(struct test_setup *ts, evalf_t func, struct test_ctx *tc, struct gkyl_array *out)
{
  gkyl_eval_on_nodes *evup = gkyl_eval_on_nodes_new(&ts->grid, &ts->basis, 1, func, tc);
  gkyl_eval_on_nodes_advance(evup, 0.0, &ts->local, out);
  gkyl_eval_on_nodes_release(evup);
}

static struct test_setup
test_setup_new(
  struct test_ctx *tc, const double *lower, const double *upper, const int *cells,
  const bool *is_periodic, int num_mu, const double *mu
)
{
  struct test_setup ts;
  gkyl_rect_grid_init(&ts.grid, tc->ndim, lower, upper, cells);
  gkyl_cart_modal_serendip(&ts.basis, tc->ndim, 1);
  int ghost[] = {1, 1};
  gkyl_create_grid_ranges(&ts.grid, ghost, &ts.local_ext, &ts.local);
  ts.is_periodic[0] = is_periodic[0];
  ts.is_periodic[1] = tc->ndim > 1 ? is_periodic[1] : false;

  ts.gxx = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  ts.gxy = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  ts.gyy = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  ts.bmag = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  eval_on_nodes(&ts, eval_gxx, tc, ts.gxx);
  eval_on_nodes(&ts, eval_gxy, tc, ts.gxy);
  eval_on_nodes(&ts, eval_gyy, tc, ts.gyy);
  eval_on_nodes(&ts, eval_bmag, tc, ts.bmag);

  struct gkyl_fem_gyroaverage_inp inp = {
    .grid = &ts.grid,
    .basis = &ts.basis,
    .local = &ts.local,
    .is_periodic = {is_periodic[0], tc->ndim > 1 ? is_periodic[1] : false},
    .gxx = ts.gxx,
    .gxy = ts.gxy,
    .gyy = ts.gyy,
    .bmag = ts.bmag,
    .mass = test_mass,
    .charge = test_charge,
    .num_mu = num_mu,
    .mu = mu,
    .use_gpu = false,
  };
  ts.gyroavg = gkyl_fem_gyroaverage_inew(&inp);
  return ts;
}

static void
test_setup_release(struct test_setup *ts)
{
  gkyl_array_release(ts->gxx);
  gkyl_array_release(ts->gxy);
  gkyl_array_release(ts->gyy);
  gkyl_array_release(ts->bmag);
  gkyl_fem_gyroaverage_release(ts->gyroavg);
}

static double
dg_inner(
  const struct test_setup *ts, const struct gkyl_array *f1, int comp1, const struct gkyl_array *f2,
  int comp2
)
{
  // L2 inner product of two DG fields (stored from components comp1 and comp2).
  double inner = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ts->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&ts->local, iter.idx);
    const double *f1_d = gkyl_array_cfetch(f1, linidx);
    const double *f2_d = gkyl_array_cfetch(f2, linidx);
    for (int k = 0; k < ts->basis.num_basis; k++) {
      inner += f1_d[comp1 + k] * f2_d[comp2 + k];
    }
  }
  return inner * ts->grid.cellVolume / pow(2.0, ts->grid.ndim);
}

static void
fill_random(const struct test_setup *ts, pcg32_random_t *rng, struct gkyl_array *out)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ts->local);
  while (gkyl_range_iter_next(&iter)) {
    double *out_d = gkyl_array_fetch(out, gkyl_range_idx(&ts->local, iter.idx));
    for (int k = 0; k < out->ncomp; k++) {
      out_d[k] = 2.0 * gkyl_pcg32_rand_double(rng) - 1.0;
    }
  }
}

static double
hat_weight(double wavenum, double dx, int alias)
{
  // Squared Fourier transform of the linear nodal basis function (over dx^2)
  // at the wavenumber wavenum + 2*pi*alias/dx.
  double arg = 0.5 * (wavenum * dx + 2.0 * M_PI * alias);
  double sinc = fabs(arg) > 0.0 ? sin(arg) / arg : 1.0;
  return pow(sinc, 4);
}

static double
symbol_1x(double wavenum, double dx, double rho_perp)
{
  // Exact response of the discrete gyroaverage to a plane wave on a uniform
  // periodic 1D grid; rho_perp is the gyroradius times sqrt(g^xx).
  int max_alias = 4000;
  double numer = 0.0, denom = 0.0;
  for (int n = -max_alias; n <= max_alias; n++) {
    double wgt = hat_weight(wavenum, dx, n);
    numer += wgt * j0((wavenum + 2.0 * M_PI * n / dx) * rho_perp);
    denom += wgt;
  }
  return numer / denom;
}

static double
symbol_2x(const double *wavenum, const double *dx, const double *metric, double rho)
{
  // Same as symbol_1x on a uniform periodic 2D grid with a constant metric.
  int max_alias = 300;
  double numer = 0.0, denom = 0.0;
  for (int n = -max_alias; n <= max_alias; n++) {
    double wgt_x = hat_weight(wavenum[0], dx[0], n);
    double kx = wavenum[0] + 2.0 * M_PI * n / dx[0];
    for (int m = -max_alias; m <= max_alias; m++) {
      double wgt = wgt_x * hat_weight(wavenum[1], dx[1], m);
      double ky = wavenum[1] + 2.0 * M_PI * m / dx[1];
      double kperp = sqrt(metric[0] * kx * kx + 2.0 * metric[1] * kx * ky + metric[2] * ky * ky);
      numer += wgt * j0(kperp * rho);
      denom += wgt;
    }
  }
  return numer / denom;
}

static void
check_scaled(
  const struct test_setup *ts, const struct gkyl_array *fin, const struct gkyl_array *fout,
  int mu_idx, double scale, double eps
)
{
  // Check that the mu_idx-th gyroaverage equals scale*fin.
  int num_basis = ts->basis.num_basis;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ts->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&ts->local, iter.idx);
    const double *fin_d = gkyl_array_cfetch(fin, linidx);
    const double *fout_d = gkyl_array_cfetch(fout, linidx);
    for (int k = 0; k < num_basis; k++) {
      TEST_CHECK(gkyl_compare_double(scale * fin_d[k], fout_d[mu_idx * num_basis + k], eps));
      TEST_MSG(
        "mu_idx=%d cell=%ld coeff=%d: expected %.12e, got %.12e", mu_idx, linidx, k,
        scale * fin_d[k], fout_d[mu_idx * num_basis + k]
      );
    }
  }
}

static void
test_constant(int ndim, bool periodic_x)
{
  // The gyroaverage of a constant is the same constant, also with a spatially
  // varying geometry and with a non-periodic boundary.
  struct test_ctx tc = {.ndim = ndim, .metric = {1.3, 0.6, 1.7}, .bmag = 1.5, .vary = 0.3};
  double lower[] = {0.0, 0.0}, upper[] = {2.0 * M_PI, 2.0 * M_PI};
  int cells[] = {12, 8};
  bool is_periodic[] = {periodic_x, true};
  double dx = (upper[0] - lower[0]) / cells[0];
  // Gyroradii from much smaller to larger than a cell.
  double rho[] = {0.0, 0.05 * dx, 0.6 * dx, 1.0 * dx, 2.7 * dx};
  int num_mu = sizeof(rho) / sizeof(rho[0]);
  double mu[num_mu];
  for (int n = 0; n < num_mu; n++) {
    mu[n] = mu_from_rho(rho[n], tc.bmag);
  }
  struct test_setup ts = test_setup_new(&tc, lower, upper, cells, is_periodic, num_mu, mu);

  struct gkyl_array *fin = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  struct gkyl_array *fout = mkarr(num_mu * ts.basis.num_basis, ts.local_ext.volume);
  eval_on_nodes(&ts, eval_one, &tc, fin);
  gkyl_fem_gyroaverage_advance(ts.gyroavg, fin, fout);
  for (int n = 0; n < num_mu; n++) {
    check_scaled(&ts, fin, fout, n, 1.0, 1e-12);
  }

  gkyl_array_release(fin);
  gkyl_array_release(fout);
  test_setup_release(&ts);
}

static void
test_linear_bounded(int ndim)
{
  // The gyroaverage of a function that is linear in a non-periodic direction
  // is the same function up to the boundary (the field is extended linearly
  // beyond it), also with a spatially varying geometry and gyroradii larger
  // than a cell.
  struct test_ctx tc = {.ndim = ndim, .metric = {1.3, 0.6, 1.7}, .bmag = 1.5, .vary = 0.3};
  double lower[] = {0.0, 0.0}, upper[] = {2.0 * M_PI, 2.0 * M_PI};
  int cells[] = {12, 8};
  bool is_periodic[] = {false, true};
  double dx = (upper[0] - lower[0]) / cells[0];
  double rho[] = {0.0, 0.05 * dx, 0.6 * dx, 1.0 * dx, 2.7 * dx};
  int num_mu = sizeof(rho) / sizeof(rho[0]);
  double mu[num_mu];
  for (int n = 0; n < num_mu; n++) {
    mu[n] = mu_from_rho(rho[n], tc.bmag);
  }
  struct test_setup ts = test_setup_new(&tc, lower, upper, cells, is_periodic, num_mu, mu);

  struct gkyl_array *fin = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  struct gkyl_array *fout = mkarr(num_mu * ts.basis.num_basis, ts.local_ext.volume);
  eval_on_nodes(&ts, eval_linear, &tc, fin);
  gkyl_fem_gyroaverage_advance(ts.gyroavg, fin, fout);
  for (int n = 0; n < num_mu; n++) {
    check_scaled(&ts, fin, fout, n, 1.0, 1e-11);
  }

  gkyl_array_release(fin);
  gkyl_array_release(fout);
  test_setup_release(&ts);
}

static void
test_symbol_1x(void)
{
  // Response to a plane wave on a uniform periodic grid: the output must be
  // the input times the exact symbol of the discrete operator.
  struct test_ctx tc = {.ndim = 1, .metric = {1.3, 0.0, 1.0}, .bmag = 1.5, .vary = 0.0};
  double lower[] = {0.0}, upper[] = {2.0 * M_PI};
  int cells[] = {16};
  bool is_periodic[] = {true};
  double dx = (upper[0] - lower[0]) / cells[0];
  tc.wavenum[0] = 2.0;

  double rho[] = {0.0, 0.05 * dx, 0.25 * dx, 0.5 * dx, 1.0 * dx, 2.0 * dx, 3.3 * dx};
  int num_mu = sizeof(rho) / sizeof(rho[0]);
  double mu[num_mu];
  for (int n = 0; n < num_mu; n++) {
    mu[n] = mu_from_rho(rho[n], tc.bmag);
  }
  struct test_setup ts = test_setup_new(&tc, lower, upper, cells, is_periodic, num_mu, mu);

  struct gkyl_array *fin = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  struct gkyl_array *fout = mkarr(num_mu * ts.basis.num_basis, ts.local_ext.volume);
  eval_on_nodes(&ts, eval_plane_wave, &tc, fin);
  gkyl_fem_gyroaverage_advance(ts.gyroavg, fin, fout);
  for (int n = 0; n < num_mu; n++) {
    double rho_perp = rho[n] * sqrt(tc.metric[0]);
    check_scaled(&ts, fin, fout, n, symbol_1x(tc.wavenum[0], dx, rho_perp), 1e-10);
  }

  // Without a gyroradius the matrix is the (tridiagonal) mass matrix.
  TEST_CHECK(gkyl_fem_gyroaverage_num_nonzeros(ts.gyroavg, 0) == 3 * cells[0]);

  gkyl_array_release(fin);
  gkyl_array_release(fout);
  test_setup_release(&ts);
}

static double
symbol_error_1x(int num_cells, double wavenum, double rho)
{
  // Difference between the response of the operator to a plane wave on a
  // uniform periodic grid and J_0(k*rho).
  struct test_ctx tc = {.ndim = 1, .metric = {1.0, 0.0, 1.0}, .bmag = 1.0, .vary = 0.0};
  double lower[] = {0.0}, upper[] = {2.0 * M_PI};
  int cells[] = {num_cells};
  bool is_periodic[] = {true};
  tc.wavenum[0] = wavenum;
  double mu[] = {mu_from_rho(rho, tc.bmag)};
  struct test_setup ts = test_setup_new(&tc, lower, upper, cells, is_periodic, 1, mu);

  struct gkyl_array *fin = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  struct gkyl_array *fout = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  eval_on_nodes(&ts, eval_plane_wave, &tc, fin);
  gkyl_fem_gyroaverage_advance(ts.gyroavg, fin, fout);
  double response = dg_inner(&ts, fin, 0, fout, 0) / dg_inner(&ts, fin, 0, fin, 0);

  gkyl_array_release(fin);
  gkyl_array_release(fout);
  test_setup_release(&ts);
  return response - j0(wavenum * rho);
}

static void
test_convergence_1x(void)
{
  // The response to a plane wave converges to J_0(k*rho) at fourth order.
  double wavenum = 2.0, rho = 0.3;
  double err_coarse = fabs(symbol_error_1x(16, wavenum, rho));
  double err_fine = fabs(symbol_error_1x(32, wavenum, rho));
  double order = log2(err_coarse / err_fine);
  TEST_CHECK(order > 3.8);
  TEST_MSG("errors %.3e (16 cells), %.3e (32 cells), order %.3f", err_coarse, err_fine, order);

  // The FLR correction 1-J_0 is also accurate for a gyroradius much smaller
  // than a cell (rho = 0.05*dx here), where sampling the field on the gyro-ring
  // would overestimate it by more than an order of magnitude.
  int num_cells = 16;
  double wavenum_small = 1.0;
  double rho_small = 0.05 * 2.0 * M_PI / num_cells;
  double flr_exact = 1.0 - j0(wavenum_small * rho_small);
  double flr_error = -symbol_error_1x(num_cells, wavenum_small, rho_small);
  TEST_CHECK(fabs(flr_error) < 0.05 * flr_exact);
  TEST_MSG("1-J_0 = %.6e, error %.3e", flr_exact, flr_error);
}

static void
test_symbol_2x(void)
{
  // Response to a plane wave on a uniform periodic 2D grid with a sheared
  // (non-orthogonal) metric, for which the gyro-ring is a tilted ellipse.
  struct test_ctx tc = {.ndim = 2, .metric = {1.3, 0.6, 1.7}, .bmag = 1.5, .vary = 0.0};
  double lower[] = {0.0, 0.0}, upper[] = {2.0 * M_PI, 2.0 * M_PI};
  int cells[] = {16, 12};
  bool is_periodic[] = {true, true};
  double dx[] = {(upper[0] - lower[0]) / cells[0], (upper[1] - lower[1]) / cells[1]};
  tc.wavenum[0] = 2.0;
  tc.wavenum[1] = -1.0;

  double rho[] = {0.0, 0.05 * dx[0], 0.4 * dx[0], 1.0 * dx[0], 2.2 * dx[0]};
  int num_mu = sizeof(rho) / sizeof(rho[0]);
  double mu[num_mu];
  for (int n = 0; n < num_mu; n++) {
    mu[n] = mu_from_rho(rho[n], tc.bmag);
  }
  struct test_setup ts = test_setup_new(&tc, lower, upper, cells, is_periodic, num_mu, mu);

  struct gkyl_array *fin = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  struct gkyl_array *fout = mkarr(num_mu * ts.basis.num_basis, ts.local_ext.volume);
  eval_on_nodes(&ts, eval_plane_wave, &tc, fin);
  gkyl_fem_gyroaverage_advance(ts.gyroavg, fin, fout);

  double kperp = sqrt(
    tc.metric[0] * tc.wavenum[0] * tc.wavenum[0] +
    2.0 * tc.metric[1] * tc.wavenum[0] * tc.wavenum[1] +
    tc.metric[2] * tc.wavenum[1] * tc.wavenum[1]
  );
  for (int n = 0; n < num_mu; n++) {
    double symbol = symbol_2x(tc.wavenum, dx, tc.metric, rho[n]);
    check_scaled(&ts, fin, fout, n, symbol, 1e-9);
    // The discrete symbol is close to the analytic gyroaverage J_0(k_perp*rho).
    TEST_CHECK(fabs(symbol - j0(kperp * rho[n])) < 2e-2);
    TEST_MSG("rho=%g: symbol %.6e, J_0 %.6e", rho[n], symbol, j0(kperp * rho[n]));
  }

  gkyl_array_release(fin);
  gkyl_array_release(fout);
  test_setup_release(&ts);
}

static void
test_adjoint(int ndim)
{
  // The adjoint identity sum_n <g_n, G_n f> = <sum_n G_n^dagger g_n, f> holds
  // to machine precision for discontinuous fields, a spatially varying
  // geometry and a non-periodic boundary.
  struct test_ctx tc = {.ndim = ndim, .metric = {1.3, 0.6, 1.7}, .bmag = 1.5, .vary = 0.3};
  double lower[] = {0.0, 0.0}, upper[] = {2.0 * M_PI, 2.0 * M_PI};
  int cells[] = {12, 8};
  bool is_periodic[] = {false, true};
  double dx = (upper[0] - lower[0]) / cells[0];
  double rho[] = {0.0, 0.3 * dx, 1.1 * dx, 2.4 * dx};
  int num_mu = sizeof(rho) / sizeof(rho[0]);
  double mu[num_mu];
  for (int n = 0; n < num_mu; n++) {
    mu[n] = mu_from_rho(rho[n], tc.bmag);
  }
  struct test_setup ts = test_setup_new(&tc, lower, upper, cells, is_periodic, num_mu, mu);
  int num_basis = ts.basis.num_basis;

  struct gkyl_array *field = mkarr(num_basis, ts.local_ext.volume);
  struct gkyl_array *field_gyroavg = mkarr(num_mu * num_basis, ts.local_ext.volume);
  struct gkyl_array *dens = mkarr(num_mu * num_basis, ts.local_ext.volume);
  struct gkyl_array *dens_gyroavg = mkarr(num_basis, ts.local_ext.volume);
  pcg32_random_t rng = gkyl_pcg32_init(false);
  fill_random(&ts, &rng, field);
  fill_random(&ts, &rng, dens);

  gkyl_fem_gyroaverage_advance(ts.gyroavg, field, field_gyroavg);
  gkyl_fem_gyroaverage_advance_adjoint(ts.gyroavg, dens, dens_gyroavg);

  double inner_fwd = 0.0;
  for (int n = 0; n < num_mu; n++) {
    inner_fwd += dg_inner(&ts, dens, n * num_basis, field_gyroavg, n * num_basis);
  }
  double inner_adj = dg_inner(&ts, dens_gyroavg, 0, field, 0);
  TEST_CHECK(fabs(inner_fwd) > 1e-3);
  TEST_CHECK(gkyl_compare_double(inner_fwd, inner_adj, 1e-12));
  TEST_MSG("forward %.15e, adjoint %.15e", inner_fwd, inner_adj);

  gkyl_array_release(field);
  gkyl_array_release(field_gyroavg);
  gkyl_array_release(dens);
  gkyl_array_release(dens_gyroavg);
  test_setup_release(&ts);
}

static double
gamma0(double b)
{
  // Gamma_0(b) = exp(-b)*I_0(b) = 1/pi int_0^pi exp(b*(cos(theta)-1)) dtheta.
  int num_pts = 2000;
  double sum = 0.0;
  for (int i = 0; i < num_pts; i++) {
    sum += exp(b * (cos(M_PI * (i + 0.5) / num_pts) - 1.0));
  }
  return sum / num_pts;
}

// Operators and buffers to apply the polarization operator.
struct test_polarization {
  struct test_setup *ts;
  int num_mu; // Number of magnetic moment cells.
  double *mu_wgt; // Integral of the Maxwellian over each magnetic moment cell.
  struct gkyl_array *dens_part; // Maxwellian density deposited at the particle position.
  struct gkyl_array *field_gyroavg, *weighted, *deposited; // Buffers.
};

static void
polarization_apply(
  struct test_polarization *pol, const struct gkyl_array *field, struct gkyl_array *out
)
{
  // Apply the polarization operator P defined by the bilinear form
  //   <g, P f> = int F_0(mu) (<g f> - <g><f>) dX dmu,
  // with a Maxwellian F_0 and one gyroaverage per magnetic moment cell:
  // P f = n_part*f - sum_j W_j G_j^dagger G_j f.
  struct test_setup *ts = pol->ts;
  int num_basis = ts->basis.num_basis, num_mu = pol->num_mu;
  gkyl_fem_gyroaverage_advance(ts->gyroavg, field, pol->field_gyroavg);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ts->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&ts->local, iter.idx);
    const double *gyroavg_d = gkyl_array_cfetch(pol->field_gyroavg, linidx);
    double *weighted_d = gkyl_array_fetch(pol->weighted, linidx);
    for (int j = 0; j < num_mu; j++) {
      for (int k = 0; k < num_basis; k++) {
        weighted_d[j * num_basis + k] = pol->mu_wgt[j] * gyroavg_d[j * num_basis + k];
      }
    }
  }
  gkyl_fem_gyroaverage_advance_adjoint(ts->gyroavg, pol->weighted, pol->deposited);

  gkyl_dg_mul_op_range(&ts->basis, 0, out, 0, pol->dens_part, 0, field, &ts->local);
  gkyl_array_accumulate_range(out, -1.0, pol->deposited, &ts->local);
}

static struct test_polarization
polarization_new(
  struct test_setup *ts, int num_mu, const double *mu, const double *mu_edges, double temp_over_bmag
)
{
  struct test_polarization pol = {.ts = ts, .num_mu = num_mu};
  int num_basis = ts->basis.num_basis;

  // W_j = int F_0 dmu over the magnetic moment cell j, with
  // F_0 = exp(-mu/(T/B))/(T/B).
  pol.mu_wgt = gkyl_malloc(sizeof(double[num_mu]));
  for (int j = 0; j < num_mu; j++) {
    pol.mu_wgt[j] = exp(-mu_edges[j] / temp_over_bmag) - exp(-mu_edges[j + 1] / temp_over_bmag);
  }

  pol.dens_part = mkarr(num_basis, ts->local_ext.volume);
  pol.field_gyroavg = mkarr(num_mu * num_basis, ts->local_ext.volume);
  pol.weighted = mkarr(num_mu * num_basis, ts->local_ext.volume);
  pol.deposited = mkarr(num_basis, ts->local_ext.volume);

  // n_part = sum_j W_j G_j^dagger 1, with the adjoint of the
  // gyroaverage that reflects the field about non-periodic boundaries: P then
  // annihilates constants away from those boundaries, and is positive
  // semi-definite for fields that vanish at them (which the gyroaverage
  // extends as odd functions, so that their products are even).
  struct gkyl_array *one = mkarr(num_basis, ts->local_ext.volume);
  eval_on_nodes(ts, eval_one, 0, one);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ts->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&ts->local, iter.idx);
    const double *one_d = gkyl_array_cfetch(one, linidx);
    double *weighted_d = gkyl_array_fetch(pol.weighted, linidx);
    for (int j = 0; j < num_mu; j++) {
      for (int k = 0; k < num_basis; k++) {
        weighted_d[j * num_basis + k] = pol.mu_wgt[j] * one_d[k];
      }
    }
  }
  struct gkyl_fem_gyroaverage_inp refl_inp = {
    .grid = &ts->grid,
    .basis = &ts->basis,
    .local = &ts->local,
    .is_periodic = {ts->is_periodic[0], ts->is_periodic[1]},
    .reflect_bnd = true,
    .gxx = ts->gxx,
    .gxy = ts->gxy,
    .gyy = ts->gyy,
    .bmag = ts->bmag,
    .mass = test_mass,
    .charge = test_charge,
    .num_mu = num_mu,
    .mu = mu,
    .use_gpu = false,
  };
  struct gkyl_fem_gyroaverage *gyroavg_refl = gkyl_fem_gyroaverage_inew(&refl_inp);
  gkyl_fem_gyroaverage_advance_adjoint(gyroavg_refl, pol.weighted, pol.dens_part);
  gkyl_fem_gyroaverage_release(gyroavg_refl);
  gkyl_array_release(one);
  return pol;
}

static void
polarization_release(struct test_polarization *pol)
{
  gkyl_free(pol->mu_wgt);
  gkyl_array_release(pol->dens_part);
  gkyl_array_release(pol->field_gyroavg);
  gkyl_array_release(pol->weighted);
  gkyl_array_release(pol->deposited);
}

static double
symmetric_min_eigenvalue(int num, double *mat, double *max_eig)
{
  // Smallest (returned) and largest eigenvalues of a symmetric matrix,
  // computed with cyclic Jacobi rotations. The matrix is overwritten.
  for (int sweep = 0; sweep < 100; sweep++) {
    double off_norm = 0.0;
    for (int p = 0; p < num; p++) {
      for (int q = p + 1; q < num; q++) {
        off_norm += mat[p * num + q] * mat[p * num + q];
      }
    }
    if (off_norm < 1e-30) {
      break;
    }
    for (int p = 0; p < num; p++) {
      for (int q = p + 1; q < num; q++) {
        if (mat[p * num + q] == 0.0) {
          continue;
        }
        double theta = 0.5 * atan2(2.0 * mat[p * num + q], mat[q * num + q] - mat[p * num + p]);
        double cs = cos(theta), sn = sin(theta);
        for (int i = 0; i < num; i++) {
          double mat_ip = mat[i * num + p], mat_iq = mat[i * num + q];
          mat[i * num + p] = cs * mat_ip - sn * mat_iq;
          mat[i * num + q] = sn * mat_ip + cs * mat_iq;
        }
        for (int i = 0; i < num; i++) {
          double mat_pi = mat[p * num + i], mat_qi = mat[q * num + i];
          mat[p * num + i] = cs * mat_pi - sn * mat_qi;
          mat[q * num + i] = sn * mat_pi + cs * mat_qi;
        }
      }
    }
  }
  double min_eig = mat[0];
  *max_eig = mat[0];
  for (int p = 1; p < num; p++) {
    min_eig = fmin(min_eig, mat[p * num + p]);
    *max_eig = fmax(*max_eig, mat[p * num + p]);
  }
  return min_eig;
}

static void
test_polarization_1x(bool periodic, double vary)
{
  // Polarization response built from the gyroaverage and its adjoint, with
  // magnetic moment cells uniformly spaced in perpendicular velocity, each
  // gyroaveraged at the average magnetic moment of the Maxwellian in the cell.
  struct test_ctx tc = {.ndim = 1, .metric = {1.0, 0.0, 1.0}, .bmag = 1.0, .vary = vary};
  double length = 8.0 * M_PI;
  double lower[] = {0.0}, upper[] = {length};
  int cells[] = {24};
  bool is_periodic[] = {periodic};

  // Temperature such that the thermal gyroradius sqrt(T*m)/(|q|*B) is 1.
  double temp_over_bmag = test_charge * test_charge * tc.bmag / test_mass;
  int num_mu = 16;
  double mu[num_mu], mu_edges[num_mu + 1];
  for (int n = 0; n < num_mu + 1; n++) {
    mu_edges[n] = 16.0 * temp_over_bmag * pow(n / (double)num_mu, 2);
  }
  for (int j = 0; j < num_mu; j++) {
    double cell_len = mu_edges[j + 1] - mu_edges[j];
    double decay = exp(-cell_len / temp_over_bmag);
    mu[j] = mu_edges[j] + temp_over_bmag - cell_len * decay / (1.0 - decay);
  }
  struct test_setup ts = test_setup_new(&tc, lower, upper, cells, is_periodic, num_mu, mu);
  struct test_polarization pol = polarization_new(&ts, num_mu, mu, mu_edges, temp_over_bmag);
  int num_basis = ts.basis.num_basis;

  struct gkyl_array *field = mkarr(num_basis, ts.local_ext.volume);
  struct gkyl_array *pol_field = mkarr(num_basis, ts.local_ext.volume);

  if (periodic) {
    // P annihilates constants, so a uniform potential drives no polarization.
    eval_on_nodes(&ts, eval_one, &tc, field);
    polarization_apply(&pol, field, pol_field);
    double norm_const = sqrt(dg_inner(&ts, pol_field, 0, pol_field, 0) / length);
    TEST_CHECK(norm_const < 1e-12);
    TEST_MSG("norm of P applied to a constant: %.3e", norm_const);
  }

  if (periodic && vary == 0.0) {
    // For a plane wave the response is 1-Gamma_0(b), b = (k*rho_t)^2. Only
    // modes with k*dx <~ 1 are checked: the error grows to a few percent at
    // k*dx = 2 because of the spatial resolution.
    double dens_tot = 0.0;
    for (int j = 0; j < num_mu; j++) {
      dens_tot += pol.mu_wgt[j];
    }
    int modes[] = {1, 2, 4};
    for (int m = 0; m < 3; m++) {
      tc.wavenum[0] = 2.0 * M_PI * modes[m] / length;
      eval_on_nodes(&ts, eval_plane_wave, &tc, field);
      polarization_apply(&pol, field, pol_field);
      double response = dg_inner(&ts, field, 0, pol_field, 0) / dg_inner(&ts, field, 0, field, 0);
      double exact = 1.0 - gamma0(tc.wavenum[0] * tc.wavenum[0]);
      TEST_CHECK(gkyl_compare_double(exact, response / dens_tot, 1e-2));
      TEST_MSG(
        "k*rho_t=%.3f: 1-Gamma_0 = %.6e, got %.6e", tc.wavenum[0], exact, response / dens_tot
      );
    }
  }

  // P restricted to continuous fields is symmetric positive semi-definite. If
  // the boundaries are not periodic this holds for fields that vanish there,
  // like the potential (the linear extension beyond the boundary then is the
  // odd reflection of the field). Build the matrix of P in the nodal basis of
  // the nodes where the field is free, and compute its smallest eigenvalue.
  int first_node = periodic ? 0 : 1;
  int num_nodes = periodic ? cells[0] : cells[0] - 1;
  double node_coords[2];
  ts.basis.node_list(node_coords);
  int lo_node = node_coords[0] < 0.0 ? 0 : 1; // Local index of the lower node of a cell.
  struct gkyl_array **nodal_fields = gkyl_malloc(sizeof(struct gkyl_array *[num_nodes]));
  struct gkyl_array **pol_nodal_fields = gkyl_malloc(sizeof(struct gkyl_array *[num_nodes]));
  for (int i = 0; i < num_nodes; i++) {
    nodal_fields[i] = mkarr(num_basis, ts.local_ext.volume);
    pol_nodal_fields[i] = mkarr(num_basis, ts.local_ext.volume);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &ts.local);
    while (gkyl_range_iter_next(&iter)) {
      int cell = iter.idx[0] - ts.local.lower[0];
      int node_lo = cell, node_up = periodic ? (cell + 1) % cells[0] : cell + 1;
      double nodal[2];
      nodal[lo_node] = node_lo == first_node + i ? 1.0 : 0.0;
      nodal[1 - lo_node] = node_up == first_node + i ? 1.0 : 0.0;
      long linidx = gkyl_range_idx(&ts.local, iter.idx);
      ts.basis.nodal_to_modal(nodal, gkyl_array_fetch(nodal_fields[i], linidx));
    }
    polarization_apply(&pol, nodal_fields[i], pol_nodal_fields[i]);
  }
  double *pol_mat = gkyl_malloc(sizeof(double[num_nodes * num_nodes]));
  double asymmetry = 0.0;
  for (int i = 0; i < num_nodes; i++) {
    for (int j = 0; j < num_nodes; j++) {
      pol_mat[i * num_nodes + j] = dg_inner(&ts, nodal_fields[i], 0, pol_nodal_fields[j], 0);
    }
  }
  for (int i = 0; i < num_nodes; i++) {
    for (int j = 0; j < i; j++) {
      asymmetry = fmax(asymmetry, fabs(pol_mat[i * num_nodes + j] - pol_mat[j * num_nodes + i]));
    }
  }
  double max_eig;
  double min_eig = symmetric_min_eigenvalue(num_nodes, pol_mat, &max_eig);
  TEST_CHECK(asymmetry < 1e-13 * max_eig);
  TEST_MSG("asymmetry %.3e", asymmetry);
  TEST_CHECK(min_eig > -1e-12 * max_eig);
  TEST_MSG("eigenvalues of P between %.3e and %.3e", min_eig, max_eig);

  for (int i = 0; i < num_nodes; i++) {
    gkyl_array_release(nodal_fields[i]);
    gkyl_array_release(pol_nodal_fields[i]);
  }
  gkyl_free(nodal_fields);
  gkyl_free(pol_nodal_fields);
  gkyl_free(pol_mat);
  gkyl_array_release(field);
  gkyl_array_release(pol_field);
  polarization_release(&pol);
  test_setup_release(&ts);
}

void
test_1x_constant_periodic(void)
{
  test_constant(1, true);
}

void
test_1x_constant_bounded(void)
{
  test_constant(1, false);
}

void
test_2x_constant_periodic(void)
{
  test_constant(2, true);
}

void
test_2x_constant_bounded(void)
{
  test_constant(2, false);
}

void
test_1x_linear_bounded(void)
{
  test_linear_bounded(1);
}

void
test_2x_linear_bounded(void)
{
  test_linear_bounded(2);
}

void
test_1x_symbol(void)
{
  test_symbol_1x();
}

void
test_1x_convergence(void)
{
  test_convergence_1x();
}

void
test_2x_symbol(void)
{
  test_symbol_2x();
}

void
test_1x_adjoint(void)
{
  test_adjoint(1);
}

void
test_2x_adjoint(void)
{
  test_adjoint(2);
}

void
test_1x_polarization_uniform(void)
{
  test_polarization_1x(true, 0.0);
}

void
test_1x_polarization_varying(void)
{
  test_polarization_1x(false, 0.3);
}

#ifdef GKYL_HAVE_CUDA
static void
check_dev_equal(
  const struct test_setup *ts, const struct gkyl_array *arr_ho, const struct gkyl_array *arr_dev
)
{
  // Check that an array on the device equals an array on the host.
  struct gkyl_array *arr_dev_ho = mkarr(arr_ho->ncomp, arr_ho->size);
  gkyl_array_copy(arr_dev_ho, arr_dev);

  double max_val = 0.0, max_diff = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ts->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&ts->local, iter.idx);
    const double *arr_ho_d = gkyl_array_cfetch(arr_ho, linidx);
    const double *arr_dev_ho_d = gkyl_array_cfetch(arr_dev_ho, linidx);
    for (int k = 0; k < arr_ho->ncomp; k++) {
      max_val = fmax(max_val, fabs(arr_ho_d[k]));
      max_diff = fmax(max_diff, fabs(arr_ho_d[k] - arr_dev_ho_d[k]));
    }
  }
  TEST_CHECK(max_val > 1e-3);
  TEST_CHECK(max_diff < 1e-12 * max_val);
  TEST_MSG("max value %.3e, max difference %.3e", max_val, max_diff);

  gkyl_array_release(arr_dev_ho);
}

static void
test_cu(int ndim)
{
  // The gyroaverage and its adjoint give the same result on the GPU as on the
  // host, for discontinuous fields, a spatially varying geometry and a
  // non-periodic boundary.
  struct test_ctx tc = {.ndim = ndim, .metric = {1.3, 0.6, 1.7}, .bmag = 1.5, .vary = 0.3};
  double lower[] = {0.0, 0.0}, upper[] = {2.0 * M_PI, 2.0 * M_PI};
  int cells[] = {12, 8};
  bool is_periodic[] = {false, true};
  double dx = (upper[0] - lower[0]) / cells[0];
  double rho[] = {0.0, 0.3 * dx, 1.1 * dx, 2.4 * dx};
  int num_mu = sizeof(rho) / sizeof(rho[0]);
  double mu[num_mu];
  for (int n = 0; n < num_mu; n++) {
    mu[n] = mu_from_rho(rho[n], tc.bmag);
  }
  struct test_setup ts = test_setup_new(&tc, lower, upper, cells, is_periodic, num_mu, mu);
  int num_basis = ts.basis.num_basis;

  struct gkyl_fem_gyroaverage_inp inp_dev = {
    .grid = &ts.grid,
    .basis = &ts.basis,
    .local = &ts.local,
    .is_periodic = {is_periodic[0], ndim > 1 ? is_periodic[1] : false},
    .gxx = ts.gxx,
    .gxy = ts.gxy,
    .gyy = ts.gyy,
    .bmag = ts.bmag,
    .mass = test_mass,
    .charge = test_charge,
    .num_mu = num_mu,
    .mu = mu,
    .use_gpu = true,
  };
  struct gkyl_fem_gyroaverage *gyroavg_dev = gkyl_fem_gyroaverage_inew(&inp_dev);

  long size = ts.local_ext.volume;
  struct gkyl_array *field = mkarr(num_basis, size);
  struct gkyl_array *field_gyroavg = mkarr(num_mu * num_basis, size);
  struct gkyl_array *dens = mkarr(num_mu * num_basis, size);
  struct gkyl_array *dens_gyroavg = mkarr(num_basis, size);
  pcg32_random_t rng = gkyl_pcg32_init(false);
  fill_random(&ts, &rng, field);
  fill_random(&ts, &rng, dens);

  struct gkyl_array *field_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_basis, size);
  struct gkyl_array *field_gyroavg_dev =
    gkyl_array_cu_dev_new(GKYL_DOUBLE, num_mu * num_basis, size);
  struct gkyl_array *dens_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_mu * num_basis, size);
  struct gkyl_array *dens_gyroavg_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_basis, size);
  gkyl_array_copy(field_dev, field);
  gkyl_array_copy(dens_dev, dens);

  gkyl_fem_gyroaverage_advance(ts.gyroavg, field, field_gyroavg);
  gkyl_fem_gyroaverage_advance_adjoint(ts.gyroavg, dens, dens_gyroavg);
  gkyl_fem_gyroaverage_advance(gyroavg_dev, field_dev, field_gyroavg_dev);
  gkyl_fem_gyroaverage_advance_adjoint(gyroavg_dev, dens_dev, dens_gyroavg_dev);

  check_dev_equal(&ts, field_gyroavg, field_gyroavg_dev);
  check_dev_equal(&ts, dens_gyroavg, dens_gyroavg_dev);

  gkyl_array_release(field);
  gkyl_array_release(field_gyroavg);
  gkyl_array_release(dens);
  gkyl_array_release(dens_gyroavg);
  gkyl_array_release(field_dev);
  gkyl_array_release(field_gyroavg_dev);
  gkyl_array_release(dens_dev);
  gkyl_array_release(dens_gyroavg_dev);
  gkyl_fem_gyroaverage_release(gyroavg_dev);
  test_setup_release(&ts);
}

void
test_1x_cu(void)
{
  test_cu(1);
}

void
test_2x_cu(void)
{
  test_cu(2);
}
#endif

TEST_LIST = {
  {"test_1x_constant_periodic", test_1x_constant_periodic},
  {"test_1x_constant_bounded", test_1x_constant_bounded},
  {"test_2x_constant_periodic", test_2x_constant_periodic},
  {"test_2x_constant_bounded", test_2x_constant_bounded},
  {"test_1x_linear_bounded", test_1x_linear_bounded},
  {"test_2x_linear_bounded", test_2x_linear_bounded},
  {"test_1x_symbol", test_1x_symbol},
  {"test_1x_convergence", test_1x_convergence},
  {"test_2x_symbol", test_2x_symbol},
  {"test_1x_adjoint", test_1x_adjoint},
  {"test_2x_adjoint", test_2x_adjoint},
  {"test_1x_polarization_uniform", test_1x_polarization_uniform},
  {"test_1x_polarization_varying", test_1x_polarization_varying},
#ifdef GKYL_HAVE_CUDA
  {"test_1x_cu", test_1x_cu},
  {"test_2x_cu", test_2x_cu},
#endif
  {NULL, NULL}
};
