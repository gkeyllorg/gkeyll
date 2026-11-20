// Test gkyl_deflated_fem_gyroaverage, which gyroaverages a 2x (x,z) or 3x
// (x,y,z) field in the planes perpendicular to the magnetic field (along z).
// Checks the response to a plane wave with a magnetic field that varies along
// z, that constants are preserved, and the adjoint identity.
//
#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_deflated_fem_gyroaverage.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#include <math.h>

#define TEST_MAX_MU 4
#define TEST_MAX_PLANES 16

// Species mass and charge used in the tests.
static const double test_mass = 2.0, test_charge = -1.5;

static struct gkyl_array *
mkarr(long nc, long size)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, size);
}

// Parameters of the fields and geometry used in the tests.
struct test_ctx {
  int cdim; // Number of configuration space dimensions.
  double wavenum; // Wavenumber of the plane wave along x.
  double metric[3]; // Reference g^xx, g^xy, g^yy.
  double bmag; // Reference magnetic field.
  double vary_perp; // Relative amplitude of the variation of the geometry in the plane.
  double vary_par; // Relative variation of the magnetic field per unit length along z.
  // Scaling applied to the plane wave on each plane (used for the expected result).
  int num_planes;
  double z_lower, dz;
  double plane_scale[TEST_MAX_PLANES];
};

static double
geo_variation(const struct test_ctx *tc, const double *xn)
{
  double y = tc->cdim > 2 ? xn[1] : 0.0;
  return 1.0 + tc->vary_perp * sin(0.7 * xn[0] + 0.3) * cos(y);
}

static double
bmag_func(const struct test_ctx *tc, const double *xn)
{
  double z = xn[tc->cdim - 1];
  return tc->bmag * (1.0 + tc->vary_par * z) / geo_variation(tc, xn);
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
  fout[0] = bmag_func(ctx, xn);
}

static void
eval_one(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 1.0;
}

static void
eval_plane_wave(double t, const double *xn, double *fout, void *ctx)
{
  // Plane wave along x with an amplitude that is linear in z.
  struct test_ctx *tc = ctx;
  double z = xn[tc->cdim - 1];
  fout[0] = cos(tc->wavenum * xn[0]) * (1.0 + 0.3 * z);
}

static void
eval_plane_wave_scaled(double t, const double *xn, double *fout, void *ctx)
{
  // Plane wave multiplied by the scaling of the plane this node is on.
  struct test_ctx *tc = ctx;
  double z = xn[tc->cdim - 1];
  int plane = (int)lround((z - tc->z_lower) / tc->dz);
  eval_plane_wave(t, xn, fout, ctx);
  fout[0] *= tc->plane_scale[plane];
}

static void
eval_smooth(double t, const double *xn, double *fout, void *ctx)
{
  // Generic smooth field.
  struct test_ctx *tc = ctx;
  double y = tc->cdim > 2 ? xn[1] : 0.0, z = xn[tc->cdim - 1];
  fout[0] = sin(1.3 * xn[0] + 0.4) * cos(2.0 * y - 0.6 * z) + 0.5 * cos(2.9 * xn[0]) * z +
            0.2 * xn[0] * sin(1.7 * z);
}

// Grid, basis, ranges, geometry and updater of a test.
struct test_setup {
  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext;
  struct gkyl_array *gxx, *gxy, *gyy, *bmag;
  struct gkyl_deflated_fem_gyroaverage *gyroavg;
};

static void
eval_on_nodes(struct test_setup *ts, evalf_t func, struct test_ctx *tc, struct gkyl_array *out)
{
  gkyl_eval_on_nodes *evup = gkyl_eval_on_nodes_new(&ts->grid, &ts->basis, 1, func, tc);
  gkyl_eval_on_nodes_advance(evup, 0.0, &ts->local, out);
  gkyl_eval_on_nodes_release(evup);
}

static struct test_setup
test_setup_new(struct test_ctx *tc, bool periodic_x, int num_mu, const double *mu)
{
  struct test_setup ts;
  int cdim = tc->cdim;
  double lower[] = {0.0, 0.0, 0.0}, upper[] = {2.0 * M_PI, 2.0 * M_PI, 2.0 * M_PI};
  int cells[] = {16, 8, 8};
  // The last direction is z.
  lower[cdim - 1] = -1.0;
  upper[cdim - 1] = 1.5;
  cells[cdim - 1] = 4;
  gkyl_rect_grid_init(&ts.grid, cdim, lower, upper, cells);
  gkyl_cart_modal_serendip(&ts.basis, cdim, 1);
  int ghost[] = {1, 1, 1};
  gkyl_create_grid_ranges(&ts.grid, ghost, &ts.local_ext, &ts.local);
  tc->num_planes = cells[cdim - 1] + 1;
  tc->z_lower = lower[cdim - 1];
  tc->dz = ts.grid.dx[cdim - 1];

  ts.gxx = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  ts.gxy = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  ts.gyy = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  ts.bmag = mkarr(ts.basis.num_basis, ts.local_ext.volume);
  eval_on_nodes(&ts, eval_gxx, tc, ts.gxx);
  eval_on_nodes(&ts, eval_gxy, tc, ts.gxy);
  eval_on_nodes(&ts, eval_gyy, tc, ts.gyy);
  eval_on_nodes(&ts, eval_bmag, tc, ts.bmag);

  struct gkyl_deflated_fem_gyroaverage_inp inp = {
    .grid = &ts.grid,
    .basis = &ts.basis,
    .local = &ts.local,
    .is_periodic = {periodic_x, true},
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
  ts.gyroavg = gkyl_deflated_fem_gyroaverage_inew(&inp);
  return ts;
}

static void
test_setup_release(struct test_setup *ts)
{
  gkyl_array_release(ts->gxx);
  gkyl_array_release(ts->gxy);
  gkyl_array_release(ts->gyy);
  gkyl_array_release(ts->bmag);
  gkyl_deflated_fem_gyroaverage_release(ts->gyroavg);
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
check_equal(
  const struct test_setup *ts, const struct gkyl_array *expected, const struct gkyl_array *fout,
  int mu_idx, double eps
)
{
  // Check that the mu_idx-th gyroaverage equals the expected field.
  int num_basis = ts->basis.num_basis;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ts->local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&ts->local, iter.idx);
    const double *expected_d = gkyl_array_cfetch(expected, linidx);
    const double *fout_d = gkyl_array_cfetch(fout, linidx);
    for (int k = 0; k < num_basis; k++) {
      TEST_CHECK(gkyl_compare_double(expected_d[k], fout_d[mu_idx * num_basis + k], eps));
      TEST_MSG(
        "mu_idx=%d cell=%ld coeff=%d: expected %.12e, got %.12e", mu_idx, linidx, k, expected_d[k],
        fout_d[mu_idx * num_basis + k]
      );
    }
  }
}

static double
symbol_1x(double wavenum, double dx, double rho_perp)
{
  // Exact response of the discrete gyroaverage to a plane wave along x on a
  // uniform periodic grid; rho_perp is the gyroradius times sqrt(g^xx).
  int max_alias = 4000;
  double numer = 0.0, denom = 0.0;
  for (int n = -max_alias; n <= max_alias; n++) {
    double arg = 0.5 * (wavenum * dx + 2.0 * M_PI * n);
    double sinc = fabs(arg) > 0.0 ? sin(arg) / arg : 1.0;
    double wgt = pow(sinc, 4);
    numer += wgt * j0((wavenum + 2.0 * M_PI * n / dx) * rho_perp);
    denom += wgt;
  }
  return numer / denom;
}

static void
test_plane_wave(int cdim)
{
  // With a magnetic field that varies along z the gyroradius, and so the
  // response to a plane wave, differs on each plane of constant z. The result
  // is the linear interpolation along z of the response on each plane.
  struct test_ctx tc = {
    .cdim = cdim,
    .wavenum = 2.0,
    .metric = {1.3, 0.6, 1.7},
    .bmag = 1.5,
    .vary_perp = 0.0,
    .vary_par = 0.2,
  };
  double rho_ref[] = {0.0, 0.2, 0.9};
  int num_mu = sizeof(rho_ref) / sizeof(rho_ref[0]);
  double mu[TEST_MAX_MU];
  for (int n = 0; n < num_mu; n++) {
    mu[n] = 0.5 * tc.bmag * pow(test_charge * rho_ref[n], 2) / test_mass;
  }
  struct test_setup ts = test_setup_new(&tc, true, num_mu, mu);
  int num_basis = ts.basis.num_basis;

  struct gkyl_array *fin = mkarr(num_basis, ts.local_ext.volume);
  struct gkyl_array *fout = mkarr(num_mu * num_basis, ts.local_ext.volume);
  struct gkyl_array *expected = mkarr(num_basis, ts.local_ext.volume);
  eval_on_nodes(&ts, eval_plane_wave, &tc, fin);
  gkyl_deflated_fem_gyroaverage_advance(ts.gyroavg, fin, fout);

  for (int n = 0; n < num_mu; n++) {
    for (int p = 0; p < tc.num_planes; p++) {
      double xn[] = {0.0, 0.0, 0.0};
      xn[cdim - 1] = tc.z_lower + p * tc.dz;
      double rho = sqrt(2.0 * test_mass * mu[n] / bmag_func(&tc, xn)) / fabs(test_charge);
      tc.plane_scale[p] = symbol_1x(tc.wavenum, ts.grid.dx[0], rho * sqrt(tc.metric[0]));
    }
    eval_on_nodes(&ts, eval_plane_wave_scaled, &tc, expected);
    check_equal(&ts, expected, fout, n, 1e-10);
  }

  gkyl_array_release(fin);
  gkyl_array_release(fout);
  gkyl_array_release(expected);
  test_setup_release(&ts);
}

static void
test_adjoint(int cdim)
{
  // With a geometry varying in all directions and a non-periodic boundary in
  // x: constants are preserved, and the adjoint identity
  // sum_n <g_n, G_n f> = <sum_n G_n^dagger g_n, f> holds to machine precision
  // for f continuous along z and g_n discontinuous.
  struct test_ctx tc = {
    .cdim = cdim,
    .metric = {1.3, 0.6, 1.7},
    .bmag = 1.5,
    .vary_perp = 0.3,
    .vary_par = 0.2,
  };
  double rho_ref[] = {0.0, 0.15, 0.5, 1.1};
  int num_mu = sizeof(rho_ref) / sizeof(rho_ref[0]);
  double mu[TEST_MAX_MU];
  for (int n = 0; n < num_mu; n++) {
    mu[n] = 0.5 * tc.bmag * pow(test_charge * rho_ref[n], 2) / test_mass;
  }
  struct test_setup ts = test_setup_new(&tc, false, num_mu, mu);
  int num_basis = ts.basis.num_basis;

  struct gkyl_array *field = mkarr(num_basis, ts.local_ext.volume);
  struct gkyl_array *field_gyroavg = mkarr(num_mu * num_basis, ts.local_ext.volume);
  struct gkyl_array *dens = mkarr(num_mu * num_basis, ts.local_ext.volume);
  struct gkyl_array *dens_gyroavg = mkarr(num_basis, ts.local_ext.volume);

  eval_on_nodes(&ts, eval_one, &tc, field);
  gkyl_deflated_fem_gyroaverage_advance(ts.gyroavg, field, field_gyroavg);
  for (int n = 0; n < num_mu; n++) {
    check_equal(&ts, field, field_gyroavg, n, 1e-12);
  }

  eval_on_nodes(&ts, eval_smooth, &tc, field);
  pcg32_random_t rng = gkyl_pcg32_init(false);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ts.local);
  while (gkyl_range_iter_next(&iter)) {
    double *dens_d = gkyl_array_fetch(dens, gkyl_range_idx(&ts.local, iter.idx));
    for (int k = 0; k < num_mu * num_basis; k++) {
      dens_d[k] = 2.0 * gkyl_pcg32_rand_double(&rng) - 1.0;
    }
  }

  gkyl_deflated_fem_gyroaverage_advance(ts.gyroavg, field, field_gyroavg);
  gkyl_deflated_fem_gyroaverage_advance_adjoint(ts.gyroavg, dens, dens_gyroavg);

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

static void
fill_random(const struct test_setup *ts, pcg32_random_t *rng, struct gkyl_array *arr)
{
  // Random DG coefficients in the local range (discontinuous in every direction).
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ts->local);
  while (gkyl_range_iter_next(&iter)) {
    double *arr_d = gkyl_array_fetch(arr, gkyl_range_idx(&ts->local, iter.idx));
    for (int k = 0; k < arr->ncomp; k++) {
      arr_d[k] = 2.0 * gkyl_pcg32_rand_double(rng) - 1.0;
    }
  }
}

static void
test_dg(int cdim)
{
  // Gyroaverage cell by cell along z (advance_dg), with a geometry varying in all directions and
  // a non-periodic boundary in x: constants are preserved, it equals advance for a field
  // continuous along z, and advance_adjoint is its adjoint for a field discontinuous along z:
  // sum_n <g_n, G_n f> = <sum_n G_n^dagger g_n, f> to machine precision.
  struct test_ctx tc = {
    .cdim = cdim,
    .metric = {1.3, 0.6, 1.7},
    .bmag = 1.5,
    .vary_perp = 0.3,
    .vary_par = 0.2,
  };
  double rho_ref[] = {0.0, 0.15, 0.5, 1.1};
  int num_mu = sizeof(rho_ref) / sizeof(rho_ref[0]);
  double mu[TEST_MAX_MU];
  for (int n = 0; n < num_mu; n++) {
    mu[n] = 0.5 * tc.bmag * pow(test_charge * rho_ref[n], 2) / test_mass;
  }
  struct test_setup ts = test_setup_new(&tc, false, num_mu, mu);
  int num_basis = ts.basis.num_basis;

  struct gkyl_array *field = mkarr(num_basis, ts.local_ext.volume);
  struct gkyl_array *field_gyroavg = mkarr(num_mu * num_basis, ts.local_ext.volume);
  struct gkyl_array *field_gyroavg_dg = mkarr(num_mu * num_basis, ts.local_ext.volume);
  struct gkyl_array *dens = mkarr(num_mu * num_basis, ts.local_ext.volume);
  struct gkyl_array *dens_gyroavg = mkarr(num_basis, ts.local_ext.volume);

  eval_on_nodes(&ts, eval_one, &tc, field);
  gkyl_deflated_fem_gyroaverage_advance_dg(ts.gyroavg, field, field_gyroavg_dg);
  for (int n = 0; n < num_mu; n++) {
    check_equal(&ts, field, field_gyroavg_dg, n, 1e-12);
  }

  // A field continuous along z (nodal interpolation of a smooth function).
  eval_on_nodes(&ts, eval_smooth, &tc, field);
  gkyl_deflated_fem_gyroaverage_advance(ts.gyroavg, field, field_gyroavg);
  gkyl_deflated_fem_gyroaverage_advance_dg(ts.gyroavg, field, field_gyroavg_dg);
  // Equal up to round-off relative to the largest coefficient (the values of the field on the z
  // surfaces of neighboring cells differ by round-off).
  double max_val = 0.0, max_diff = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ts.local);
  while (gkyl_range_iter_next(&iter)) {
    long linidx = gkyl_range_idx(&ts.local, iter.idx);
    const double *fwd_d = gkyl_array_cfetch(field_gyroavg, linidx);
    const double *dg_d = gkyl_array_cfetch(field_gyroavg_dg, linidx);
    for (int k = 0; k < num_mu * num_basis; k++) {
      max_val = fmax(max_val, fabs(fwd_d[k]));
      max_diff = fmax(max_diff, fabs(fwd_d[k] - dg_d[k]));
    }
  }
  TEST_CHECK(max_val > 1e-3);
  TEST_CHECK(max_diff < 1e-13 * max_val);
  TEST_MSG("max value %.3e, max difference %.3e", max_val, max_diff);

  // A field discontinuous along z.
  pcg32_random_t rng = gkyl_pcg32_init(false);
  fill_random(&ts, &rng, field);
  fill_random(&ts, &rng, dens);
  gkyl_deflated_fem_gyroaverage_advance_dg(ts.gyroavg, field, field_gyroavg_dg);
  gkyl_deflated_fem_gyroaverage_advance_adjoint(ts.gyroavg, dens, dens_gyroavg);
  double inner_fwd = 0.0;
  for (int n = 0; n < num_mu; n++) {
    inner_fwd += dg_inner(&ts, dens, n * num_basis, field_gyroavg_dg, n * num_basis);
  }
  double inner_adj = dg_inner(&ts, dens_gyroavg, 0, field, 0);
  TEST_CHECK(fabs(inner_fwd) > 1e-3);
  TEST_CHECK(gkyl_compare_double(inner_fwd, inner_adj, 1e-12));
  TEST_MSG("forward %.15e, adjoint %.15e", inner_fwd, inner_adj);

  gkyl_array_release(field);
  gkyl_array_release(field_gyroavg);
  gkyl_array_release(field_gyroavg_dg);
  gkyl_array_release(dens);
  gkyl_array_release(dens_gyroavg);
  test_setup_release(&ts);
}

void
test_2x_dg(void)
{
  test_dg(2);
}

void
test_3x_dg(void)
{
  test_dg(3);
}

void
test_2x_plane_wave(void)
{
  test_plane_wave(2);
}

void
test_3x_plane_wave(void)
{
  test_plane_wave(3);
}

void
test_2x_adjoint(void)
{
  test_adjoint(2);
}

void
test_3x_adjoint(void)
{
  test_adjoint(3);
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
test_cu(int cdim)
{
  // The gyroaverage and its adjoint give the same result on the GPU as on the
  // host, with a geometry varying in all directions (given on the device) and
  // a non-periodic boundary in x.
  struct test_ctx tc = {
    .cdim = cdim,
    .metric = {1.3, 0.6, 1.7},
    .bmag = 1.5,
    .vary_perp = 0.3,
    .vary_par = 0.2,
  };
  double rho_ref[] = {0.0, 0.15, 0.5, 1.1};
  int num_mu = sizeof(rho_ref) / sizeof(rho_ref[0]);
  double mu[TEST_MAX_MU];
  for (int n = 0; n < num_mu; n++) {
    mu[n] = 0.5 * tc.bmag * pow(test_charge * rho_ref[n], 2) / test_mass;
  }
  struct test_setup ts = test_setup_new(&tc, false, num_mu, mu);
  int num_basis = ts.basis.num_basis;
  long size = ts.local_ext.volume;

  struct gkyl_array *gxx_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_basis, size);
  struct gkyl_array *gxy_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_basis, size);
  struct gkyl_array *gyy_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_basis, size);
  struct gkyl_array *bmag_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_basis, size);
  gkyl_array_copy(gxx_dev, ts.gxx);
  gkyl_array_copy(gxy_dev, ts.gxy);
  gkyl_array_copy(gyy_dev, ts.gyy);
  gkyl_array_copy(bmag_dev, ts.bmag);
  struct gkyl_deflated_fem_gyroaverage_inp inp_dev = {
    .grid = &ts.grid,
    .basis = &ts.basis,
    .local = &ts.local,
    .is_periodic = {false, true},
    .gxx = gxx_dev,
    .gxy = gxy_dev,
    .gyy = gyy_dev,
    .bmag = bmag_dev,
    .mass = test_mass,
    .charge = test_charge,
    .num_mu = num_mu,
    .mu = mu,
    .use_gpu = true,
  };
  struct gkyl_deflated_fem_gyroaverage *gyroavg_dev = gkyl_deflated_fem_gyroaverage_inew(&inp_dev);

  struct gkyl_array *field = mkarr(num_basis, size);
  struct gkyl_array *field_gyroavg = mkarr(num_mu * num_basis, size);
  struct gkyl_array *dens = mkarr(num_mu * num_basis, size);
  struct gkyl_array *dens_gyroavg = mkarr(num_basis, size);
  eval_on_nodes(&ts, eval_smooth, &tc, field);
  pcg32_random_t rng = gkyl_pcg32_init(false);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ts.local);
  while (gkyl_range_iter_next(&iter)) {
    double *dens_d = gkyl_array_fetch(dens, gkyl_range_idx(&ts.local, iter.idx));
    for (int k = 0; k < num_mu * num_basis; k++) {
      dens_d[k] = 2.0 * gkyl_pcg32_rand_double(&rng) - 1.0;
    }
  }

  struct gkyl_array *field_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_basis, size);
  struct gkyl_array *field_gyroavg_dev =
    gkyl_array_cu_dev_new(GKYL_DOUBLE, num_mu * num_basis, size);
  struct gkyl_array *dens_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_mu * num_basis, size);
  struct gkyl_array *dens_gyroavg_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, num_basis, size);
  gkyl_array_copy(field_dev, field);
  gkyl_array_copy(dens_dev, dens);

  gkyl_deflated_fem_gyroaverage_advance(ts.gyroavg, field, field_gyroavg);
  gkyl_deflated_fem_gyroaverage_advance_adjoint(ts.gyroavg, dens, dens_gyroavg);
  gkyl_deflated_fem_gyroaverage_advance(gyroavg_dev, field_dev, field_gyroavg_dev);
  gkyl_deflated_fem_gyroaverage_advance_adjoint(gyroavg_dev, dens_dev, dens_gyroavg_dev);

  check_dev_equal(&ts, field_gyroavg, field_gyroavg_dev);
  check_dev_equal(&ts, dens_gyroavg, dens_gyroavg_dev);

  // Gyroaverage cell by cell along z of a field discontinuous along z.
  fill_random(&ts, &rng, field);
  gkyl_array_copy(field_dev, field);
  gkyl_deflated_fem_gyroaverage_advance_dg(ts.gyroavg, field, field_gyroavg);
  gkyl_deflated_fem_gyroaverage_advance_dg(gyroavg_dev, field_dev, field_gyroavg_dev);
  check_dev_equal(&ts, field_gyroavg, field_gyroavg_dev);

  gkyl_array_release(gxx_dev);
  gkyl_array_release(gxy_dev);
  gkyl_array_release(gyy_dev);
  gkyl_array_release(bmag_dev);
  gkyl_array_release(field);
  gkyl_array_release(field_gyroavg);
  gkyl_array_release(dens);
  gkyl_array_release(dens_gyroavg);
  gkyl_array_release(field_dev);
  gkyl_array_release(field_gyroavg_dev);
  gkyl_array_release(dens_dev);
  gkyl_array_release(dens_gyroavg_dev);
  gkyl_deflated_fem_gyroaverage_release(gyroavg_dev);
  test_setup_release(&ts);
}

void
test_2x_cu(void)
{
  test_cu(2);
}

void
test_3x_cu(void)
{
  test_cu(3);
}
#endif

TEST_LIST = {
  {"test_2x_plane_wave", test_2x_plane_wave},
  {"test_3x_plane_wave", test_3x_plane_wave},
  {"test_2x_adjoint", test_2x_adjoint},
  {"test_3x_adjoint", test_3x_adjoint},
  {"test_2x_dg", test_2x_dg},
  {"test_3x_dg", test_3x_dg},
#ifdef GKYL_HAVE_CUDA
  {"test_2x_cu", test_2x_cu},
  {"test_3x_cu", test_3x_cu},
#endif
  {NULL, NULL}
};
