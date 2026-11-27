// Tests for the updater of DG fluid variables (gkyl_dg_calc_fluid_vars) and the surface kernels
// of the DG Euler equation with the Lax solver:
//   em_source  : the source of a charged fluid coupled explicitly to the field. Without a magnetic
//                field it is the source of an applied acceleration; for uniform fields and a
//                uniform flow it is rho a + rho u x (q/m B) and rho u . a, also for a density that
//                varies in the cell. The CFL rate of a cell becomes the larger of the rate it had
//                and the rate of the source.
//   positivity : the limiter scales the slopes of a cell with an inadmissible vertex density or
//                internal energy until the vertices are at the floor, leaves the cell averages
//                alone, and does not touch an admissible cell.
//   limiter    : a slope along the sweep direction of the momentum component along it that the
//                limiter removes takes with it the higher coefficients of the density, of that
//                momentum component and of the energy (the components of the acoustic waves), in
//                every sweep direction, and leaves those of the tangential momenta.
//   lax_flux   : the surface kernels with the Lax solver of the wave equation object leave the
//                cell averages of a uniform state alone, and across a jump the rates of the cell
//                averages on the two sides of an interface are those of the Lax-Friedrichs flux.
#include <acutest.h>

#include <float.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_calc_fluid_vars.h>
#include <gkyl_dg_euler.h>
#include <gkyl_dg_euler_priv.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_iso_euler.h>

// Grid, ranges, geometry and updater of a cdim-dimensional p = 1 test.
struct fluid_vars_test {
  struct gkyl_basis basis;
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext;
  struct gkyl_wave_geom *geom;
  struct gkyl_wv_eqn *wv_eqn;
  struct gkyl_dg_calc_fluid_vars *up;
};

// The test of an equation object (which the test then owns).
static struct fluid_vars_test
fluid_vars_test_from_eqn(int cdim, struct gkyl_wv_eqn *wv_eqn)
{
  struct fluid_vars_test ft;
  double lower[3] = {0.0, 0.0, 0.0}, upper[3] = {1.0, 1.0, 1.0};
  int cells[3] = {4, 4, 4}, nghost[3] = {1, 1, 1};
  gkyl_cart_modal_serendip(&ft.basis, cdim, 1);
  gkyl_rect_grid_init(&ft.grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&ft.grid, nghost, &ft.local_ext, &ft.local);
  ft.geom = gkyl_wave_geom_new(&ft.grid, &ft.local_ext, NULL, NULL, false);
  ft.wv_eqn = wv_eqn;
  ft.up = gkyl_dg_calc_fluid_vars_new(ft.wv_eqn, ft.geom, &ft.basis, &ft.local_ext, 0.0, false);
  return ft;
}

static struct fluid_vars_test
fluid_vars_test_new(int cdim, enum gkyl_wv_euler_rp rp_type)
{
  return fluid_vars_test_from_eqn(
    cdim, gkyl_wv_euler_inew(&(struct gkyl_wv_euler_inp){
            .gas_gamma = 5.0 / 3.0,
            .rp_type = rp_type,
            .use_gpu = false,
          })
  );
}

// Isothermal Euler, [rho, rho ux, rho uy, rho uz] with p = vt^2 rho.
static struct fluid_vars_test
fluid_vars_iso_test_new(int cdim, double vt)
{
  return fluid_vars_test_from_eqn(cdim, gkyl_wv_iso_euler_new(vt, false));
}

static void
fluid_vars_test_release(struct fluid_vars_test *ft)
{
  gkyl_dg_calc_fluid_vars_release(ft->up);
  gkyl_wv_eqn_release(ft->wv_eqn);
  gkyl_wave_geom_release(ft->geom);
}

// Coefficient of the constant basis function of a cell average (the basis is orthonormal on
// [-1, 1]^cdim).
static double
avg_coeff(int cdim)
{
  return sqrt((double)(1 << cdim));
}

static void
test_em_source(int cdim)
{
  struct fluid_vars_test ft = fluid_vars_test_new(cdim, WV_EULER_RP_ROE);
  int nb = ft.basis.num_basis;
  double c0 = avg_coeff(cdim);

  struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, 5 * nb, ft.local_ext.volume);
  struct gkyl_array *qmem = gkyl_array_new(GKYL_DOUBLE, 6 * nb, ft.local_ext.volume);
  struct gkyl_array *accel = gkyl_array_new(GKYL_DOUBLE, 3 * nb, ft.local_ext.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, 5 * nb, ft.local_ext.volume);
  struct gkyl_array *rhs_accel = gkyl_array_new(GKYL_DOUBLE, 5 * nb, ft.local_ext.volume);
  struct gkyl_array *omega_p_sq = gkyl_array_new(GKYL_DOUBLE, 1, ft.local_ext.volume);
  struct gkyl_array *cflrate = gkyl_array_new(GKYL_DOUBLE, 1, ft.local_ext.volume);
  gkyl_array_clear(fluid, 0.0);
  gkyl_array_clear(qmem, 0.0);
  gkyl_array_clear(accel, 0.0);
  // Plasma frequency 2; an advective rate of 1 in the first cells along x and 100 in the others.
  double omega_p = 2.0, rate_lo = 1.0, rate_up = 100.0;
  gkyl_array_clear(omega_p_sq, omega_p * omega_p);

  // A uniform flow u in a density that varies linearly in x, and uniform a and q/m B.
  double rho0 = 1.3, rho1 = 0.2, uvec[3] = {0.3, -0.5, 0.7};
  double avec[3] = {0.4, 0.1, -0.6}, bvec[3] = {0.9, -0.2, 0.5};
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ft.local, iter.idx);
    double *fluid_d = gkyl_array_fetch(fluid, loc);
    double *qmem_d = gkyl_array_fetch(qmem, loc);
    double *accel_d = gkyl_array_fetch(accel, loc);
    fluid_d[0] = rho0 * c0;
    fluid_d[1] = rho1;
    fluid_d[4 * nb] = 0.9 * c0; // internal energy; the bulk kinetic energy is added below
    for (int d = 0; d < 3; ++d) {
      fluid_d[4 * nb] += 0.5 * rho0 * uvec[d] * uvec[d] * c0;
      fluid_d[(1 + d) * nb] = uvec[d] * rho0 * c0;
      fluid_d[(1 + d) * nb + 1] = uvec[d] * rho1;
      qmem_d[d * nb] = avec[d] * c0;
      qmem_d[(3 + d) * nb] = bvec[d] * c0;
      accel_d[d * nb] = avec[d] * c0;
    }
  }

  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    double *cflrate_d = gkyl_array_fetch(cflrate, gkyl_range_idx(&ft.local, iter.idx));
    cflrate_d[0] = iter.idx[0] == 1 ? rate_lo : rate_up;
  }

  // rho (a + u x (q/m B)) and rho u . a, with the density expansion of the cell.
  double force[3] = {
    avec[0] + uvec[1] * bvec[2] - uvec[2] * bvec[1],
    avec[1] + uvec[2] * bvec[0] - uvec[0] * bvec[2], avec[2] + uvec[0] * bvec[1] - uvec[1] * bvec[0]
  };
  double work = uvec[0] * avec[0] + uvec[1] * avec[1] + uvec[2] * avec[2];
  // Rate of the source from the cell averages: the larger of omega_p + |q/m B| and
  // 2 |a + u x (q/m B)| / sqrt(2 e / rho), with the energy set below to an internal energy of 0.9.
  double e_int = 0.9;
  double rate_freq = omega_p + sqrt(bvec[0] * bvec[0] + bvec[1] * bvec[1] + bvec[2] * bvec[2]);
  double rate_accel =
    sqrt(2.0 * (force[0] * force[0] + force[1] * force[1] + force[2] * force[2]) * rho0 / e_int);
  double rate_src = fmax(rate_freq, rate_accel);
  TEST_CHECK(rate_src > rate_lo && rate_src < rate_up);

  gkyl_array_clear(rhs, 0.0);
  gkyl_dg_calc_fluid_vars_em_source(ft.up, &ft.local, qmem, omega_p_sq, fluid, cflrate, rhs);
  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ft.local, iter.idx);
    const double *rhs_d = gkyl_array_cfetch(rhs, loc);
    // The CFL rate is the larger of the advective rate and the rate of the source.
    const double *cflrate_d = gkyl_array_cfetch(cflrate, loc);
    TEST_CHECK(gkyl_compare(cflrate_d[0], iter.idx[0] == 1 ? rate_src : rate_up, 1e-13));
    for (int k = 0; k < nb; ++k) {
      double rho_k = k == 0 ? rho0 * c0 : (k == 1 ? rho1 : 0.0);
      TEST_CHECK(gkyl_compare(rhs_d[k], 0.0, 1e-14));
      for (int d = 0; d < 3; ++d) {
        TEST_CHECK(gkyl_compare(rhs_d[(1 + d) * nb + k], rho_k * force[d], 1e-14));
      }
      TEST_CHECK(gkyl_compare(rhs_d[4 * nb + k], rho_k * work, 1e-14));
    }
  }

  // Without a magnetic field: the source of an applied acceleration.
  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    double *qmem_d = gkyl_array_fetch(qmem, gkyl_range_idx(&ft.local, iter.idx));
    for (int d = 0; d < 3; ++d) {
      qmem_d[(3 + d) * nb] = 0.0;
    }
  }
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(rhs_accel, 0.0);
  gkyl_dg_calc_fluid_vars_em_source(ft.up, &ft.local, qmem, omega_p_sq, fluid, cflrate, rhs);
  gkyl_dg_calc_fluid_vars_source(ft.up, &ft.local, accel, fluid, rhs_accel);
  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ft.local, iter.idx);
    const double *rhs_d = gkyl_array_cfetch(rhs, loc);
    const double *rhs_accel_d = gkyl_array_cfetch(rhs_accel, loc);
    for (int k = 0; k < 5 * nb; ++k) {
      TEST_CHECK(gkyl_compare(rhs_d[k], rhs_accel_d[k], 1e-14));
    }
  }

  gkyl_array_release(fluid);
  gkyl_array_release(qmem);
  gkyl_array_release(accel);
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_accel);
  gkyl_array_release(omega_p_sq);
  gkyl_array_release(cflrate);
  fluid_vars_test_release(&ft);
}

// Smallest vertex density and internal energy of a p = 1 cell.
static void
vertex_min(int cdim, int nb, const double *fluid_d, double *rho_min, double *e_min)
{
  double c0 = avg_coeff(cdim);
  *rho_min = DBL_MAX;
  *e_min = DBL_MAX;
  for (int v = 0; v < (1 << cdim); ++v) {
    double q[5];
    for (int m = 0; m < 5; ++m) {
      // Serendipity p = 1: the constant, then the linear functions of each direction.
      q[m] = fluid_d[m * nb] / c0;
      for (int d = 0; d < cdim; ++d) {
        double sgn = (v >> d) & 1 ? 1.0 : -1.0;
        q[m] += sgn * sqrt(3.0) * fluid_d[m * nb + 1 + d] / c0;
      }
    }
    double e = q[4] - 0.5 * (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]) / q[0];
    *rho_min = fmin(*rho_min, q[0]);
    *e_min = fmin(*e_min, e);
  }
}

static void
test_positivity(int cdim)
{
  struct fluid_vars_test ft = fluid_vars_test_new(cdim, WV_EULER_RP_ROE);
  int nb = ft.basis.num_basis;
  double c0 = avg_coeff(cdim), eps = 1.0e-2;

  struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, 5 * nb, ft.local_ext.volume);
  struct gkyl_array *fluid_in = gkyl_array_new(GKYL_DOUBLE, 5 * nb, ft.local_ext.volume);
  gkyl_array_clear(fluid, 0.0);

  // Cell averages rho = 1, E = 1 and M = (0.5, 0, 0) (internal energy 0.875), or M = 0 in the
  // cells of the second kind. Three kinds of cells by the first index: admissible slopes (1), a
  // negative vertex density in a fluid at rest (2), and a negative vertex internal energy with a
  // positive density (3, 4). Only the cross terms of the basis are zero.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    double *fluid_d = gkyl_array_fetch(fluid, gkyl_range_idx(&ft.local, iter.idx));
    fluid_d[0] = c0;
    fluid_d[nb] = iter.idx[0] == 2 ? 0.0 : 0.5 * c0;
    fluid_d[4 * nb] = c0;
    for (int d = 0; d < cdim; ++d) {
      if (iter.idx[0] == 1) {
        fluid_d[1 + d] = 0.1 * c0 / sqrt(3.0);
        fluid_d[4 * nb + 1 + d] = 0.05 * c0 / sqrt(3.0);
      } else if (iter.idx[0] == 2) {
        fluid_d[1 + d] = 1.2 * c0 / sqrt(3.0) / cdim;
      } else {
        fluid_d[4 * nb + 1 + d] = 1.1 * c0 / sqrt(3.0) / cdim;
        fluid_d[nb + 1 + d] = -0.3 * c0 / sqrt(3.0) / cdim;
      }
    }
  }
  gkyl_array_copy(fluid_in, fluid);

  gkyl_dg_calc_fluid_vars_positivity(ft.up, &ft.local, eps, fluid);

  double e_avg = 1.0 - 0.5 * 0.25;
  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ft.local, iter.idx);
    const double *fluid_d = gkyl_array_cfetch(fluid, loc);
    const double *fluid_in_d = gkyl_array_cfetch(fluid_in, loc);
    // The cell averages are untouched.
    for (int m = 0; m < 5; ++m) {
      TEST_CHECK(fluid_d[m * nb] == fluid_in_d[m * nb]);
    }
    double rho_min, e_min, rho_min_in, e_min_in;
    vertex_min(cdim, nb, fluid_d, &rho_min, &e_min);
    vertex_min(cdim, nb, fluid_in_d, &rho_min_in, &e_min_in);
    if (iter.idx[0] == 1) {
      // An admissible cell is not touched.
      for (int k = 0; k < 5 * nb; ++k) {
        TEST_CHECK(fluid_d[k] == fluid_in_d[k]);
      }
    } else if (iter.idx[0] == 2) {
      // The vertex density was negative and is now at the floor (the internal energy of the
      // fluid at rest is its energy, uniform in the cell).
      TEST_CHECK(rho_min_in < 0.0);
      TEST_CHECK(gkyl_compare(rho_min, eps, 1e-12));
      TEST_CHECK(gkyl_compare(e_min, 1.0, 1e-12));
    } else {
      // The vertex internal energy was negative and is now at the floor.
      TEST_CHECK(rho_min_in > eps);
      TEST_CHECK(e_min_in < 0.0);
      TEST_CHECK(gkyl_compare(e_min, eps * e_avg, 1e-12));
      TEST_CHECK(rho_min > eps);
    }
  }

  gkyl_array_release(fluid);
  gkyl_array_release(fluid_in);
  fluid_vars_test_release(&ft);
}

// The surface kernels in 1x with the Lax solver: a uniform state, and a jump between two cells.
static void
test_lax_flux_1x()
{
  struct fluid_vars_test ft = fluid_vars_test_new(1, WV_EULER_RP_LAX);
  int nb = ft.basis.num_basis;
  double c0 = avg_coeff(1), gas_gamma = 5.0 / 3.0;

  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = {0};
  struct gkyl_vlasov_position_map *pos_map =
    gkyl_vlasov_position_map_new(&ft.grid, &ft.local, &ft.local_ext, &ft.basis, inp_pmap, false);
  struct gkyl_array *u = gkyl_array_new(GKYL_DOUBLE, 3 * nb, ft.local_ext.volume);
  struct gkyl_array *p = gkyl_array_new(GKYL_DOUBLE, nb, ft.local_ext.volume);
  struct gkyl_dg_eqn *eqn = gkyl_dg_euler_inew(&(struct gkyl_dg_euler_inp){
    .cbasis = &ft.basis,
    .conf_range = &ft.local,
    .wv_eqn = ft.wv_eqn,
    .geom = ft.geom,
    .pos_map = pos_map,
    .u = u,
    .p = p,
    .use_gpu = false,
  });
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  TEST_CHECK(gkyl_wv_euler_rp_type(ft.wv_eqn) == WV_EULER_RP_LAX);
  TEST_CHECK(euler->surf[0] == euler_surfx_1x_ser_p1);

  // Piecewise-constant states.
  double ql[5] = {1.0, 0.3, -0.2, 0.1, 2.5}, qr[5] = {0.4, -0.1, 0.25, 0.0, 1.1};
  double fluid_l[10] = {0.0}, fluid_r[10] = {0.0};
  for (int m = 0; m < 5; ++m) {
    fluid_l[m * nb] = ql[m] * c0;
    fluid_r[m * nb] = qr[m] * c0;
  }
  double xc[1] = {0.5}, dx[1] = {0.25}, jac[2] = {1.0, 1.0};
  const struct gkyl_wave_cell_geom *cgeom = gkyl_wave_geom_get(ft.geom, (int[]){2});

  // A uniform state: no rate of the cell averages.
  double out[10] = {0.0};
  euler_surfx_1x_ser_p1(xc, dx, jac, ft.wv_eqn, cgeom, cgeom, fluid_l, fluid_l, fluid_l, out);
  for (int m = 0; m < 5; ++m) {
    TEST_CHECK(gkyl_compare(out[m * nb], 0.0, 1e-13));
  }

  // States l | l | r and l | r | r: the interface between the two middle cells.
  double out_l[10] = {0.0}, out_r[10] = {0.0};
  euler_surfx_1x_ser_p1(xc, dx, jac, ft.wv_eqn, cgeom, cgeom, fluid_l, fluid_l, fluid_r, out_l);
  euler_surfx_1x_ser_p1(xc, dx, jac, ft.wv_eqn, cgeom, cgeom, fluid_l, fluid_r, fluid_r, out_r);

  // Lax-Friedrichs flux (F(ql) + F(qr))/2 - max(|u| + c_s) (qr - ql)/2.
  double flux_l[5], flux_r[5], amax = 0.0;
  const double *qs[2] = {ql, qr};
  double *fs[2] = {flux_l, flux_r};
  for (int s = 0; s < 2; ++s) {
    const double *q = qs[s];
    double ux = q[1] / q[0];
    double pr = (gas_gamma - 1.0) * (q[4] - 0.5 * (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]) / q[0]);
    fs[s][0] = q[1];
    fs[s][1] = q[1] * ux + pr;
    fs[s][2] = q[2] * ux;
    fs[s][3] = q[3] * ux;
    fs[s][4] = (q[4] + pr) * ux;
    amax = fmax(amax, fabs(ux) + sqrt(gas_gamma * pr / q[0]));
  }
  for (int m = 0; m < 5; ++m) {
    double fhat = 0.5 * (flux_l[m] + flux_r[m]) - 0.5 * amax * (qr[m] - ql[m]);
    // Rate of the cell averages: -(fhat - F(ql))/dx on the left, (fhat - F(qr))/dx on the right.
    TEST_CHECK(gkyl_compare(out_l[m * nb] / c0, -(fhat - flux_l[m]) / dx[0], 1e-12));
    TEST_CHECK(gkyl_compare(out_r[m * nb] / c0, (fhat - flux_r[m]) / dx[0], 1e-12));
  }

  gkyl_dg_eqn_release(eqn);
  gkyl_array_release(u);
  gkyl_array_release(p);
  gkyl_vlasov_position_map_release(pos_map);
  fluid_vars_test_release(&ft);
}

void
test_em_source_1x()
{
  test_em_source(1);
}

void
test_em_source_2x()
{
  test_em_source(2);
}

void
test_em_source_3x()
{
  test_em_source(3);
}

void
test_positivity_1x()
{
  test_positivity(1);
}

void
test_positivity_2x()
{
  test_positivity(2);
}

void
test_positivity_3x()
{
  test_positivity(3);
}

static void
test_limiter(int cdim)
{
  int nb = 1 << cdim;
  double c0 = avg_coeff(cdim);
  int idx_cell[3] = {2, 2, 2};

  for (int dir = 0; dir < cdim; ++dir) {
    struct fluid_vars_test ft = fluid_vars_test_new(cdim, WV_EULER_RP_ROE);
    struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, 5 * nb, ft.local_ext.volume);
    gkyl_array_clear(fluid, 0.0);

    // A uniform gas at rest, rho = 1 and E = 1.
    for (long linc = 0; linc < ft.local_ext.volume; ++linc) {
      double *fluid_d = gkyl_array_fetch(fluid, linc);
      fluid_d[0] = c0;
      fluid_d[4 * nb] = c0;
    }

    // In one cell: a slope along dir of the momentum along dir, which its neighbours do not
    // support, and the same higher coefficients in all five components.
    double *fluid_d = gkyl_array_fetch(fluid, gkyl_range_idx(&ft.local_ext, idx_cell));
    fluid_d[(1 + dir) * nb + 1 + dir] = 0.05;
    for (int c = 0; c < 5; ++c) {
      for (int k = cdim + 1; k < nb; ++k) {
        fluid_d[c * nb + k] = 0.1;
      }
    }

    gkyl_dg_calc_fluid_vars_limiter(ft.up, &ft.local, fluid);

    TEST_CHECK(fluid_d[(1 + dir) * nb + 1 + dir] == 0.0);
    for (int c = 0; c < 5; ++c) {
      bool is_acoustic = c == 0 || c == 1 + dir || c == 4;
      for (int k = cdim + 1; k < nb; ++k) {
        TEST_CHECK(fluid_d[c * nb + k] == (is_acoustic ? 0.0 : 0.1));
        TEST_MSG("sweep %d, component %d, coefficient %d: %g", dir, c, k, fluid_d[c * nb + k]);
      }
    }
    // The cell averages are untouched.
    TEST_CHECK(fluid_d[0] == c0 && fluid_d[4 * nb] == c0);

    gkyl_array_release(fluid);
    fluid_vars_test_release(&ft);
  }
}

static void
test_limiter_2x()
{
  test_limiter(2);
}

static void
test_limiter_3x()
{
  test_limiter(3);
}

// Isothermal Euler: the pressure is vt^2 rho, the sources are those of Euler on the density and
// the momentum, the rate of the source is omega_p + |q/m B|, and nothing is written past the
// four blocks of the state (the arrays have a fifth block holding a sentinel).
static void
test_iso_sources(int cdim)
{
  double vt = 0.7, sentinel = -3.25;
  struct fluid_vars_test ft = fluid_vars_test_new(cdim, WV_EULER_RP_ROE);
  struct fluid_vars_test ft_iso = fluid_vars_iso_test_new(cdim, vt);
  int nb = ft.basis.num_basis;
  double c0 = avg_coeff(cdim);
  long vol = ft.local_ext.volume;

  struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, 5 * nb, vol);
  struct gkyl_array *qmem = gkyl_array_new(GKYL_DOUBLE, 6 * nb, vol);
  struct gkyl_array *accel = gkyl_array_new(GKYL_DOUBLE, 3 * nb, vol);
  struct gkyl_array *vel = gkyl_array_new(GKYL_DOUBLE, 3 * nb, vol);
  struct gkyl_array *pressure = gkyl_array_new(GKYL_DOUBLE, nb, vol);
  struct gkyl_array *pressure_surf = gkyl_array_new(GKYL_DOUBLE, 2 * cdim * nb / 2, vol);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, 5 * nb, vol);
  struct gkyl_array *rhs_iso = gkyl_array_new(GKYL_DOUBLE, 5 * nb, vol);
  struct gkyl_array *omega_p_sq = gkyl_array_new(GKYL_DOUBLE, 1, vol);
  struct gkyl_array *cflrate = gkyl_array_new(GKYL_DOUBLE, 1, vol);
  struct gkyl_array *cflrate_iso = gkyl_array_new(GKYL_DOUBLE, 1, vol);
  gkyl_array_clear(fluid, 0.0);
  gkyl_array_clear(qmem, 0.0);
  gkyl_array_clear(accel, 0.0);
  gkyl_array_clear(vel, 0.0);
  double omega_p = 2.0, rate_lo = 1.0, rate_up = 100.0;
  gkyl_array_clear(omega_p_sq, omega_p * omega_p);

  // A uniform flow u in a density that varies linearly in x, and uniform a and q/m B.
  double rho0 = 1.3, rho1 = 0.2, uvec[3] = {0.3, -0.5, 0.7};
  double avec[3] = {0.4, 0.1, -0.6}, bvec[3] = {0.9, -0.2, 0.5};
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ft.local, iter.idx);
    double *fluid_d = gkyl_array_fetch(fluid, loc);
    double *qmem_d = gkyl_array_fetch(qmem, loc);
    double *accel_d = gkyl_array_fetch(accel, loc);
    double *cflrate_d = gkyl_array_fetch(cflrate, loc);
    fluid_d[0] = rho0 * c0;
    fluid_d[1] = rho1;
    fluid_d[4 * nb] = 5.0 * c0; // the energy of the Euler reference
    for (int d = 0; d < 3; ++d) {
      fluid_d[(1 + d) * nb] = uvec[d] * rho0 * c0;
      fluid_d[(1 + d) * nb + 1] = uvec[d] * rho1;
      qmem_d[d * nb] = avec[d] * c0;
      qmem_d[(3 + d) * nb] = bvec[d] * c0;
      accel_d[d * nb] = avec[d] * c0;
    }
    cflrate_d[0] = iter.idx[0] == 1 ? rate_lo : rate_up;
  }
  gkyl_array_copy(cflrate_iso, cflrate);

  gkyl_dg_calc_fluid_vars_pressure(ft_iso.up, &ft.local, fluid, vel, pressure, pressure_surf);

  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(rhs_iso, 0.0);
  gkyl_array_shiftc(rhs_iso, sentinel, 4 * nb);
  gkyl_dg_calc_fluid_vars_em_source(ft.up, &ft.local, qmem, omega_p_sq, fluid, cflrate, rhs);
  gkyl_dg_calc_fluid_vars_em_source(
    ft_iso.up, &ft.local, qmem, omega_p_sq, fluid, cflrate_iso, rhs_iso
  );
  double rate_freq = omega_p + sqrt(bvec[0] * bvec[0] + bvec[1] * bvec[1] + bvec[2] * bvec[2]);
  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ft.local, iter.idx);
    const double *fluid_d = gkyl_array_cfetch(fluid, loc);
    const double *pressure_d = gkyl_array_cfetch(pressure, loc);
    const double *rhs_d = gkyl_array_cfetch(rhs, loc);
    const double *rhs_iso_d = gkyl_array_cfetch(rhs_iso, loc);
    const double *cflrate_iso_d = gkyl_array_cfetch(cflrate_iso, loc);
    for (int k = 0; k < nb; ++k) {
      TEST_CHECK(gkyl_compare(pressure_d[k], vt * vt * fluid_d[k], 1e-15));
    }
    for (int k = 0; k < 4 * nb; ++k) {
      TEST_CHECK(rhs_iso_d[k] == rhs_d[k]);
    }
    TEST_CHECK(rhs_iso_d[4 * nb] == sentinel);
    for (int k = 1; k < nb; ++k) {
      TEST_CHECK(rhs_iso_d[4 * nb + k] == 0.0);
    }
    TEST_CHECK(gkyl_compare(cflrate_iso_d[0], iter.idx[0] == 1 ? rate_freq : rate_up, 1e-13));
  }

  // The source of an applied acceleration.
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(rhs_iso, 0.0);
  gkyl_array_shiftc(rhs_iso, sentinel, 4 * nb);
  gkyl_dg_calc_fluid_vars_source(ft.up, &ft.local, accel, fluid, rhs);
  gkyl_dg_calc_fluid_vars_source(ft_iso.up, &ft.local, accel, fluid, rhs_iso);
  double max_rhs = 0.0;
  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ft.local, iter.idx);
    const double *rhs_d = gkyl_array_cfetch(rhs, loc);
    const double *rhs_iso_d = gkyl_array_cfetch(rhs_iso, loc);
    for (int k = 0; k < 4 * nb; ++k) {
      TEST_CHECK(rhs_iso_d[k] == rhs_d[k]);
      max_rhs = fmax(max_rhs, fabs(rhs_iso_d[k]));
    }
    TEST_CHECK(rhs_iso_d[4 * nb] == sentinel);
  }
  TEST_CHECK(max_rhs > 0.0);

  gkyl_array_release(fluid);
  gkyl_array_release(qmem);
  gkyl_array_release(accel);
  gkyl_array_release(vel);
  gkyl_array_release(pressure);
  gkyl_array_release(pressure_surf);
  gkyl_array_release(rhs);
  gkyl_array_release(rhs_iso);
  gkyl_array_release(omega_p_sq);
  gkyl_array_release(cflrate);
  gkyl_array_release(cflrate_iso);
  fluid_vars_test_release(&ft);
  fluid_vars_test_release(&ft_iso);
}

// Isothermal Euler: the positivity limiter scales the density slopes only, to the floor.
static void
test_iso_positivity(int cdim)
{
  struct fluid_vars_test ft = fluid_vars_iso_test_new(cdim, 0.7);
  int nb = ft.basis.num_basis;
  double c0 = avg_coeff(cdim), eps = 1.0e-2;

  struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, 4 * nb, ft.local_ext.volume);
  struct gkyl_array *fluid_in = gkyl_array_new(GKYL_DOUBLE, 4 * nb, ft.local_ext.volume);
  gkyl_array_clear(fluid, 0.0);

  // Cell averages rho = 1 and M = (0.5, 0, 0) with a momentum slope; admissible density slopes
  // in the cells of the first kind (by the first index), a negative vertex density in the others.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    double *fluid_d = gkyl_array_fetch(fluid, gkyl_range_idx(&ft.local, iter.idx));
    fluid_d[0] = c0;
    fluid_d[nb] = 0.5 * c0;
    for (int d = 0; d < cdim; ++d) {
      fluid_d[1 + d] = (iter.idx[0] == 1 ? 0.1 : 1.2 / cdim) * c0 / sqrt(3.0);
      fluid_d[nb + 1 + d] = -0.3 * c0 / sqrt(3.0);
    }
  }
  gkyl_array_copy(fluid_in, fluid);

  gkyl_dg_calc_fluid_vars_positivity(ft.up, &ft.local, eps, fluid);

  gkyl_range_iter_init(&iter, &ft.local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&ft.local, iter.idx);
    const double *fluid_d = gkyl_array_cfetch(fluid, loc);
    const double *fluid_in_d = gkyl_array_cfetch(fluid_in, loc);
    // The cell average of the density and the momentum are untouched.
    TEST_CHECK(fluid_d[0] == fluid_in_d[0]);
    for (int k = nb; k < 4 * nb; ++k) {
      TEST_CHECK(fluid_d[k] == fluid_in_d[k]);
    }
    double rho_min = DBL_MAX, rho_min_in = DBL_MAX;
    for (int v = 0; v < (1 << cdim); ++v) {
      double rho_v = fluid_d[0] / c0, rho_v_in = fluid_in_d[0] / c0;
      for (int d = 0; d < cdim; ++d) {
        double sgn = (v >> d) & 1 ? 1.0 : -1.0;
        rho_v += sgn * sqrt(3.0) * fluid_d[1 + d] / c0;
        rho_v_in += sgn * sqrt(3.0) * fluid_in_d[1 + d] / c0;
      }
      rho_min = fmin(rho_min, rho_v);
      rho_min_in = fmin(rho_min_in, rho_v_in);
    }
    if (iter.idx[0] == 1) {
      for (int k = 0; k < nb; ++k) {
        TEST_CHECK(fluid_d[k] == fluid_in_d[k]);
      }
    } else {
      TEST_CHECK(rho_min_in < 0.0);
      TEST_CHECK(gkyl_compare(rho_min, eps, 1e-12));
    }
  }

  gkyl_array_release(fluid);
  gkyl_array_release(fluid_in);
  fluid_vars_test_release(&ft);
}

// Isothermal Euler: the slope limiter removes a momentum slope its neighbours do not support
// and drops the higher coefficients of the components of the limited (sound) waves; a density
// ramp its neighbours support is kept.
static void
test_iso_limiter(int cdim)
{
  int nb = 1 << cdim;
  double c0 = avg_coeff(cdim);
  int idx_cell[3] = {2, 2, 2};

  for (int dir = 0; dir < cdim; ++dir) {
    struct fluid_vars_test ft = fluid_vars_iso_test_new(cdim, 0.7);
    struct gkyl_array *fluid = gkyl_array_new(GKYL_DOUBLE, 4 * nb, ft.local_ext.volume);
    gkyl_array_clear(fluid, 0.0);

    // A uniform gas at rest, rho = 1.
    for (long linc = 0; linc < ft.local_ext.volume; ++linc) {
      double *fluid_d = gkyl_array_fetch(fluid, linc);
      fluid_d[0] = c0;
    }

    // In one cell: a slope along dir of the momentum along dir, which its neighbours do not
    // support, and the same higher coefficients in all four components.
    double *fluid_d = gkyl_array_fetch(fluid, gkyl_range_idx(&ft.local_ext, idx_cell));
    fluid_d[(1 + dir) * nb + 1 + dir] = 0.05;
    for (int c = 0; c < 4; ++c) {
      for (int k = cdim + 1; k < nb; ++k) {
        fluid_d[c * nb + k] = 0.1;
      }
    }

    gkyl_dg_calc_fluid_vars_limiter(ft.up, &ft.local, fluid);

    TEST_CHECK(fluid_d[(1 + dir) * nb + 1 + dir] == 0.0);
    for (int c = 0; c < 4; ++c) {
      bool is_acoustic = c == 0 || c == 1 + dir;
      for (int k = cdim + 1; k < nb; ++k) {
        TEST_CHECK(fluid_d[c * nb + k] == (is_acoustic ? 0.0 : 0.1));
        TEST_MSG("sweep %d, component %d, coefficient %d: %g", dir, c, k, fluid_d[c * nb + k]);
      }
    }
    // The cell averages are untouched.
    TEST_CHECK(fluid_d[0] == c0);

    // A density ramp along dir across all cells (slope = the difference of the neighbouring
    // averages over sqrt(3), the limiter's own estimate) is kept.
    double ramp = 0.02;
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &ft.local_ext);
    while (gkyl_range_iter_next(&iter)) {
      double *ramp_d = gkyl_array_fetch(fluid, gkyl_range_idx(&ft.local_ext, iter.idx));
      for (int k = 0; k < 4 * nb; ++k) {
        ramp_d[k] = 0.0;
      }
      ramp_d[0] = (1.0 + ramp * iter.idx[dir]) * c0;
      ramp_d[1 + dir] = ramp * c0 / (2.0 * sqrt(3.0));
    }
    gkyl_dg_calc_fluid_vars_limiter(ft.up, &ft.local, fluid);
    TEST_CHECK(gkyl_compare(fluid_d[1 + dir], ramp * c0 / (2.0 * sqrt(3.0)), 1e-13));
    TEST_MSG("sweep %d: slope %g of %g", dir, fluid_d[1 + dir], ramp * c0 / (2.0 * sqrt(3.0)));

    gkyl_array_release(fluid);
    fluid_vars_test_release(&ft);
  }
}

static void
test_iso_sources_1x()
{
  test_iso_sources(1);
}

static void
test_iso_sources_2x()
{
  test_iso_sources(2);
}

static void
test_iso_sources_3x()
{
  test_iso_sources(3);
}

static void
test_iso_positivity_1x()
{
  test_iso_positivity(1);
}

static void
test_iso_positivity_2x()
{
  test_iso_positivity(2);
}

static void
test_iso_positivity_3x()
{
  test_iso_positivity(3);
}

static void
test_iso_limiter_1x()
{
  test_iso_limiter(1);
}

static void
test_iso_limiter_2x()
{
  test_iso_limiter(2);
}

static void
test_iso_limiter_3x()
{
  test_iso_limiter(3);
}

TEST_LIST = {
  {"em_source_1x", test_em_source_1x},
  {"em_source_2x", test_em_source_2x},
  {"em_source_3x", test_em_source_3x},
  {"positivity_1x", test_positivity_1x},
  {"positivity_2x", test_positivity_2x},
  {"positivity_3x", test_positivity_3x},
  {"lax_flux_1x", test_lax_flux_1x},
  {"limiter_2x", test_limiter_2x},
  {"limiter_3x", test_limiter_3x},
  {"iso_sources_1x", test_iso_sources_1x},
  {"iso_sources_2x", test_iso_sources_2x},
  {"iso_sources_3x", test_iso_sources_3x},
  {"iso_positivity_1x", test_iso_positivity_1x},
  {"iso_positivity_2x", test_iso_positivity_2x},
  {"iso_positivity_3x", test_iso_positivity_3x},
  {"iso_limiter_1x", test_iso_limiter_1x},
  {"iso_limiter_2x", test_iso_limiter_2x},
  {"iso_limiter_3x", test_iso_limiter_3x},
  {NULL, NULL}
};
