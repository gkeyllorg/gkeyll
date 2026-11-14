// The PKPM LBO in the configuration-nodal, velocity-modal form (velocity-only kernels at each node
// through hyper_dg). On each of the 8 blocks h = [F, G, M0_perp (3), N_perp (3)] the operator is
//   d_w[nu (w - u) h + nu vtsq d_w h] + S,
// with zero flux at the ends of the velocity grid and the couplings S_G = 2 nu (vtsq F - G),
// S_M0 = -nu M0_perp and S_N = 4 nu vtsq M0_perp - 3 nu N_perp. Checked here:
// 1. The moments of the right-hand side at every node, which are exact in the weak form:
//      int rhs dw     = int S dw,
//      int w rhs dw   = -nu (M1 - u M0) - nu vtsq [h] + int w S dw,
//      int w^2 rhs dw = -2 nu (M2 - u M1) + 2 nu vtsq (M0 - [w h]) + int w^2 S dw,
//    with M0, M1, M2 the moments of h and [g] = g(w_max) - g(w_min), for random states; and the CFL
//    rate, the largest rate of the nodes of a cell.
// 2. In 2x and 3x, the right-hand side at a node is that of the 1x operator for the same node
//    values: the nodes of a cell do not couple.
// 3. In 1x, convergence to the operator applied to smooth distributions as the velocity grid is
//    refined.
#include <math.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_updater_lbo_pkpm_nodal.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#define NB 8 // blocks of the distribution functions
#define NVM 3 // velocity modes of a node

static double
rnd(void)
{
  return 2.0 * rand() / (double)RAND_MAX - 1.0;
}

// Orthonormal Legendre modes of a velocity cell at the logical coordinate eta.
static void
vel_modes(double eta, double *modes)
{
  modes[0] = sqrt(0.5);
  modes[1] = sqrt(1.5) * eta;
  modes[2] = sqrt(0.625) * (3.0 * eta * eta - 1.0);
}

// Coupling S of block c at a point from the eight block values there.
static double
coupling(int c, double nu, double vtsq, const double *h)
{
  if (c == 0) {
    return 0.0;
  }
  if (c == 1) {
    return 2.0 * nu * (vtsq * h[0] - h[1]);
  }
  if (c < 5) {
    return -nu * h[c];
  }
  return 4.0 * nu * vtsq * h[c - 3] - 3.0 * nu * h[c];
}

// One application of the operator on a grid of num_cells^cdim configuration cells and num_vel
// velocity cells: the node values of nu and [u, vtsq], the nodal form of the distributions and the
// outputs.
struct lbo_case {
  int cdim, num_nodes, num_basis, num_vel;
  struct gkyl_basis cbasis, pbasis;
  struct gkyl_rect_grid phase_grid;
  struct gkyl_range crange, crange_ext, prange, prange_ext;
  struct gkyl_array *nu, *prim, *fnodal, *rhs, *cfl;
};

static void
case_init(struct lbo_case *s, int cdim, int num_cells, int num_vel, double w_lo, double w_up)
{
  s->cdim = cdim;
  s->num_nodes = 1 << cdim;
  s->num_basis = NVM * s->num_nodes;
  s->num_vel = num_vel;
  gkyl_cart_modal_serendip(&s->cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&s->pbasis, cdim, 1);

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    lower[d] = 0.0;
    upper[d] = 1.0 + d;
    cells[d] = num_cells;
    ghost[d] = 1;
  }
  lower[cdim] = w_lo;
  upper[cdim] = w_up;
  cells[cdim] = num_vel;
  ghost[cdim] = 0;
  struct gkyl_rect_grid conf_grid;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&s->phase_grid, cdim + 1, lower, upper, cells);
  gkyl_create_grid_ranges(&conf_grid, ghost, &s->crange_ext, &s->crange);
  gkyl_create_grid_ranges(&s->phase_grid, ghost, &s->prange_ext, &s->prange);

  s->nu = gkyl_array_new(GKYL_DOUBLE, s->num_nodes, s->crange_ext.volume);
  s->prim = gkyl_array_new(GKYL_DOUBLE, 2 * s->num_nodes, s->crange_ext.volume);
  s->fnodal = gkyl_array_new(GKYL_DOUBLE, NB * s->num_basis, s->prange_ext.volume);
  s->rhs = gkyl_array_new(GKYL_DOUBLE, NB * s->num_basis, s->prange_ext.volume);
  s->cfl = gkyl_array_new(GKYL_DOUBLE, 1, s->prange_ext.volume);
  gkyl_array_clear(s->nu, 0.0);
  gkyl_array_clear(s->prim, 0.0);
  gkyl_array_clear(s->fnodal, 0.0);
}

static void
case_advance(struct lbo_case *s)
{
  struct gkyl_dg_lbo_pkpm_nodal_auxfields aux = {.nu = s->nu, .prim_moms = s->prim};
  struct gkyl_dg_updater_lbo_pkpm_nodal *up = gkyl_dg_updater_lbo_pkpm_nodal_new(
    &s->phase_grid, &s->cbasis, &s->pbasis, &s->crange, &aux, false
  );
  gkyl_array_clear(s->rhs, 0.0);
  gkyl_array_clear(s->cfl, 0.0);
  gkyl_dg_updater_lbo_pkpm_nodal_advance(up, &s->prange, s->fnodal, s->cfl, s->rhs);
  gkyl_dg_updater_lbo_pkpm_nodal_release(up);
}

static void
case_release(struct lbo_case *s)
{
  gkyl_array_release(s->nu);
  gkyl_array_release(s->prim);
  gkyl_array_release(s->fnodal);
  gkyl_array_release(s->rhs);
  gkyl_array_release(s->cfl);
}

// A frequency near 2, a flow inside the velocity grid (so the drag changes direction inside it)
// and a thermal speed squared near 1.5, different at every node; random distributions.
static void
case_randomize(struct lbo_case *s)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s->crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&s->crange, iter.idx);
    double *nu_d = gkyl_array_fetch(s->nu, loc), *prim_d = gkyl_array_fetch(s->prim, loc);
    for (int q = 0; q < s->num_nodes; ++q) {
      nu_d[q] = 2.0 + 0.3 * rnd();
      prim_d[q] = 0.3 + 0.5 * rnd();
      prim_d[s->num_nodes + q] = 1.5 + 0.2 * rnd();
    }
  }
  gkyl_range_iter_init(&iter, &s->prange);
  while (gkyl_range_iter_next(&iter)) {
    double *f_d = gkyl_array_fetch(s->fnodal, gkyl_range_idx(&s->prange, iter.idx));
    for (int i = 0; i < NB * s->num_basis; ++i) {
      f_d[i] = rnd();
    }
  }
}

// Moments [int g dw, int w g dw, int w^2 g dw] over the velocity grid of the eight blocks of an
// array in the nodal form at node q of a configuration cell, and the block values at the two ends
// of the velocity grid.
static void
node_moments(
  const struct lbo_case *s, const struct gkyl_array *arr, const int *conf_idx, int q,
  double moms[NB][3], double ends[NB][2]
)
{
  int cdim = s->cdim, idx[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    idx[d] = conf_idx[d];
  }
  for (int c = 0; c < NB; ++c) {
    moms[c][0] = moms[c][1] = moms[c][2] = 0.0;
  }
  double dw = s->phase_grid.dx[cdim], xc[GKYL_MAX_DIM], modes[NVM];
  for (int k = s->prange.lower[cdim]; k <= s->prange.upper[cdim]; ++k) {
    idx[cdim] = k;
    gkyl_rect_grid_cell_center(&s->phase_grid, idx, xc);
    const double *arr_d = gkyl_array_cfetch(arr, gkyl_range_idx(&s->prange, idx));
    for (int c = 0; c < NB; ++c) {
      const double *coef = &arr_d[c * s->num_basis + q * NVM];
      for (int g = 0; g < 3; ++g) {
        double eta = gkyl_gauss_ordinates_3[g], w = xc[cdim] + 0.5 * dw * eta;
        vel_modes(eta, modes);
        double val = coef[0] * modes[0] + coef[1] * modes[1] + coef[2] * modes[2];
        double wt = 0.5 * dw * gkyl_gauss_weights_3[g];
        moms[c][0] += wt * val;
        moms[c][1] += wt * w * val;
        moms[c][2] += wt * w * w * val;
      }
      if (k == s->prange.lower[cdim]) {
        vel_modes(-1.0, modes);
        ends[c][0] = coef[0] * modes[0] + coef[1] * modes[1] + coef[2] * modes[2];
      }
      if (k == s->prange.upper[cdim]) {
        vel_modes(1.0, modes);
        ends[c][1] = coef[0] * modes[0] + coef[1] * modes[1] + coef[2] * modes[2];
      }
    }
  }
}

static void
test_moments(int cdim)
{
  struct lbo_case s;
  double w_lo = -2.0, w_up = 2.5;
  case_init(&s, cdim, 2, 6, w_lo, w_up);
  case_randomize(&s);
  case_advance(&s);
  int nn = s.num_nodes;

  double err_mom = 0.0, max_mom = 0.0;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.crange);
  while (gkyl_range_iter_next(&iter)) {
    long cloc = gkyl_range_idx(&s.crange, iter.idx);
    const double *nu_d = gkyl_array_cfetch(s.nu, cloc), *prim_d = gkyl_array_cfetch(s.prim, cloc);
    for (int q = 0; q < nn; ++q) {
      double nu = nu_d[q], u = prim_d[q], vtsq = prim_d[nn + q];
      double mf[NB][3], mr[NB][3], ef[NB][2], er[NB][2];
      node_moments(&s, s.fnodal, iter.idx, q, mf, ef);
      node_moments(&s, s.rhs, iter.idx, q, mr, er);
      for (int c = 0; c < NB; ++c) {
        // the couplings are linear in the blocks, so their moments are the couplings of the moments
        double src[3];
        for (int n = 0; n < 3; ++n) {
          double h[NB];
          for (int b = 0; b < NB; ++b) {
            h[b] = mf[b][n];
          }
          src[n] = coupling(c, nu, vtsq, h);
        }
        double jump = ef[c][1] - ef[c][0], jump_w = w_up * ef[c][1] - w_lo * ef[c][0];
        double expected[3] = {
          src[0], -nu * (mf[c][1] - u * mf[c][0]) - nu * vtsq * jump + src[1],
          -2.0 * nu * (mf[c][2] - u * mf[c][1]) + 2.0 * nu * vtsq * (mf[c][0] - jump_w) + src[2]
        };
        for (int n = 0; n < 3; ++n) {
          err_mom = fmax(err_mom, fabs(mr[c][n] - expected[n]));
          max_mom = fmax(max_mom, fabs(expected[n]));
        }
      }
    }
  }
  TEST_CHECK(max_mom > 1.0);
  TEST_CHECK(err_mom < 1.0e-13 * max_mom);
  TEST_MSG("cdim %d: moments of the operator off by %g (largest %g)", cdim, err_mom, max_mom);

  // largest node rate: 2.5 (2/dv) nu |u - w_c| + 3 nu + 9 (2/dv)^2 nu vtsq
  double err_cfl = 0.0, max_cfl = 0.0;
  gkyl_range_iter_init(&iter, &s.prange);
  while (gkyl_range_iter_next(&iter)) {
    long cloc = gkyl_range_idx(&s.crange, iter.idx);
    const double *nu_d = gkyl_array_cfetch(s.nu, cloc), *prim_d = gkyl_array_cfetch(s.prim, cloc);
    double xc[GKYL_MAX_DIM], rdv = 2.0 / s.phase_grid.dx[cdim], rate_max = 0.0;
    gkyl_rect_grid_cell_center(&s.phase_grid, iter.idx, xc);
    for (int q = 0; q < nn; ++q) {
      double rate = 2.5 * rdv * nu_d[q] * fabs(prim_d[q] - xc[cdim]) + 3.0 * nu_d[q] +
                    9.0 * rdv * rdv * nu_d[q] * prim_d[nn + q];
      rate_max = fmax(rate_max, rate);
    }
    const double *cfl_d = gkyl_array_cfetch(s.cfl, gkyl_range_idx(&s.prange, iter.idx));
    err_cfl = fmax(err_cfl, fabs(cfl_d[0] - rate_max));
    max_cfl = fmax(max_cfl, rate_max);
  }
  TEST_CHECK(err_cfl < 1.0e-13 * max_cfl);
  TEST_MSG("CFL rate differs from the largest node rate by %g (largest %g)", err_cfl, max_cfl);

  case_release(&s);
}

// The nodes of the cells of a cdim-dimensional case, numbered cell by cell, laid out as the nodes
// of a 1x case with the same velocity grid: node g is node g % 2 of the 1x cell g / 2.
static void
test_nodes_uncoupled(int cdim)
{
  struct lbo_case s, s1;
  double w_lo = -2.0, w_up = 2.5;
  int num_vel = 6;
  case_init(&s, cdim, 2, num_vel, w_lo, w_up);
  case_randomize(&s);
  int nn = s.num_nodes, num_cells_1x = (int)s.crange.volume * nn / 2;
  case_init(&s1, 1, num_cells_1x, num_vel, w_lo, w_up);

  struct gkyl_range_iter iter;
  long cell = 0;
  gkyl_range_iter_init(&iter, &s.crange);
  while (gkyl_range_iter_next(&iter)) {
    long cloc = gkyl_range_idx(&s.crange, iter.idx);
    const double *nu_d = gkyl_array_cfetch(s.nu, cloc), *prim_d = gkyl_array_cfetch(s.prim, cloc);
    for (int q = 0; q < nn; ++q) {
      long node = cell * nn + q;
      int idx1[2] = {s1.crange.lower[0] + (int)(node / 2), 0}, q1 = (int)(node % 2);
      long cloc1 = gkyl_range_idx(&s1.crange, idx1);
      ((double *)gkyl_array_fetch(s1.nu, cloc1))[q1] = nu_d[q];
      double *prim1_d = gkyl_array_fetch(s1.prim, cloc1);
      prim1_d[q1] = prim_d[q];
      prim1_d[2 + q1] = prim_d[nn + q];
      int idx[GKYL_MAX_DIM];
      for (int d = 0; d < cdim; ++d) {
        idx[d] = iter.idx[d];
      }
      for (int k = 0; k < num_vel; ++k) {
        idx[cdim] = s.prange.lower[cdim] + k;
        idx1[1] = s1.prange.lower[1] + k;
        const double *f_d = gkyl_array_cfetch(s.fnodal, gkyl_range_idx(&s.prange, idx));
        double *f1_d = gkyl_array_fetch(s1.fnodal, gkyl_range_idx(&s1.prange, idx1));
        for (int c = 0; c < NB; ++c) {
          for (int j = 0; j < NVM; ++j) {
            f1_d[c * s1.num_basis + q1 * NVM + j] = f_d[c * s.num_basis + q * NVM + j];
          }
        }
      }
    }
    cell += 1;
  }
  case_advance(&s);
  case_advance(&s1);

  double err_rhs = 0.0, max_rhs = 0.0, err_cfl = 0.0;
  cell = 0;
  gkyl_range_iter_init(&iter, &s.crange);
  while (gkyl_range_iter_next(&iter)) {
    int idx[GKYL_MAX_DIM];
    for (int d = 0; d < cdim; ++d) {
      idx[d] = iter.idx[d];
    }
    for (int k = 0; k < num_vel; ++k) {
      idx[cdim] = s.prange.lower[cdim] + k;
      long ploc = gkyl_range_idx(&s.prange, idx);
      const double *rhs_d = gkyl_array_cfetch(s.rhs, ploc);
      double cfl_1x = 0.0;
      for (int q = 0; q < nn; ++q) {
        long node = cell * nn + q;
        int idx1[2] = {s1.crange.lower[0] + (int)(node / 2), s1.prange.lower[1] + k};
        int q1 = (int)(node % 2);
        long ploc1 = gkyl_range_idx(&s1.prange, idx1);
        const double *rhs1_d = gkyl_array_cfetch(s1.rhs, ploc1);
        for (int c = 0; c < NB; ++c) {
          for (int j = 0; j < NVM; ++j) {
            double val = rhs_d[c * s.num_basis + q * NVM + j];
            err_rhs = fmax(err_rhs, fabs(val - rhs1_d[c * s1.num_basis + q1 * NVM + j]));
            max_rhs = fmax(max_rhs, fabs(val));
          }
        }
        cfl_1x = fmax(cfl_1x, ((const double *)gkyl_array_cfetch(s1.cfl, ploc1))[0]);
      }
      err_cfl = fmax(err_cfl, fabs(((const double *)gkyl_array_cfetch(s.cfl, ploc))[0] - cfl_1x));
    }
    cell += 1;
  }
  TEST_CHECK(max_rhs > 1.0);
  TEST_CHECK(err_rhs < 1.0e-13 * max_rhs);
  TEST_MSG("cdim %d: operator differs from 1x by %g (largest %g)", cdim, err_rhs, max_rhs);
  TEST_CHECK(err_cfl == 0.0);
  TEST_MSG("cdim %d: CFL rate differs from 1x by %g", cdim, err_cfl);

  case_release(&s);
  case_release(&s1);
}

// Block c of the smooth state: a Gaussian of its own amplitude, center and width, with its first
// two derivatives.
static void
smooth_block(int c, double w, double *h, double *dh, double *ddh)
{
  double amp = 1.0 + 0.2 * c, center = 0.15 * (c - 3), var = 0.7 + 0.1 * c;
  double val = amp * exp(-0.5 * (w - center) * (w - center) / var);
  *h = val;
  *dh = -(w - center) / var * val;
  *ddh = ((w - center) * (w - center) / (var * var) - 1.0 / var) * val;
}

// Largest difference, over the velocity modes of all cells, nodes and blocks, between the operator
// and the projection of d_w[nu (w - u) h + nu vtsq d_w h] + S for the smooth state.
static double
smooth_error(int num_vel, double *max_val)
{
  struct lbo_case s;
  double w_lo = -7.0, w_up = 7.0;
  case_init(&s, 1, 1, num_vel, w_lo, w_up);
  const double nu_n[2] = {2.0, 1.3}, u_n[2] = {0.3, -0.4}, vtsq_n[2] = {1.2, 0.8};
  int cidx[1] = {s.crange.lower[0]};
  long cloc = gkyl_range_idx(&s.crange, cidx);
  double *nu_d = gkyl_array_fetch(s.nu, cloc), *prim_d = gkyl_array_fetch(s.prim, cloc);
  for (int q = 0; q < 2; ++q) {
    nu_d[q] = nu_n[q];
    prim_d[q] = u_n[q];
    prim_d[2 + q] = vtsq_n[q];
  }

  struct gkyl_array *ref = gkyl_array_new(GKYL_DOUBLE, NB * s.num_basis, s.prange_ext.volume);
  gkyl_array_clear(ref, 0.0);
  double dw = s.phase_grid.dx[1], xc[2], modes[NVM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &s.prange);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(&s.prange, iter.idx);
    double *f_d = gkyl_array_fetch(s.fnodal, ploc), *ref_d = gkyl_array_fetch(ref, ploc);
    gkyl_rect_grid_cell_center(&s.phase_grid, iter.idx, xc);
    for (int g = 0; g < 5; ++g) {
      double eta = gkyl_gauss_ordinates_5[g], w = xc[1] + 0.5 * dw * eta;
      vel_modes(eta, modes);
      double h[NB], dh[NB], ddh[NB];
      for (int c = 0; c < NB; ++c) {
        smooth_block(c, w, &h[c], &dh[c], &ddh[c]);
      }
      for (int q = 0; q < 2; ++q) {
        for (int c = 0; c < NB; ++c) {
          double op = nu_n[q] * (h[c] + (w - u_n[q]) * dh[c]) + nu_n[q] * vtsq_n[q] * ddh[c] +
                      coupling(c, nu_n[q], vtsq_n[q], h);
          for (int j = 0; j < NVM; ++j) {
            f_d[c * s.num_basis + q * NVM + j] += gkyl_gauss_weights_5[g] * modes[j] * h[c];
            ref_d[c * s.num_basis + q * NVM + j] += gkyl_gauss_weights_5[g] * modes[j] * op;
          }
        }
      }
    }
  }
  case_advance(&s);

  double err = 0.0;
  *max_val = 0.0;
  gkyl_range_iter_init(&iter, &s.prange);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(&s.prange, iter.idx);
    const double *rhs_d = gkyl_array_cfetch(s.rhs, ploc), *ref_d = gkyl_array_cfetch(ref, ploc);
    for (int i = 0; i < NB * s.num_basis; ++i) {
      err = fmax(err, fabs(rhs_d[i] - ref_d[i]));
      *max_val = fmax(*max_val, fabs(ref_d[i]));
    }
  }
  gkyl_array_release(ref);
  case_release(&s);
  return err;
}

static void
test_smooth(void)
{
  double max_val[3];
  double err[3] = {
    smooth_error(56, &max_val[0]), smooth_error(112, &max_val[1]), smooth_error(224, &max_val[2])
  };
  TEST_CHECK(err[0] / err[1] > 3.0 && err[1] / err[2] > 3.0);
  TEST_MSG(
    "difference from the smooth operator at 56, 112, 224 velocity cells: %g %g %g (largest %g)",
    err[0], err[1], err[2], max_val[2]
  );
  TEST_CHECK(err[2] < 1.0e-3 * max_val[2]);
}

static void
test_moments_1x(void)
{
  test_moments(1);
}

static void
test_moments_2x(void)
{
  test_moments(2);
}

static void
test_moments_3x(void)
{
  test_moments(3);
}

static void
test_nodes_uncoupled_2x(void)
{
  test_nodes_uncoupled(2);
}

static void
test_nodes_uncoupled_3x(void)
{
  test_nodes_uncoupled(3);
}

TEST_LIST = {
  {"moments_1x", test_moments_1x},
  {"moments_2x", test_moments_2x},
  {"moments_3x", test_moments_3x},
  {"nodes_uncoupled_2x", test_nodes_uncoupled_2x},
  {"nodes_uncoupled_3x", test_nodes_uncoupled_3x},
  {"smooth", test_smooth},
  {NULL, NULL}
};
