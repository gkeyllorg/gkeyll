// The primitive moments [u, vtsq] of the PKPM LBO at the configuration-space nodes: at every node
// they satisfy
//   M0 u - cM vtsq = M1,
//   M1 u + (3 M0 - cE) vtsq = M2 + 2 M0_G,
// with cM = F(w_max) - F(w_min) and cE = w_max F(w_max) - w_min F(w_min). The coefficients are
// evaluated here from the modal expansions at the node, whose position in the cell is read off the
// nodal form of the coordinates.
#include <math.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_conf_nodal.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_pkpm_lbo_moms.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#define NB 8 // blocks of the distribution functions
#define NM 6 // moments [M0, M2, M0_G, M1, M3, M1_G]
#define NCELL 3 // configuration cells per direction
#define NVEL 4 // velocity cells

static double
rnd(void)
{
  return 2.0 * rand() / (double)RAND_MAX - 1.0;
}

static void
test_moms(int cdim)
{
  int pdim = cdim + 1, nc = 1 << cdim, np = 3 * nc;
  double mean = pow(sqrt(2.0), cdim); // expansion coefficient of the constant 1
  double vel_bounds[2] = {-3.0, 4.0};

  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, 1);

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM], clo[GKYL_MAX_DIM], cup[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    lower[d] = 0.0;
    upper[d] = 1.0 + d;
    cells[d] = NCELL;
    clo[d] = 1;
    cup[d] = NCELL;
  }
  lower[cdim] = vel_bounds[0];
  upper[cdim] = vel_bounds[1];
  cells[cdim] = NVEL;
  clo[cdim] = 1;
  cup[cdim] = NVEL;
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, pdim, lower, upper, cells);
  struct gkyl_range crange, prange;
  gkyl_range_init(&crange, cdim, clo, cup);
  gkyl_range_init(&prange, pdim, clo, cup);

  struct gkyl_array *moms = gkyl_array_new(GKYL_DOUBLE, NM * nc, crange.volume);
  struct gkyl_array *fin = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *fnodal = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *prim = gkyl_array_new(GKYL_DOUBLE, 2 * nc, crange.volume);
  struct gkyl_array *coord = gkyl_array_new(GKYL_DOUBLE, cdim * nc, crange.volume);
  struct gkyl_array *coord_nodal = gkyl_array_new(GKYL_DOUBLE, cdim * nc, crange.volume);

  // Expansion of each logical coordinate of a cell (2-point Gauss-Legendre is exact for it).
  double coord_modal[GKYL_MAX_CDIM * 8] = {0.0};
  int node_lower[GKYL_MAX_CDIM] = {0, 0, 0}, node_upper[GKYL_MAX_CDIM] = {1, 1, 1};
  struct gkyl_range node_range;
  gkyl_range_init(&node_range, cdim, node_lower, node_upper);
  struct gkyl_range_iter node_iter;
  gkyl_range_iter_init(&node_iter, &node_range);
  while (gkyl_range_iter_next(&node_iter)) {
    double z[GKYL_MAX_CDIM], basis[8];
    for (int d = 0; d < cdim; ++d) {
      z[d] = gkyl_gauss_ordinates_2[node_iter.idx[d]];
    }
    cbasis.eval(z, basis);
    for (int d = 0; d < cdim; ++d) {
      for (int k = 0; k < nc; ++k) {
        coord_modal[d * nc + k] += z[d] * basis[k];
      }
    }
  }

  // Moments near [1, 1, 1, 0] with slopes that keep every corner value positive, and
  // distributions small at the velocity edges.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&crange, iter.idx);
    double *moms_d = gkyl_array_fetch(moms, loc);
    for (int i = 0; i < nc; ++i) {
      double slope = 0.1 * rnd();
      moms_d[i] = (i == 0 ? mean : 0.0) + slope; // M0
      moms_d[nc + i] = (i == 0 ? mean : 0.0) + 0.1 * rnd(); // M2
      moms_d[2 * nc + i] = (i == 0 ? mean : 0.0) + 0.1 * rnd(); // M0_G
      moms_d[3 * nc + i] = 0.2 * rnd(); // M1
      moms_d[4 * nc + i] = rnd(); // M3 (not read)
      moms_d[5 * nc + i] = rnd(); // M1_G (not read)
    }
    double *coord_d = gkyl_array_fetch(coord, loc);
    for (int i = 0; i < cdim * nc; ++i) {
      coord_d[i] = coord_modal[i];
    }
  }
  gkyl_range_iter_init(&iter, &prange);
  while (gkyl_range_iter_next(&iter)) {
    double *f_d = gkyl_array_fetch(fin, gkyl_range_idx(&prange, iter.idx));
    for (int i = 0; i < NB * np; ++i) {
      f_d[i] = 0.02 * rnd();
    }
  }

  struct gkyl_pkpm_lbo_moms *up = gkyl_pkpm_lbo_moms_new(&cbasis, &pbasis, vel_bounds, false);
  struct gkyl_dg_conf_nodal *conf_nodal = gkyl_dg_conf_nodal_new(&cbasis, &pbasis, false);
  gkyl_dg_conf_nodal_phase_to_nodal(conf_nodal, &prange, fin, fnodal);
  gkyl_pkpm_lbo_moms_advance(up, &crange, &prange, moms, fnodal, prim);
  gkyl_dg_conf_nodal_conf_to_nodal(conf_nodal, &crange, coord, coord_nodal);

  double err_res = 0.0, max_corr = 0.0;
  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&crange, iter.idx);
    const double *moms_d = gkyl_array_cfetch(moms, loc);
    const double *prim_d = gkyl_array_cfetch(prim, loc);
    const double *coord_d = gkyl_array_cfetch(coord_nodal, loc);
    // F in the lowest and the highest velocity cell of this configuration cell
    int pidx[GKYL_MAX_DIM];
    for (int d = 0; d < cdim; ++d) {
      pidx[d] = iter.idx[d];
    }
    pidx[cdim] = prange.lower[cdim];
    const double *f_lo = gkyl_array_cfetch(fin, gkyl_range_idx(&prange, pidx));
    pidx[cdim] = prange.upper[cdim];
    const double *f_up = gkyl_array_cfetch(fin, gkyl_range_idx(&prange, pidx));
    for (int q = 0; q < nc; ++q) {
      double z[GKYL_MAX_DIM];
      for (int d = 0; d < cdim; ++d) {
        z[d] = coord_d[d * nc + q];
      }
      double m0 = cbasis.eval_expand(z, &moms_d[0]), m2 = cbasis.eval_expand(z, &moms_d[nc]);
      double m0_g = cbasis.eval_expand(z, &moms_d[2 * nc]);
      double m1 = cbasis.eval_expand(z, &moms_d[3 * nc]);
      z[cdim] = -1.0;
      double edge_lo = pbasis.eval_expand(z, f_lo);
      z[cdim] = 1.0;
      double edge_up = pbasis.eval_expand(z, f_up);
      double corr_m = edge_up - edge_lo;
      double corr_e = vel_bounds[1] * edge_up - vel_bounds[0] * edge_lo;
      double u = prim_d[q], vtsq = prim_d[nc + q];
      err_res = fmax(err_res, fabs(m0 * u - corr_m * vtsq - m1));
      err_res = fmax(err_res, fabs(m1 * u + (3.0 * m0 - corr_e) * vtsq - (m2 + 2.0 * m0_g)));
      max_corr = fmax(max_corr, fmax(fabs(corr_m), fabs(corr_e)));
    }
  }
  TEST_CHECK(max_corr > 1.0e-3);
  TEST_MSG("largest correction for the velocity extent %g", max_corr);
  TEST_CHECK(err_res < 1.0e-13);
  TEST_MSG("residual of the node systems %g", err_res);

  gkyl_pkpm_lbo_moms_release(up);
  gkyl_dg_conf_nodal_release(conf_nodal);
  gkyl_array_release(fnodal);
  gkyl_array_release(moms);
  gkyl_array_release(fin);
  gkyl_array_release(prim);
  gkyl_array_release(coord);
  gkyl_array_release(coord_nodal);
}

static void
test_moms_1x(void)
{
  test_moms(1);
}

static void
test_moms_2x(void)
{
  test_moms(2);
}

static void
test_moms_3x(void)
{
  test_moms(3);
}

TEST_LIST =
  {{"moms_1x", test_moms_1x}, {"moms_2x", test_moms_2x}, {"moms_3x", test_moms_3x}, {NULL, NULL}};
