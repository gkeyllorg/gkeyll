// The primitive moments of the PKPM LBO at the configuration-space nodes against the Galerkin
// solve (boundary-correction and primitive-moment calculators): at p=1 the node values of the
// Galerkin solution are the solutions of the 2x2 systems at the nodes.
#include <math.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_mom_bcorr_lbo_pkpm.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_dg_conf_nodal.h>
#include <gkyl_pkpm_lbo_moms.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_pkpm.h>
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
  struct gkyl_array *prim_ref_nodal = gkyl_array_new(GKYL_DOUBLE, 2 * nc, crange.volume);
  struct gkyl_array *nu_modal = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct gkyl_array *bcorr = gkyl_array_new(GKYL_DOUBLE, 2 * nc, crange.volume);
  struct gkyl_array *prim_ref = gkyl_array_new(GKYL_DOUBLE, 2 * nc, crange.volume);
  struct gkyl_array *prim = gkyl_array_new(GKYL_DOUBLE, 2 * nc, crange.volume);

  // Moments near [1, 1, 1, 0] with slopes that keep every corner value positive, a frequency
  // near 2, and distributions small at the velocity edges.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&crange, iter.idx);
    double *moms_d = gkyl_array_fetch(moms, loc);
    double *nu_d = gkyl_array_fetch(nu_modal, loc);
    for (int i = 0; i < nc; ++i) {
      double slope = 0.1 * rnd();
      moms_d[i] = (i == 0 ? mean : 0.0) + slope; // M0
      moms_d[nc + i] = (i == 0 ? mean : 0.0) + 0.1 * rnd(); // M2
      moms_d[2 * nc + i] = (i == 0 ? mean : 0.0) + 0.1 * rnd(); // M0_G
      moms_d[3 * nc + i] = 0.2 * rnd(); // M1
      moms_d[4 * nc + i] = rnd(); // M3 (not read)
      moms_d[5 * nc + i] = rnd(); // M1_G (not read)
      nu_d[i] = (i == 0 ? 2.0 * mean : 0.0) + 0.3 * rnd();
    }
  }
  gkyl_range_iter_init(&iter, &prange);
  while (gkyl_range_iter_next(&iter)) {
    double *f_d = gkyl_array_fetch(fin, gkyl_range_idx(&prange, iter.idx));
    for (int i = 0; i < NB * np; ++i) {
      f_d[i] = 0.02 * rnd();
    }
  }

  // Galerkin reference
  struct gkyl_mom_calc_bcorr *bcorr_calc =
    gkyl_mom_calc_bcorr_lbo_pkpm_new(&grid, &cbasis, &pbasis, vel_bounds, 1.0, false);
  struct gkyl_prim_lbo_calc *prim_calc =
    gkyl_prim_lbo_pkpm_calc_new(&grid, &cbasis, &pbasis, &crange, false);
  gkyl_mom_calc_bcorr_advance(bcorr_calc, &prange, &crange, fin, bcorr);
  gkyl_prim_lbo_calc_advance(prim_calc, &crange, moms, bcorr, nu_modal, prim_ref);

  struct gkyl_pkpm_lbo_moms *up = gkyl_pkpm_lbo_moms_new(&cbasis, &pbasis, vel_bounds, false);
  struct gkyl_dg_conf_nodal *conf_nodal = gkyl_dg_conf_nodal_new(&cbasis, &pbasis, false);
  gkyl_dg_conf_nodal_phase_to_nodal(conf_nodal, &prange, fin, fnodal);
  gkyl_pkpm_lbo_moms_advance(up, &crange, &prange, moms, fnodal, prim);
  gkyl_dg_conf_nodal_conf_to_nodal(conf_nodal, &crange, prim_ref, prim_ref_nodal);

  double err_prim = 0.0;
  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&crange, iter.idx);
    const double *ref_nodal = gkyl_array_cfetch(prim_ref_nodal, loc);
    const double *prim_d = gkyl_array_cfetch(prim, loc);
    for (int i = 0; i < 2 * nc; ++i) {
      err_prim = fmax(err_prim, fabs(prim_d[i] - ref_nodal[i]));
    }
  }
  TEST_CHECK(err_prim < 1.0e-13);
  TEST_MSG("primitive moments differ from the Galerkin solve by %g", err_prim);

  gkyl_pkpm_lbo_moms_release(up);
  gkyl_dg_conf_nodal_release(conf_nodal);
  gkyl_array_release(fnodal);
  gkyl_array_release(prim_ref_nodal);
  gkyl_mom_calc_bcorr_release(bcorr_calc);
  gkyl_prim_lbo_calc_release(prim_calc);
  gkyl_array_release(moms);
  gkyl_array_release(fin);
  gkyl_array_release(nu_modal);
  gkyl_array_release(bcorr);
  gkyl_array_release(prim_ref);
  gkyl_array_release(prim);
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
