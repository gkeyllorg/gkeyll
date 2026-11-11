// The PKPM LBO in the configuration-nodal, velocity-modal form (conversion, velocity-only kernels
// at each node through hyper_dg, conversion back) against the Galerkin drag and diffusion
// updaters with the same collision frequency and primitive moments: at p=1 every weak product is
// a pointwise product at the Gauss-Legendre nodes, so the two agree to round-off on all 8 blocks,
// at interior and at edge velocity cells. The CFL rate is the largest node rate, which is not
// below the Galerkin updaters' rate from the cell averages.
#include <math.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_conf_nodal.h>
#include <gkyl_dg_updater_lbo_pkpm.h>
#include <gkyl_dg_updater_lbo_pkpm_nodal.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

#define NB 8 // blocks of the distribution functions
#define NCELL 2 // configuration cells per direction
#define NVEL 6 // velocity cells

static double
rnd(void)
{
  return 2.0 * rand() / (double)RAND_MAX - 1.0;
}

static void
test_lbo(int cdim)
{
  int pdim = cdim + 1, nc = 1 << cdim, np = 3 * nc;
  double mean = pow(sqrt(2.0), cdim); // expansion coefficient of the constant 1

  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, 1);

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  for (int d = 0; d < cdim; ++d) {
    lower[d] = 0.0;
    upper[d] = 1.0 + d;
    cells[d] = NCELL;
    ghost[d] = 1;
  }
  lower[cdim] = -2.0;
  upper[cdim] = 2.5;
  cells[cdim] = NVEL;
  ghost[cdim] = 0;
  struct gkyl_rect_grid conf_grid, phase_grid;
  gkyl_rect_grid_init(&conf_grid, cdim, lower, upper, cells);
  gkyl_rect_grid_init(&phase_grid, pdim, lower, upper, cells);
  struct gkyl_range crange, crange_ext, prange, prange_ext;
  gkyl_create_grid_ranges(&conf_grid, ghost, &crange_ext, &crange);
  gkyl_create_grid_ranges(&phase_grid, ghost, &prange_ext, &prange);

  struct gkyl_array *nu_modal = gkyl_array_new(GKYL_DOUBLE, nc, crange_ext.volume);
  struct gkyl_array *prim_modal = gkyl_array_new(GKYL_DOUBLE, 2 * nc, crange_ext.volume);
  struct gkyl_array *nu_prim_modal = gkyl_array_new(GKYL_DOUBLE, 2 * nc, crange_ext.volume);
  struct gkyl_array *nu = gkyl_array_new(GKYL_DOUBLE, nc, crange_ext.volume);
  struct gkyl_array *prim = gkyl_array_new(GKYL_DOUBLE, 2 * nc, crange_ext.volume);
  struct gkyl_array *fin = gkyl_array_new(GKYL_DOUBLE, NB * np, prange_ext.volume);
  struct gkyl_array *fnodal = gkyl_array_new(GKYL_DOUBLE, NB * np, prange_ext.volume);
  struct gkyl_array *rhs_nodal = gkyl_array_new(GKYL_DOUBLE, NB * np, prange_ext.volume);
  struct gkyl_array *rhs_ref = gkyl_array_new(GKYL_DOUBLE, NB * np, prange_ext.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, NB * np, prange_ext.volume);
  struct gkyl_array *cfl_ref = gkyl_array_new(GKYL_DOUBLE, 1, prange_ext.volume);
  struct gkyl_array *cfl = gkyl_array_new(GKYL_DOUBLE, 1, prange_ext.volume);

  // A frequency near 2, a flow inside the velocity grid (so the drag changes direction inside
  // it) and a thermal speed squared near 1.5, all varying in the cell; random distributions.
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&crange, iter.idx);
    double *nu_d = gkyl_array_fetch(nu_modal, loc);
    double *prim_d = gkyl_array_fetch(prim_modal, loc);
    for (int i = 0; i < nc; ++i) {
      nu_d[i] = (i == 0 ? 2.0 * mean : 0.0) + 0.3 * rnd();
      prim_d[i] = (i == 0 ? 0.3 * mean : 0.0) + 0.5 * rnd();
      prim_d[nc + i] = (i == 0 ? 1.5 * mean : 0.0) + 0.2 * rnd();
    }
  }
  struct gkyl_dg_conf_nodal *conf_nodal = gkyl_dg_conf_nodal_new(&cbasis, &pbasis, false);
  gkyl_dg_conf_nodal_conf_to_nodal(conf_nodal, &crange, nu_modal, nu);
  gkyl_dg_conf_nodal_conf_to_nodal(conf_nodal, &crange, prim_modal, prim);
  gkyl_dg_mul_op_range(&cbasis, 0, nu_prim_modal, 0, prim_modal, 0, nu_modal, &crange);
  gkyl_dg_mul_op_range(&cbasis, 1, nu_prim_modal, 1, prim_modal, 0, nu_modal, &crange);
  gkyl_range_iter_init(&iter, &prange);
  while (gkyl_range_iter_next(&iter)) {
    double *f_d = gkyl_array_fetch(fin, gkyl_range_idx(&prange, iter.idx));
    for (int i = 0; i < NB * np; ++i) {
      f_d[i] = rnd();
    }
  }

  // Galerkin reference: the drag and diffusion updaters.
  struct gkyl_dg_lbo_pkpm_drag_auxfields drag_inp = {
    .nuSum = nu_modal,
    .nuPrimMomsSum = nu_prim_modal,
  };
  struct gkyl_dg_lbo_pkpm_diff_auxfields diff_inp = {
    .nuSum = nu_modal,
    .nuPrimMomsSum = nu_prim_modal,
  };
  struct gkyl_dg_updater_collisions *slvr = gkyl_dg_updater_lbo_pkpm_new(
    &phase_grid, &cbasis, &pbasis, &crange, &drag_inp, &diff_inp, false
  );
  gkyl_array_clear(rhs_ref, 0.0);
  gkyl_array_clear(cfl_ref, 0.0);
  gkyl_dg_updater_lbo_pkpm_advance(slvr, &prange, fin, cfl_ref, rhs_ref);

  // Nodal form: convert, apply the operator at the nodes, convert back.
  struct gkyl_dg_lbo_pkpm_nodal_auxfields nodal_inp = {.nu = nu, .prim_moms = prim};
  struct gkyl_dg_updater_lbo_pkpm_nodal *up =
    gkyl_dg_updater_lbo_pkpm_nodal_new(&phase_grid, &cbasis, &pbasis, &crange, &nodal_inp, false);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(rhs_nodal, 0.0);
  gkyl_array_clear(cfl, 0.0);
  gkyl_dg_conf_nodal_phase_to_nodal(conf_nodal, &prange, fin, fnodal);
  gkyl_dg_updater_lbo_pkpm_nodal_advance(up, &prange, fnodal, cfl, rhs_nodal);
  gkyl_dg_conf_nodal_phase_to_modal_accumulate(conf_nodal, &prange, rhs_nodal, rhs);

  double err_rhs = 0.0, max_rhs = 0.0, err_cfl = 0.0, max_cfl = 0.0;
  gkyl_range_iter_init(&iter, &prange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&prange, iter.idx);
    const double *ref_d = gkyl_array_cfetch(rhs_ref, loc);
    const double *rhs_d = gkyl_array_cfetch(rhs, loc);
    for (int i = 0; i < NB * np; ++i) {
      err_rhs = fmax(err_rhs, fabs(rhs_d[i] - ref_d[i]));
      max_rhs = fmax(max_rhs, fabs(ref_d[i]));
    }
    // largest node rate: 2.5 (2/dv) nu |u - w_c| + 3 nu + 9 (2/dv)^2 nu vtsq
    long cloc = gkyl_range_idx(&crange, iter.idx);
    const double *nu_d = gkyl_array_cfetch(nu, cloc);
    const double *prim_d = gkyl_array_cfetch(prim, cloc);
    double xc[GKYL_MAX_DIM], rdv = 2.0 / phase_grid.dx[cdim], rate_max = 0.0;
    gkyl_rect_grid_cell_center(&phase_grid, iter.idx, xc);
    for (int q = 0; q < nc; ++q) {
      double rate = 2.5 * rdv * nu_d[q] * fabs(prim_d[q] - xc[cdim]) + 3.0 * nu_d[q] +
                    9.0 * rdv * rdv * nu_d[q] * prim_d[nc + q];
      rate_max = fmax(rate_max, rate);
    }
    const double *cfl_ref_d = gkyl_array_cfetch(cfl_ref, loc);
    const double *cfl_d = gkyl_array_cfetch(cfl, loc);
    err_cfl = fmax(err_cfl, fabs(cfl_d[0] - rate_max));
    max_cfl = fmax(max_cfl, rate_max);
    TEST_CHECK(cfl_d[0] >= cfl_ref_d[0] * (1.0 - 1.0e-13));
  }
  TEST_CHECK(max_rhs > 1.0);
  TEST_CHECK(err_rhs < 1.0e-13 * max_rhs);
  TEST_MSG("operator differs from the Galerkin updaters by %g (largest value %g)", err_rhs, max_rhs);
  TEST_CHECK(err_cfl < 1.0e-13 * max_cfl);
  TEST_MSG("CFL rate differs from the largest node rate by %g (largest value %g)", err_cfl, max_cfl);

  gkyl_dg_updater_lbo_pkpm_nodal_release(up);
  gkyl_dg_conf_nodal_release(conf_nodal);
  gkyl_dg_updater_lbo_pkpm_release(slvr);
  gkyl_array_release(nu_modal);
  gkyl_array_release(prim_modal);
  gkyl_array_release(nu_prim_modal);
  gkyl_array_release(nu);
  gkyl_array_release(prim);
  gkyl_array_release(fin);
  gkyl_array_release(fnodal);
  gkyl_array_release(rhs_nodal);
  gkyl_array_release(rhs_ref);
  gkyl_array_release(rhs);
  gkyl_array_release(cfl_ref);
  gkyl_array_release(cfl);
}

static void
test_lbo_1x(void)
{
  test_lbo(1);
}

static void
test_lbo_2x(void)
{
  test_lbo(2);
}

static void
test_lbo_3x(void)
{
  test_lbo(3);
}

TEST_LIST =
  {{"lbo_1x", test_lbo_1x}, {"lbo_2x", test_lbo_2x}, {"lbo_3x", test_lbo_3x}, {NULL, NULL}};
