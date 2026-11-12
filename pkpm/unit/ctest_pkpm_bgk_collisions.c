// Tests for the nodal PKPM BGK operator: the explicit term must equal the Galerkin product of the
// collision frequency and (J f_t - f) block by block, relax only F and G toward the target, keep
// int M0_perp dw = 0, and add the largest node rate to the CFL rate. The implicit step must change
// f by the Galerkin product of (1 - S) and (J f_t - f).
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_conf_nodal.h>
#include <gkyl_pkpm_bgk_collisions.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#define NB 8
#define NV 3

static double
rnd(void)
{
  return 2.0 * rand() / (double)RAND_MAX - 1.0;
}

// With implicit set, the node values of rate are the shape factor S, and the product is taken with
// 1 - S.
static void
test_bgk(int cdim, bool implicit)
{
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, 1);
  int nc = cbasis.num_basis, np = pbasis.num_basis;

  int clo[3] = {1, 1, 1}, cup[3] = {2, 2, 2}, vlo[1] = {1}, vup[1] = {4};
  int plo[4] = {1, 1, 1, 1}, pup[4] = {2, 2, 2, 2};
  pup[cdim] = vup[0];
  struct gkyl_range crange, vrange, prange;
  gkyl_range_init(&crange, cdim, clo, cup);
  gkyl_range_init(&vrange, 1, vlo, vup);
  gkyl_range_init(&prange, cdim + 1, plo, pup);

  struct gkyl_array *rate = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct gkyl_array *nu_modal = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct gkyl_array *factor = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct gkyl_array *jacob = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct gkyl_array *f_lte = gkyl_array_new(GKYL_DOUBLE, NV, vrange.volume);
  struct gkyl_array *fin = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *rhs = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *fnodal = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *rhs_nodal = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *target_nodal = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *target = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *cfl = gkyl_array_new(GKYL_DOUBLE, 1, prange.volume);
  struct gkyl_array *block = gkyl_array_new(GKYL_DOUBLE, np, prange.volume);
  struct gkyl_array *prod = gkyl_array_new(GKYL_DOUBLE, np, prange.volume);
  gkyl_array_clear(rhs, 0.0);
  gkyl_array_clear(rhs_nodal, 0.0);
  gkyl_array_clear(target_nodal, 0.0);
  gkyl_array_clear(target, 0.0);
  gkyl_array_clear(cfl, 0.0);

  struct gkyl_dg_conf_nodal *conf_nodal = gkyl_dg_conf_nodal_new(&cbasis, &pbasis, false);

  double J = 1.7;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&crange, iter.idx);
    double *r = gkyl_array_fetch(rate, loc);
    double *fac = gkyl_array_fetch(factor, loc);
    for (int q = 0; q < nc; ++q) {
      r[q] = implicit ? 0.5 * (1.0 + rnd()) : 1.0 + 0.5 * fabs(rnd());
      fac[q] = implicit ? 1.0 - r[q] : r[q];
    }
    double *j = gkyl_array_fetch(jacob, loc);
    j[0] = J;
  }
  gkyl_dg_conf_nodal_conf_to_modal(conf_nodal, &crange, factor, nu_modal);
  gkyl_range_iter_init(&iter, &vrange);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(f_lte, gkyl_range_idx(&vrange, iter.idx));
    for (int j = 0; j < NV; ++j) {
      d[j] = rnd();
    }
  }
  gkyl_range_iter_init(&iter, &prange);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(fin, gkyl_range_idx(&prange, iter.idx));
    for (int k = 0; k < NB * np; ++k) {
      d[k] = rnd();
    }
  }

  // Explicit: convert, accumulate the increment at the nodes, convert back.
  struct gkyl_pkpm_bgk_collisions *up = gkyl_pkpm_bgk_collisions_new(&cbasis, &pbasis, false);
  gkyl_dg_conf_nodal_phase_to_nodal(conf_nodal, &prange, fin, fnodal);
  if (implicit) {
    // in place on the expansions: the change of f is the right-hand side
    gkyl_array_copy(rhs, fin);
    gkyl_pkpm_bgk_collisions_implicit_advance(
      up, &crange, &vrange, &prange, rate, jacob, f_lte, rhs
    );
    gkyl_array_accumulate(rhs, -1.0, fin);
  } else {
    gkyl_pkpm_bgk_collisions_advance(
      up, &crange, &vrange, &prange, rate, jacob, f_lte, fnodal, rhs_nodal, cfl
    );
    gkyl_dg_conf_nodal_phase_to_modal_accumulate(conf_nodal, &prange, rhs_nodal, rhs);
  }

  // The modal target [J f_t, J f_t, 0, ...], uniform in configuration space.
  gkyl_range_iter_init(&iter, &prange);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(&prange, iter.idx);
    const double *ft = gkyl_array_cfetch(f_lte, gkyl_range_idx(&vrange, &iter.idx[cdim]));
    double *tn = gkyl_array_fetch(target_nodal, ploc);
    for (int b = 0; b < 2; ++b) {
      for (int q = 0; q < nc; ++q) {
        for (int j = 0; j < NV; ++j) {
          tn[b * np + q * NV + j] = J * ft[j];
        }
      }
    }
  }
  gkyl_dg_conf_nodal_phase_to_modal_accumulate(conf_nodal, &prange, target_nodal, target);

  // The reference: block by block, the Galerkin product of nu (or 1 - S) and (J f_t - f).
  for (int b = 0; b < NB; ++b) {
    gkyl_range_iter_init(&iter, &prange);
    while (gkyl_range_iter_next(&iter)) {
      long ploc = gkyl_range_idx(&prange, iter.idx);
      const double *td = gkyl_array_cfetch(target, ploc);
      const double *fd = gkyl_array_cfetch(fin, ploc);
      double *bd = gkyl_array_fetch(block, ploc);
      for (int k = 0; k < np; ++k) {
        bd[k] = td[b * np + k] - fd[b * np + k];
      }
    }
    gkyl_dg_mul_conf_phase_op_range(&cbasis, &pbasis, prod, nu_modal, block, &crange, &prange);

    gkyl_range_iter_init(&iter, &prange);
    while (gkyl_range_iter_next(&iter)) {
      long ploc = gkyl_range_idx(&prange, iter.idx);
      const double *rd = gkyl_array_cfetch(rhs, ploc), *pd = gkyl_array_cfetch(prod, ploc);
      for (int k = 0; k < np; ++k) {
        TEST_CHECK(fabs(rd[b * np + k] - pd[k]) < 1e-13);
      }
    }
  }

  gkyl_range_iter_init(&iter, &prange);
  while (gkyl_range_iter_next(&iter)) {
    const double *r = gkyl_array_cfetch(rate, gkyl_range_idx(&crange, iter.idx));
    double rmax = 0.0;
    for (int q = 0; q < nc; ++q) {
      rmax = fmax(rmax, r[q]);
    }
    const double *c = gkyl_array_cfetch(cfl, gkyl_range_idx(&prange, iter.idx));
    TEST_CHECK(implicit ? c[0] == 0.0 : gkyl_compare_double(c[0], rmax, 1e-15));
  }

  gkyl_pkpm_bgk_collisions_release(up);
  gkyl_dg_conf_nodal_release(conf_nodal);
  gkyl_array_release(factor);
  gkyl_array_release(fnodal);
  gkyl_array_release(rhs_nodal);
  gkyl_array_release(target_nodal);
  gkyl_array_release(target);
  gkyl_array_release(rate);
  gkyl_array_release(nu_modal);
  gkyl_array_release(jacob);
  gkyl_array_release(f_lte);
  gkyl_array_release(fin);
  gkyl_array_release(rhs);
  gkyl_array_release(cfl);
  gkyl_array_release(block);
  gkyl_array_release(prod);
}

void
test_explicit_1x()
{
  test_bgk(1, false);
}

void
test_explicit_2x()
{
  test_bgk(2, false);
}

void
test_explicit_3x()
{
  test_bgk(3, false);
}

void
test_implicit_1x()
{
  test_bgk(1, true);
}

void
test_implicit_2x()
{
  test_bgk(2, true);
}

void
test_implicit_3x()
{
  test_bgk(3, true);
}

TEST_LIST = {
  {"explicit_1x", test_explicit_1x},
  {"explicit_2x", test_explicit_2x},
  {"explicit_3x", test_explicit_3x},
  {"implicit_1x", test_implicit_1x},
  {"implicit_2x", test_implicit_2x},
  {"implicit_3x", test_implicit_3x},
  {NULL, NULL}
};
