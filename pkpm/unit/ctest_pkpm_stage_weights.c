// The weighted stage combination of the normalized PKPM distributions: against the same
// combination formed at the configuration-space nodes, with a temperature factor that differs
// at every node; the temperature factor of the combined state; the SSP-RK3 combination left
// untouched when the factor is 1; and the update of the factor by an implicit step.
#include <math.h>
#include <stdlib.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_basis.h>
#include <gkyl_dg_conf_nodal.h>
#include <gkyl_pkpm_stage_weights.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#define NB 8 // blocks of the distributions
#define NV 3 // velocity modes per node

static double
rnd(void)
{
  return 2.0 * rand() / (double)RAND_MAX - 1.0;
}

static void
test_weights(int cdim)
{
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, 1);
  int nc = cbasis.num_basis, np = pbasis.num_basis;

  int clo[3] = {1, 1, 1}, cup[3] = {2, 2, 2}, plo[4] = {1, 1, 1, 1}, pup[4] = {2, 2, 2, 2};
  pup[cdim] = 3;
  struct gkyl_range crange, prange;
  gkyl_range_init(&crange, cdim, clo, cup);
  gkyl_range_init(&prange, cdim + 1, plo, pup);

  struct gkyl_array *f_old = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *f_fe = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *f = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *f_plain = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *old_nodal = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *fe_nodal = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *ref_nodal = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *ref = gkyl_array_new(GKYL_DOUBLE, NB * np, prange.volume);
  struct gkyl_array *lam = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct gkyl_array *lam0 = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);
  struct gkyl_array *ratio = gkyl_array_new(GKYL_DOUBLE, nc, crange.volume);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &prange);
  while (gkyl_range_iter_next(&iter)) {
    long ploc = gkyl_range_idx(&prange, iter.idx);
    double *old_d = gkyl_array_fetch(f_old, ploc), *fe_d = gkyl_array_fetch(f_fe, ploc);
    for (int k = 0; k < NB * np; ++k) {
      old_d[k] = rnd();
      fe_d[k] = rnd();
    }
  }
  // temperature factors from strong cooling (0.05) to strong heating (150), node by node
  gkyl_range_iter_init(&iter, &crange);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&crange, iter.idx);
    double *lam_d = gkyl_array_fetch(lam0, loc), *ratio_d = gkyl_array_fetch(ratio, loc);
    for (int q = 0; q < nc; ++q) {
      lam_d[q] = pow(10.0, 1.75 * rnd() + 0.4);
      ratio_d[q] = 0.5 + 0.4 * rnd();
    }
  }

  struct gkyl_dg_conf_nodal *conf_nodal = gkyl_dg_conf_nodal_new(&cbasis, &pbasis, false);
  struct gkyl_pkpm_stage_weights *up = gkyl_pkpm_stage_weights_new(&cbasis, &pbasis, false);
  gkyl_dg_conf_nodal_phase_to_nodal(conf_nodal, &prange, f_old, old_nodal);
  gkyl_dg_conf_nodal_phase_to_nodal(conf_nodal, &prange, f_fe, fe_nodal);

  double coefs[2] = {1.0 / 4.0, 2.0 / 3.0};
  for (int ic = 0; ic < 2; ++ic) {
    double coef_fe = coefs[ic], coef_old = 1.0 - coef_fe;
    gkyl_array_set(f_plain, coef_old, f_old);
    gkyl_array_accumulate(f_plain, coef_fe, f_fe);

    // lam = 1: the SSP-RK3 combination is not changed, and the factor stays 1
    gkyl_array_copy(f, f_plain);
    gkyl_array_clear(lam, 1.0);
    gkyl_pkpm_stage_weights_advance(up, &crange, &prange, coef_fe, f_old, f_fe, lam, f);
    bool is_same = true;
    gkyl_range_iter_init(&iter, &prange);
    while (gkyl_range_iter_next(&iter)) {
      long ploc = gkyl_range_idx(&prange, iter.idx);
      const double *f_d = gkyl_array_cfetch(f, ploc), *plain_d = gkyl_array_cfetch(f_plain, ploc);
      for (int k = 0; k < NB * np; ++k) {
        is_same = is_same && f_d[k] == plain_d[k];
      }
    }
    TEST_CHECK(is_same);
    gkyl_range_iter_init(&iter, &crange);
    while (gkyl_range_iter_next(&iter)) {
      const double *lam_d = gkyl_array_cfetch(lam, gkyl_range_idx(&crange, iter.idx));
      for (int q = 0; q < nc; ++q) {
        TEST_CHECK(lam_d[q] == 1.0);
      }
    }

    // the weighted combination against the one formed at the nodes
    gkyl_array_copy(f, f_plain);
    gkyl_array_copy(lam, lam0);
    gkyl_pkpm_stage_weights_advance(up, &crange, &prange, coef_fe, f_old, f_fe, lam, f);
    gkyl_range_iter_init(&iter, &prange);
    while (gkyl_range_iter_next(&iter)) {
      long ploc = gkyl_range_idx(&prange, iter.idx);
      const double *lam_d = gkyl_array_cfetch(lam0, gkyl_range_idx(&crange, iter.idx));
      const double *old_d = gkyl_array_cfetch(old_nodal, ploc),
                   *fe_d = gkyl_array_cfetch(fe_nodal, ploc);
      double *ref_d = gkyl_array_fetch(ref_nodal, ploc);
      for (int b = 0; b < NB; ++b) {
        for (int q = 0; q < nc; ++q) {
          double w_fe = coef_fe * lam_d[q] / (coef_old + coef_fe * lam_d[q]);
          TEST_CHECK(w_fe >= 0.0 && w_fe <= 1.0);
          for (int j = 0; j < NV; ++j) {
            int k = b * np + q * NV + j;
            ref_d[k] = (1.0 - w_fe) * old_d[k] + w_fe * fe_d[k];
          }
        }
      }
    }
    gkyl_array_clear(ref, 0.0);
    gkyl_dg_conf_nodal_phase_to_modal_accumulate(conf_nodal, &prange, ref_nodal, ref);
    double err = 0.0;
    gkyl_range_iter_init(&iter, &prange);
    while (gkyl_range_iter_next(&iter)) {
      long ploc = gkyl_range_idx(&prange, iter.idx);
      const double *f_d = gkyl_array_cfetch(f, ploc), *ref_d = gkyl_array_cfetch(ref, ploc);
      for (int k = 0; k < NB * np; ++k) {
        err = fmax(err, fabs(f_d[k] - ref_d[k]));
      }
    }
    TEST_CHECK(err < 1.0e-13);
    TEST_MSG(
      "cdim %d, c_Y %g: weighted combination differs from the nodal one by %g", cdim, coef_fe, err
    );

    // temperature factor of the combined state, then an implicit step
    gkyl_pkpm_stage_weights_scale(up, &crange, ratio, lam);
    gkyl_range_iter_init(&iter, &crange);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&crange, iter.idx);
      const double *lam_d = gkyl_array_cfetch(lam, loc), *lam0_d = gkyl_array_cfetch(lam0, loc);
      const double *ratio_d = gkyl_array_cfetch(ratio, loc);
      for (int q = 0; q < nc; ++q) {
        TEST_CHECK(
          gkyl_compare_double(lam_d[q], (coef_old + coef_fe * lam0_d[q]) / ratio_d[q], 1.0e-14)
        );
      }
    }
  }

  gkyl_pkpm_stage_weights_release(up);
  gkyl_dg_conf_nodal_release(conf_nodal);
  gkyl_array_release(f_old);
  gkyl_array_release(f_fe);
  gkyl_array_release(f);
  gkyl_array_release(f_plain);
  gkyl_array_release(old_nodal);
  gkyl_array_release(fe_nodal);
  gkyl_array_release(ref_nodal);
  gkyl_array_release(ref);
  gkyl_array_release(lam);
  gkyl_array_release(lam0);
  gkyl_array_release(ratio);
}

static void
test_weights_1x(void)
{
  test_weights(1);
}

static void
test_weights_2x(void)
{
  test_weights(2);
}

static void
test_weights_3x(void)
{
  test_weights(3);
}

TEST_LIST = {
  {"weights_1x", test_weights_1x},
  {"weights_2x", test_weights_2x},
  {"weights_3x", test_weights_3x},
  {NULL, NULL}
};
