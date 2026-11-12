// Tests for the collision frequency at the configuration-space Gauss-Legendre nodes.
#include <math.h>

#include <acutest.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_pkpm_coll_nu.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

static void
fill(struct gkyl_array *arr, const struct gkyl_range *range, double a, double b)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    double *d = gkyl_array_fetch(arr, gkyl_range_idx(range, iter.idx));
    for (int k = 0; k < arr->ncomp; ++k) {
      d[k] = a + b * k;
    }
  }
}

// Reference frequency, Spitzer scaling, the floor and non-positive densities.
static void
test_advance(int cdim)
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  int nc = cbasis.num_basis;
  int lower[GKYL_MAX_CDIM] = {1, 1, 1}, upper[GKYL_MAX_CDIM] = {1, 1, 1};
  struct gkyl_range range;
  gkyl_range_init(&range, cdim, lower, upper);

  struct gkyl_array *ref = gkyl_array_new(GKYL_DOUBLE, nc, range.volume);
  // node fields [n, u (3), vtsq] of the species and of its partner
  struct gkyl_array *prim_s = gkyl_array_new(GKYL_DOUBLE, 5 * nc, range.volume);
  struct gkyl_array *prim_r = gkyl_array_new(GKYL_DOUBLE, 5 * nc, range.volume);
  struct gkyl_array *nu = gkyl_array_new(GKYL_DOUBLE, nc, range.volume);
  fill(ref, &range, 2.0, 0.1);
  fill(prim_s, &range, 0.5, 0.01);
  fill(prim_r, &range, 1.5, 0.02);

  struct gkyl_pkpm_coll_nu *up_ref = gkyl_pkpm_coll_nu_new(&cbasis, false, 0.0, false);
  gkyl_pkpm_coll_nu_advance(up_ref, &range, ref, prim_s, prim_r, nu);
  long loc = gkyl_range_idx(&range, lower);
  const double *nu_d = gkyl_array_cfetch(nu, loc), *ref_d = gkyl_array_cfetch(ref, loc);
  for (int q = 0; q < nc; ++q) {
    TEST_CHECK(gkyl_compare_double(nu_d[q], ref_d[q], 1e-15));
  }

  double vtsq_min = 0.75;
  struct gkyl_pkpm_coll_nu *up_sp = gkyl_pkpm_coll_nu_new(&cbasis, true, vtsq_min, false);
  double *n_w = gkyl_array_fetch(prim_r, loc);
  n_w[nc - 1] = -1.0; // a node with a negative density has no collisions
  gkyl_pkpm_coll_nu_advance(up_sp, &range, ref, prim_s, prim_r, nu);
  const double *vs_d = (const double *)gkyl_array_cfetch(prim_s, loc) + 4 * nc;
  const double *vr_d = (const double *)gkyl_array_cfetch(prim_r, loc) + 4 * nc;
  for (int q = 0; q < nc; ++q) {
    double vsum = fmax(vs_d[q], vtsq_min) + fmax(vr_d[q], vtsq_min);
    double expected = q == nc - 1 ? 0.0 : ref_d[q] * n_w[q] / pow(vsum, 1.5);
    TEST_CHECK(gkyl_compare_double(nu_d[q], expected, 1e-14));
    TEST_MSG("node %d: got %g, expected %g", q, nu_d[q], expected);
  }

  gkyl_pkpm_coll_nu_release(up_ref);
  gkyl_pkpm_coll_nu_release(up_sp);
  gkyl_array_release(ref);
  gkyl_array_release(prim_s);
  gkyl_array_release(prim_r);
  gkyl_array_release(nu);
}

void
test_advance_1x()
{
  test_advance(1);
}

void
test_advance_3x()
{
  test_advance(3);
}

TEST_LIST = {{"advance_1x", test_advance_1x}, {"advance_3x", test_advance_3x}, {NULL, NULL}};
