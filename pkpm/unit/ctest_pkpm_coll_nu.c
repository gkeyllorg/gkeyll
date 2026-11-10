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

// The node values of an expansion are those of the expansion evaluated at the nodes.
static void
test_to_nodal(int cdim)
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  int nc = cbasis.num_basis, num_comp = 2;
  int lower[GKYL_MAX_CDIM] = {1, 1, 1}, upper[GKYL_MAX_CDIM] = {2, 1, 1};
  struct gkyl_range range;
  gkyl_range_init(&range, cdim, lower, upper);

  struct gkyl_array *modal = gkyl_array_new(GKYL_DOUBLE, num_comp * nc, range.volume);
  struct gkyl_array *nodal = gkyl_array_new(GKYL_DOUBLE, num_comp * nc, range.volume);
  fill(modal, &range, 1.3, 0.17);

  struct gkyl_pkpm_coll_nu *up = gkyl_pkpm_coll_nu_new(&cbasis, false, 0.0, false);
  gkyl_pkpm_coll_nu_to_nodal(up, &range, modal, nodal);

  double gx = 1.0 / sqrt(3.0);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);
    const double *m = gkyl_array_cfetch(modal, loc), *n = gkyl_array_cfetch(nodal, loc);
    for (int c = 0; c < num_comp; ++c) {
      for (int q = 0; q < nc; ++q) {
        double x[3];
        for (int d = 0; d < cdim; ++d) {
          x[d] = ((q >> d) & 1) ? gx : -gx;
        }
        TEST_CHECK(gkyl_compare_double(n[c * nc + q], cbasis.eval_expand(x, &m[c * nc]), 1e-14));
      }
    }
  }

  gkyl_pkpm_coll_nu_release(up);
  gkyl_array_release(modal);
  gkyl_array_release(nodal);
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
  struct gkyl_array *n = gkyl_array_new(GKYL_DOUBLE, nc, range.volume);
  struct gkyl_array *vs = gkyl_array_new(GKYL_DOUBLE, nc, range.volume);
  struct gkyl_array *vr = gkyl_array_new(GKYL_DOUBLE, nc, range.volume);
  struct gkyl_array *nu = gkyl_array_new(GKYL_DOUBLE, nc, range.volume);
  fill(ref, &range, 2.0, 0.1);
  fill(n, &range, 1.5, 0.2);
  fill(vs, &range, 0.7, 0.05);
  fill(vr, &range, 0.3, 0.02);

  struct gkyl_pkpm_coll_nu *up_ref = gkyl_pkpm_coll_nu_new(&cbasis, false, 0.0, false);
  gkyl_pkpm_coll_nu_advance(up_ref, &range, ref, n, vs, vr, nu);
  long loc = gkyl_range_idx(&range, lower);
  const double *nu_d = gkyl_array_cfetch(nu, loc), *ref_d = gkyl_array_cfetch(ref, loc);
  for (int q = 0; q < nc; ++q) {
    TEST_CHECK(gkyl_compare_double(nu_d[q], ref_d[q], 1e-15));
  }

  double vtsq_min = 0.75;
  struct gkyl_pkpm_coll_nu *up_sp = gkyl_pkpm_coll_nu_new(&cbasis, true, vtsq_min, false);
  double *n_w = gkyl_array_fetch(n, loc);
  n_w[nc - 1] = -1.0; // a node with a negative density has no collisions
  gkyl_pkpm_coll_nu_advance(up_sp, &range, ref, n, vs, vr, nu);
  const double *vs_d = gkyl_array_cfetch(vs, loc), *vr_d = gkyl_array_cfetch(vr, loc);
  for (int q = 0; q < nc; ++q) {
    double vsum = fmax(vs_d[q], vtsq_min) + fmax(vr_d[q], vtsq_min);
    double expected = q == nc - 1 ? 0.0 : ref_d[q] * n_w[q] / pow(vsum, 1.5);
    TEST_CHECK(gkyl_compare_double(nu_d[q], expected, 1e-14));
    TEST_MSG("node %d: got %g, expected %g", q, nu_d[q], expected);
  }

  gkyl_pkpm_coll_nu_release(up_ref);
  gkyl_pkpm_coll_nu_release(up_sp);
  gkyl_array_release(ref);
  gkyl_array_release(n);
  gkyl_array_release(vs);
  gkyl_array_release(vr);
  gkyl_array_release(nu);
}

void
test_to_nodal_1x()
{
  test_to_nodal(1);
}

void
test_to_nodal_2x()
{
  test_to_nodal(2);
}

void
test_to_nodal_3x()
{
  test_to_nodal(3);
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

TEST_LIST = {{"to_nodal_1x", test_to_nodal_1x}, {"to_nodal_2x", test_to_nodal_2x},
             {"to_nodal_3x", test_to_nodal_3x}, {"advance_1x", test_advance_1x},
             {"advance_3x", test_advance_3x},   {NULL, NULL}};
