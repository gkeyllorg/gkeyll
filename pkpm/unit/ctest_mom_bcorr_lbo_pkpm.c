// Test construction of PKPM LBO boundary-correction moment objects.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_mom_bcorr_lbo_pkpm.h>
#include <gkyl_mom_bcorr_lbo_pkpm_kernels.h>
#include <gkyl_mom_type.h>

void
test_bcorr_1x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order);

  double vBoundary[] = {-6.0, 6.0};
  struct gkyl_mom_type *m = gkyl_mom_bcorr_lbo_pkpm_new(&cbasis, &pbasis, vBoundary, 1.0, false);

  TEST_CHECK(m->cdim == 1);
  TEST_CHECK(m->pdim == 2);
  TEST_CHECK(m->poly_order == 1);
  TEST_CHECK(m->num_config == cbasis.num_basis);
  TEST_CHECK(m->num_phase == pbasis.num_basis);
  // pkpm in local rest frame: only energy correction => 2 components.
  TEST_CHECK(m->num_mom == 2);
  TEST_CHECK(m->kernel != NULL);

  gkyl_mom_type_release(m);
}

void
test_bcorr_2x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 3, poly_order);

  double vBoundary[] = {-6.0, 6.0};
  struct gkyl_mom_type *m = gkyl_mom_bcorr_lbo_pkpm_new(&cbasis, &pbasis, vBoundary, 1.0, false);
  TEST_CHECK(m->cdim == 2);
  TEST_CHECK(m->pdim == 3);
  TEST_CHECK(m->num_mom == 2);
  gkyl_mom_type_release(m);
}

// The corrections of the two velocity edges accumulate independently: each edge weights its own
// trace by its own vBoundary, in whichever order the edges are visited.
void
test_bcorr_edges_additive()
{
  double f_up[6] = {0.3, 0.01, 0.05, 0.002, 0.01, 0.001};
  double f_lo[6] = {0.1, 0.02, -0.03, 0.001, 0.02, 0.002};
  double v_bnd[2] = {-6.0, 6.0}, dxv[2] = {1.0, 0.5};
  int idx[2] = {0, 0};

  double both[4] = {0.0}, up[4] = {0.0}, lo[4] = {0.0};
  mom_bcorr_lbo_pkpm_1x1v_ser_p1(idx, GKYL_VX_UPPER, v_bnd, dxv, 1.0, f_up, both);
  mom_bcorr_lbo_pkpm_1x1v_ser_p1(idx, GKYL_VX_LOWER, v_bnd, dxv, 1.0, f_lo, both);
  mom_bcorr_lbo_pkpm_1x1v_ser_p1(idx, GKYL_VX_UPPER, v_bnd, dxv, 1.0, f_up, up);
  mom_bcorr_lbo_pkpm_1x1v_ser_p1(idx, GKYL_VX_LOWER, v_bnd, dxv, 1.0, f_lo, lo);

  for (int c = 0; c < 4; ++c) {
    TEST_CHECK(gkyl_compare_double(both[c], up[c] + lo[c], 1e-14));
    TEST_MSG("component %d: got %g, expected %g", c, both[c], up[c] + lo[c]);
  }
  // The energy corrections are the number corrections weighted by the edge velocity.
  TEST_CHECK(gkyl_compare_double(up[2], v_bnd[1] * up[0], 1e-14));
  TEST_CHECK(gkyl_compare_double(lo[2], v_bnd[0] * lo[0], 1e-14));
}

TEST_LIST = {
  {"bcorr_1x1v_p1", test_bcorr_1x1v_p1},
  {"bcorr_2x1v_p1", test_bcorr_2x1v_p1},
  {"bcorr_edges_additive", test_bcorr_edges_additive},
  {NULL, NULL}
};
