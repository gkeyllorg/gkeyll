// Test construction of PKPM LBO boundary-correction moment objects.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_mom_bcorr_lbo_pkpm.h>
#include <gkyl_mom_type.h>

void
test_bcorr_1x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order);

  double vBoundary[] = {-6.0, 6.0};
  struct gkyl_mom_type *m =
    gkyl_mom_bcorr_lbo_pkpm_new(&cbasis, &pbasis, vBoundary, 1.0, false);

  TEST_CHECK( m->cdim == 1 );
  TEST_CHECK( m->pdim == 2 );
  TEST_CHECK( m->poly_order == 1 );
  TEST_CHECK( m->num_config == cbasis.num_basis );
  TEST_CHECK( m->num_phase == pbasis.num_basis );
  // pkpm in local rest frame: only energy correction => 2 components.
  TEST_CHECK( m->num_mom == 2 );
  TEST_CHECK( m->kernel != NULL );

  gkyl_mom_type_release(m);
}

void
test_bcorr_1x1v_p2()
{
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order);

  double vBoundary[] = {-8.0, 8.0};
  struct gkyl_mom_type *m =
    gkyl_mom_bcorr_lbo_pkpm_new(&cbasis, &pbasis, vBoundary, 2.0, false);
  TEST_CHECK( m->poly_order == 2 );
  TEST_CHECK( m->num_mom == 2 );
  TEST_CHECK( m->kernel != NULL );
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
  struct gkyl_mom_type *m =
    gkyl_mom_bcorr_lbo_pkpm_new(&cbasis, &pbasis, vBoundary, 1.0, false);
  TEST_CHECK( m->cdim == 2 );
  TEST_CHECK( m->pdim == 3 );
  TEST_CHECK( m->num_mom == 2 );
  gkyl_mom_type_release(m);
}

TEST_LIST = {
  { "bcorr_1x1v_p1", test_bcorr_1x1v_p1 },
  { "bcorr_1x1v_p2", test_bcorr_1x1v_p2 },
  { "bcorr_2x1v_p1", test_bcorr_2x1v_p1 },
  { NULL, NULL },
};
