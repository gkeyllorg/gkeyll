// Tests for the Vlasov LBO boundary-correction moment-type constructor.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_mom_type.h>
#include <gkyl_mom_bcorr_lbo_vlasov.h>
#include <gkyl_mom_bcorr_lbo_vlasov_priv.h>

void
test_bcorr_1x1v()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  // vBoundary holds [lower_v0, ..., upper_v0, ...] (2*vdim entries)
  double vBoundary[] = { -3.0, 3.0 };
  struct gkyl_mom_type *momt = gkyl_mom_bcorr_lbo_vlasov_new(&cbasis, &pbasis, vBoundary, false);

  TEST_CHECK( momt->cdim == cdim );
  TEST_CHECK( momt->pdim == pdim );
  TEST_CHECK( momt->poly_order == poly_order );
  TEST_CHECK( momt->num_config == cbasis.num_basis );
  TEST_CHECK( momt->num_phase == pbasis.num_basis );
  // boundary correction produces (vdim+1) moments (momentum + energy)
  TEST_CHECK( momt->num_mom == vdim+1 );
  TEST_CHECK( momt->on_dev == momt );

  struct mom_type_bcorr_lbo_vlasov *bcorr = container_of(momt, struct mom_type_bcorr_lbo_vlasov, momt);
  TEST_CHECK( bcorr->vBoundary[0] == -3.0 );
  TEST_CHECK( bcorr->vBoundary[vdim] == 3.0 );
  TEST_CHECK( bcorr->kernel != 0 );

  gkyl_mom_type_release(momt);
}

void
test_bcorr_1x2v()
{
  int cdim = 1, vdim = 2, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  double vBoundary[] = { -2.0, -4.0, 2.0, 4.0 };
  struct gkyl_mom_type *momt = gkyl_mom_bcorr_lbo_vlasov_new(&cbasis, &pbasis, vBoundary, false);

  TEST_CHECK( momt->num_mom == vdim+1 );  // 3

  struct mom_type_bcorr_lbo_vlasov *bcorr = container_of(momt, struct mom_type_bcorr_lbo_vlasov, momt);
  TEST_CHECK( bcorr->vBoundary[0] == -2.0 );
  TEST_CHECK( bcorr->vBoundary[1] == -4.0 );
  TEST_CHECK( bcorr->vBoundary[vdim+0] == 2.0 );
  TEST_CHECK( bcorr->vBoundary[vdim+1] == 4.0 );

  gkyl_mom_type_release(momt);
}

void
test_bcorr_1x3v()
{
  int cdim = 1, vdim = 3, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  double vBoundary[] = { -1.0, -2.0, -3.0, 1.0, 2.0, 3.0 };
  struct gkyl_mom_type *momt = gkyl_mom_bcorr_lbo_vlasov_new(&cbasis, &pbasis, vBoundary, false);

  TEST_CHECK( momt->num_mom == vdim+1 );  // 4

  struct mom_type_bcorr_lbo_vlasov *bcorr = container_of(momt, struct mom_type_bcorr_lbo_vlasov, momt);
  TEST_CHECK( bcorr->vBoundary[2] == -3.0 );
  TEST_CHECK( bcorr->vBoundary[vdim+2] == 3.0 );

  gkyl_mom_type_release(momt);
}

TEST_LIST = {
  { "bcorr_1x1v", test_bcorr_1x1v },
  { "bcorr_1x2v", test_bcorr_1x2v },
  { "bcorr_1x3v", test_bcorr_1x3v },
  { NULL, NULL },
};
