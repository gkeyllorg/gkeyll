// Tests for the special-relativistic Vlasov moment-type objects (num_mom etc.).
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_mom_vlasov_sr.h>

void
test_mom_sr_1x3v()
{
  int cdim = 1, vdim = 3, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  struct gkyl_range crange, vrange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 8 });
  gkyl_range_init_from_shape(&vrange, vdim, (int[]) { 4, 4, 4 });

  struct gkyl_mom_type *m0 = gkyl_mom_vlasov_sr_new(&cbasis, &pbasis, &crange, &vrange, GKYL_F_MOMENT_M0, false);
  struct gkyl_mom_type *m1 = gkyl_mom_vlasov_sr_new(&cbasis, &pbasis, &crange, &vrange, GKYL_F_MOMENT_M1, false);
  struct gkyl_mom_type *m2 = gkyl_mom_vlasov_sr_new(&cbasis, &pbasis, &crange, &vrange, GKYL_F_MOMENT_M2, false);
  struct gkyl_mom_type *m3 = gkyl_mom_vlasov_sr_new(&cbasis, &pbasis, &crange, &vrange, GKYL_F_MOMENT_M3, false);
  struct gkyl_mom_type *ni = gkyl_mom_vlasov_sr_new(&cbasis, &pbasis, &crange, &vrange, GKYL_F_MOMENT_NI, false);
  struct gkyl_mom_type *tij = gkyl_mom_vlasov_sr_new(&cbasis, &pbasis, &crange, &vrange, GKYL_F_MOMENT_TIJ, false);

  // shared metadata
  TEST_CHECK( m0->cdim == cdim );
  TEST_CHECK( m0->pdim == pdim );
  TEST_CHECK( m0->poly_order == poly_order );
  TEST_CHECK( m0->num_config == cbasis.num_basis );
  TEST_CHECK( m0->num_phase == pbasis.num_basis );

  // number of components for each moment type
  TEST_CHECK( m0->num_mom == 1 );
  TEST_CHECK( m1->num_mom == vdim );           // 3
  TEST_CHECK( m2->num_mom == 1 );
  TEST_CHECK( m3->num_mom == vdim );           // 3
  TEST_CHECK( ni->num_mom == 1+vdim );         // 4
  TEST_CHECK( tij->num_mom == 1+vdim+(vdim*(vdim+1))/2 ); // 1+3+6 = 10

  gkyl_mom_type_release(m0);
  gkyl_mom_type_release(m1);
  gkyl_mom_type_release(m2);
  gkyl_mom_type_release(m3);
  gkyl_mom_type_release(ni);
  gkyl_mom_type_release(tij);
}

void
test_int_mom_sr_1x3v()
{
  int cdim = 1, vdim = 3, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  struct gkyl_range crange, vrange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 8 });
  gkyl_range_init_from_shape(&vrange, vdim, (int[]) { 4, 4, 4 });

  struct gkyl_mom_type *intm = gkyl_int_mom_vlasov_sr_new(&cbasis, &pbasis, &crange, &vrange,
    GKYL_F_MOMENT_M0ENERGYM3, false);

  // integrated SR moments are (M0, M2, M3i) = 2 + vdim
  TEST_CHECK( intm->num_mom == 2+vdim );  // 5
  TEST_CHECK( intm->cdim == cdim );
  TEST_CHECK( intm->pdim == pdim );

  gkyl_mom_type_release(intm);
}

TEST_LIST = {
  { "mom_sr_1x3v", test_mom_sr_1x3v },
  { "int_mom_sr_1x3v", test_int_mom_sr_1x3v },
  { NULL, NULL },
};
