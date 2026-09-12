// Tests for the (non-relativistic) Vlasov moment-type objects, covering
// component counts (num_mom) across moment types and velocity dimensions.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_mom_type.h>
#include <gkyl_mom_vlasov.h>

void
test_mom_1x1v()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  struct gkyl_mom_type *m0 = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M0, false);
  struct gkyl_mom_type *m1 = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M1, false);
  struct gkyl_mom_type *m2 = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M2, false);
  struct gkyl_mom_type *m2ij = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M2IJ, false);
  struct gkyl_mom_type *m3 = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M3, false);
  struct gkyl_mom_type *m3ijk = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M3IJK, false);

  TEST_CHECK( m0->cdim == cdim );
  TEST_CHECK( m0->pdim == pdim );
  TEST_CHECK( m0->poly_order == poly_order );
  TEST_CHECK( m0->num_config == cbasis.num_basis );
  TEST_CHECK( m0->num_phase == pbasis.num_basis );

  TEST_CHECK( m0->num_mom == 1 );
  TEST_CHECK( m1->num_mom == vdim );                  // 1
  TEST_CHECK( m2->num_mom == 1 );
  TEST_CHECK( m2ij->num_mom == vdim*(vdim+1)/2 );     // 1
  TEST_CHECK( m3->num_mom == vdim );                  // 1
  TEST_CHECK( m3ijk->num_mom == 1 );                  // m3ijk_count[0]

  gkyl_mom_type_release(m0);
  gkyl_mom_type_release(m1);
  gkyl_mom_type_release(m2);
  gkyl_mom_type_release(m2ij);
  gkyl_mom_type_release(m3);
  gkyl_mom_type_release(m3ijk);
}

void
test_mom_1x2v()
{
  int cdim = 1, vdim = 2, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  struct gkyl_mom_type *m1 = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M1, false);
  struct gkyl_mom_type *m2ij = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M2IJ, false);
  struct gkyl_mom_type *m3 = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M3, false);
  struct gkyl_mom_type *m3ijk = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M3IJK, false);

  TEST_CHECK( m1->num_mom == vdim );               // 2
  TEST_CHECK( m2ij->num_mom == vdim*(vdim+1)/2 );  // 3
  TEST_CHECK( m3->num_mom == vdim );               // 2
  TEST_CHECK( m3ijk->num_mom == 4 );               // m3ijk_count[1]

  gkyl_mom_type_release(m1);
  gkyl_mom_type_release(m2ij);
  gkyl_mom_type_release(m3);
  gkyl_mom_type_release(m3ijk);
}

void
test_mom_1x3v()
{
  int cdim = 1, vdim = 3, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  struct gkyl_mom_type *m1 = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M1, false);
  struct gkyl_mom_type *m2ij = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M2IJ, false);
  struct gkyl_mom_type *m3 = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M3, false);
  struct gkyl_mom_type *m3ijk = gkyl_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M3IJK, false);

  TEST_CHECK( m1->num_mom == vdim );               // 3
  TEST_CHECK( m2ij->num_mom == vdim*(vdim+1)/2 );  // 6
  TEST_CHECK( m3->num_mom == vdim );               // 3
  TEST_CHECK( m3ijk->num_mom == 10 );              // m3ijk_count[2]

  gkyl_mom_type_release(m1);
  gkyl_mom_type_release(m2ij);
  gkyl_mom_type_release(m3);
  gkyl_mom_type_release(m3ijk);
}

void
test_int_mom()
{
  int cdim = 1, vdim = 2, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  // integrated (M0, M1i, M2) => 2 + vdim components
  struct gkyl_mom_type *intm = gkyl_int_mom_vlasov_new(&cbasis, &pbasis, GKYL_F_MOMENT_M0M1M2, false);
  TEST_CHECK( intm->num_mom == 2+vdim );  // 4
  TEST_CHECK( intm->cdim == cdim );
  TEST_CHECK( intm->pdim == pdim );

  gkyl_mom_type_release(intm);
}

TEST_LIST = {
  { "mom_1x1v", test_mom_1x1v },
  { "mom_1x2v", test_mom_1x2v },
  { "mom_1x3v", test_mom_1x3v },
  { "int_mom", test_int_mom },
  { NULL, NULL },
};
