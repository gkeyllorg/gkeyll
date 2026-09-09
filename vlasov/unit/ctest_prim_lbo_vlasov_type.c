// Tests for the Vlasov LBO primitive-moment type object constructor.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_vlasov.h>

static void
check_prim(struct gkyl_prim_lbo_type *prim, int cdim, int pdim, int poly_order,
  int num_config, int num_phase, int udim)
{
  TEST_CHECK( prim->cdim == cdim );
  TEST_CHECK( prim->pdim == pdim );
  TEST_CHECK( prim->poly_order == poly_order );
  TEST_CHECK( prim->num_config == num_config );
  TEST_CHECK( prim->num_phase == num_phase );
  TEST_CHECK( prim->udim == udim );
  TEST_CHECK( prim->self_prim != 0 );
  TEST_CHECK( prim->cross_prim != 0 );
  TEST_CHECK( prim->on_dev == prim );
}

void
test_prim_1x1v()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  struct gkyl_prim_lbo_type *prim = gkyl_prim_lbo_vlasov_new(&cbasis, &pbasis, false);
  // udim is the velocity dimension
  check_prim(prim, cdim, pdim, poly_order, cbasis.num_basis, pbasis.num_basis, vdim);

  gkyl_prim_lbo_type_release(prim);
}

void
test_prim_1x2v()
{
  int cdim = 1, vdim = 2, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  struct gkyl_prim_lbo_type *prim = gkyl_prim_lbo_vlasov_new(&cbasis, &pbasis, false);
  check_prim(prim, cdim, pdim, poly_order, cbasis.num_basis, pbasis.num_basis, vdim);

  gkyl_prim_lbo_type_release(prim);
}

void
test_prim_1x3v()
{
  int cdim = 1, vdim = 3, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  struct gkyl_prim_lbo_type *prim = gkyl_prim_lbo_vlasov_new(&cbasis, &pbasis, false);
  check_prim(prim, cdim, pdim, poly_order, cbasis.num_basis, pbasis.num_basis, vdim);

  gkyl_prim_lbo_type_release(prim);
}

void
test_prim_2x2v()
{
  int cdim = 2, vdim = 2, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  struct gkyl_prim_lbo_type *prim = gkyl_prim_lbo_vlasov_new(&cbasis, &pbasis, false);
  check_prim(prim, cdim, pdim, poly_order, cbasis.num_basis, pbasis.num_basis, vdim);

  gkyl_prim_lbo_type_release(prim);
}

TEST_LIST = {
  { "prim_1x1v", test_prim_1x1v },
  { "prim_1x2v", test_prim_1x2v },
  { "prim_1x3v", test_prim_1x3v },
  { "prim_2x2v", test_prim_2x2v },
  { NULL, NULL },
};
