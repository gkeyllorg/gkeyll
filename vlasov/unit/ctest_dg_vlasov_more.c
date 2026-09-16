// Additional tests for the Vlasov DG equation object constructor,
// covering configurations not exercised by ctest_dg_vlasov.c
// (streaming-only field, higher dimensions, internal pointers).
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_priv.h>

void
test_vlasov_null_field()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 2);
  gkyl_cart_modal_serendip(&pbasis, pdim, 2);

  struct gkyl_range crange, prange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 20 });
  gkyl_range_init_from_shape(&prange, pdim, (int[]) { 20, 10 });

  // Streaming-only: field_id = NULL uses the streaming volume kernel.
  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_new(&cbasis, &pbasis, &crange, &prange,
    GKYL_MODEL_DEFAULT, GKYL_FIELD_NULL, false);

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != 0 );
  TEST_CHECK( eqn->surf_term != 0 );
  TEST_CHECK( eqn->on_dev == eqn );

  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);
  TEST_CHECK( vlasov->cdim == cdim );
  TEST_CHECK( vlasov->pdim == pdim );
  TEST_CHECK( vlasov->conf_range.volume == 20 );

  gkyl_dg_eqn_release(eqn);
}

void
test_vlasov_2x3v()
{
  int cdim = 2, vdim = 3, pdim = cdim+vdim;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);

  struct gkyl_range crange, prange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 6, 4 });
  gkyl_range_init_from_shape(&prange, pdim, (int[]) { 6, 4, 2, 2, 2 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_new(&cbasis, &pbasis, &crange, &prange,
    GKYL_MODEL_DEFAULT, GKYL_FIELD_E_B, false);

  TEST_CHECK( eqn->num_equations == 1 );

  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);
  TEST_CHECK( vlasov->cdim == 2 );
  TEST_CHECK( vlasov->pdim == 5 );
  TEST_CHECK( vlasov->conf_range.volume == 24 );

  gkyl_dg_eqn_release(eqn);
}

void
test_vlasov_tensor_basis()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_tensor(&cbasis, cdim, 2);
  gkyl_cart_modal_tensor(&pbasis, pdim, 2);

  struct gkyl_range crange, prange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 12 });
  gkyl_range_init_from_shape(&prange, pdim, (int[]) { 12, 8 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_new(&cbasis, &pbasis, &crange, &prange,
    GKYL_MODEL_DEFAULT, GKYL_FIELD_E_B, false);

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != 0 );

  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);
  TEST_CHECK( vlasov->cdim == cdim );
  TEST_CHECK( vlasov->pdim == pdim );

  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "vlasov_null_field", test_vlasov_null_field },
  { "vlasov_2x3v", test_vlasov_2x3v },
  { "vlasov_tensor_basis", test_vlasov_tensor_basis },
  { NULL, NULL },
};
