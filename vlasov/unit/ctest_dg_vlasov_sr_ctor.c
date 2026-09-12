// Tests for the special-relativistic Vlasov DG equation object constructor.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_dg_vlasov_sr.h>
#include <gkyl_dg_vlasov_sr_priv.h>

void
test_vlasov_sr_1x1v_p2()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 2);
  gkyl_cart_modal_serendip(&pbasis, pdim, 2);

  struct gkyl_range crange, vrange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 10 });
  gkyl_range_init_from_shape(&vrange, vdim, (int[]) { 20 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_sr_new(&cbasis, &pbasis, &crange, &vrange,
    GKYL_FIELD_E_B, false);

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != 0 );
  TEST_CHECK( eqn->surf_term != 0 );
  TEST_CHECK( eqn->boundary_surf_term != 0 );
  TEST_CHECK( eqn->on_dev == eqn );

  struct dg_vlasov_sr *sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  TEST_CHECK( sr->cdim == cdim );
  TEST_CHECK( sr->pdim == pdim );
  TEST_CHECK( sr->conf_range.volume == 10 );
  TEST_CHECK( sr->vel_range.volume == 20 );
  TEST_CHECK( sr->auxfields.qmem == 0 );
  TEST_CHECK( sr->auxfields.gamma == 0 );

  gkyl_dg_eqn_release(eqn);
}

void
test_vlasov_sr_1x1v_null_field()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 2);
  gkyl_cart_modal_serendip(&pbasis, pdim, 2);

  struct gkyl_range crange, vrange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 6 });
  gkyl_range_init_from_shape(&vrange, vdim, (int[]) { 12 });

  // With a NULL field the volume term uses the streaming-only kernel,
  // but should still be a valid non-NULL pointer.
  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_sr_new(&cbasis, &pbasis, &crange, &vrange,
    GKYL_FIELD_NULL, false);

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != 0 );

  struct dg_vlasov_sr *sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  TEST_CHECK( sr->cdim == cdim );
  TEST_CHECK( sr->pdim == pdim );

  gkyl_dg_eqn_release(eqn);
}

void
test_vlasov_sr_1x2v_hybrid()
{
  int cdim = 1, vdim = 2, pdim = cdim+vdim;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);

  struct gkyl_range crange, vrange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 8 });
  gkyl_range_init_from_shape(&vrange, vdim, (int[]) { 4, 4 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_vlasov_sr_new(&cbasis, &pbasis, &crange, &vrange,
    GKYL_FIELD_E_B, false);

  TEST_CHECK( eqn->num_equations == 1 );

  struct dg_vlasov_sr *sr = container_of(eqn, struct dg_vlasov_sr, eqn);
  TEST_CHECK( sr->cdim == cdim );
  TEST_CHECK( sr->pdim == pdim );
  TEST_CHECK( sr->conf_range.volume == 8 );
  TEST_CHECK( sr->vel_range.volume == 16 );

  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "vlasov_sr_1x1v_p2", test_vlasov_sr_1x1v_p2 },
  { "vlasov_sr_1x1v_null_field", test_vlasov_sr_1x1v_null_field },
  { "vlasov_sr_1x2v_hybrid", test_vlasov_sr_1x2v_hybrid },
  { NULL, NULL },
};
