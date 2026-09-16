// Tests for the Vlasov LBO drag DG equation object constructor.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_lbo_vlasov_drag.h>
#include <gkyl_dg_lbo_vlasov_drag_priv.h>

void
test_lbo_drag_1x1v()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  double lower[] = {0.0, -3.0};
  double upper[] = {1.0,  3.0};
  int cells[] = {8, 16};
  struct gkyl_rect_grid pgrid;
  gkyl_rect_grid_init(&pgrid, pdim, lower, upper, cells);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 8 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_lbo_vlasov_drag_new(&cbasis, &pbasis, &crange, &pgrid, false);

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != 0 );
  TEST_CHECK( eqn->on_dev == eqn );

  struct dg_lbo_vlasov_drag *drag = container_of(eqn, struct dg_lbo_vlasov_drag, eqn);
  TEST_CHECK( drag->cdim == cdim );
  TEST_CHECK( drag->vdim == vdim );
  TEST_CHECK( drag->pdim == pdim );
  TEST_CHECK( drag->num_cbasis == cbasis.num_basis );
  TEST_CHECK( drag->conf_range.volume == 8 );
  // viMax pulled from upper velocity bounds of the grid
  TEST_CHECK( drag->viMax[0] == 3.0 );
  // vMaxSq is max over vdim of (upper)^2
  TEST_CHECK( drag->vMaxSq == 9.0 );
  // aux fields start NULL
  TEST_CHECK( drag->auxfields.nuSum == 0 );
  TEST_CHECK( drag->auxfields.nuPrimMomsSum == 0 );

  gkyl_dg_eqn_release(eqn);
}

void
test_lbo_drag_1x2v()
{
  int cdim = 1, vdim = 2, pdim = cdim+vdim;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);

  double lower[] = {0.0, -2.0, -4.0};
  double upper[] = {1.0,  2.0,  4.0};
  int cells[] = {4, 8, 8};
  struct gkyl_rect_grid pgrid;
  gkyl_rect_grid_init(&pgrid, pdim, lower, upper, cells);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 4 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_lbo_vlasov_drag_new(&cbasis, &pbasis, &crange, &pgrid, false);

  TEST_CHECK( eqn->num_equations == 1 );

  struct dg_lbo_vlasov_drag *drag = container_of(eqn, struct dg_lbo_vlasov_drag, eqn);
  TEST_CHECK( drag->cdim == cdim );
  TEST_CHECK( drag->vdim == vdim );
  TEST_CHECK( drag->pdim == pdim );
  TEST_CHECK( drag->viMax[0] == 2.0 );
  TEST_CHECK( drag->viMax[1] == 4.0 );
  TEST_CHECK( drag->vMaxSq == 16.0 );

  gkyl_dg_eqn_release(eqn);
}

void
test_lbo_drag_2x2v()
{
  int cdim = 2, vdim = 2, pdim = cdim+vdim;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  gkyl_cart_modal_hybrid(&pbasis, cdim, vdim);

  double lower[] = {0.0, 0.0, -5.0, -6.0};
  double upper[] = {1.0, 1.0,  5.0,  6.0};
  int cells[] = {4, 4, 8, 8};
  struct gkyl_rect_grid pgrid;
  gkyl_rect_grid_init(&pgrid, pdim, lower, upper, cells);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, cdim, (int[]) { 4, 4 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_lbo_vlasov_drag_new(&cbasis, &pbasis, &crange, &pgrid, false);

  struct dg_lbo_vlasov_drag *drag = container_of(eqn, struct dg_lbo_vlasov_drag, eqn);
  TEST_CHECK( drag->cdim == cdim );
  TEST_CHECK( drag->vdim == vdim );
  TEST_CHECK( drag->pdim == pdim );
  TEST_CHECK( drag->conf_range.volume == 16 );
  TEST_CHECK( drag->viMax[0] == 5.0 );
  TEST_CHECK( drag->viMax[1] == 6.0 );
  TEST_CHECK( drag->vMaxSq == 36.0 );

  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "lbo_drag_1x1v", test_lbo_drag_1x1v },
  { "lbo_drag_1x2v", test_lbo_drag_1x2v },
  { "lbo_drag_2x2v", test_lbo_drag_2x2v },
  { NULL, NULL },
};
