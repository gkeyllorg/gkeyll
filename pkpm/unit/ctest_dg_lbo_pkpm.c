// Test construction of the PKPM LBO drag and diffusion DG equation objects.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_lbo_pkpm_diff.h>
#include <gkyl_dg_lbo_pkpm_drag.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static void
mk_setup(int cdim, struct gkyl_range *conf, struct gkyl_range *conf_ext,
  struct gkyl_rect_grid *pgrid)
{
  int pdim = cdim+1;
  double clower[GKYL_MAX_DIM], cupper[GKYL_MAX_DIM];
  int ccells[GKYL_MAX_DIM], cghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) { clower[d] = -1.0; cupper[d] = 1.0; ccells[d] = 4; cghost[d] = 1; }
  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, cdim, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cghost, conf_ext, conf);

  double plower[GKYL_MAX_DIM], pupper[GKYL_MAX_DIM];
  int pcells[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) { plower[d] = -1.0; pupper[d] = 1.0; pcells[d] = 4; }
  plower[cdim] = -6.0; pupper[cdim] = 6.0; pcells[cdim] = 8;
  gkyl_rect_grid_init(pgrid, pdim, plower, pupper, pcells);
}

void
test_lbo_pkpm_drag_1x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order);

  struct gkyl_range conf, conf_ext;
  struct gkyl_rect_grid pgrid;
  mk_setup(1, &conf, &conf_ext, &pgrid);

  struct gkyl_dg_eqn *eqn =
    gkyl_dg_lbo_pkpm_drag_new(&cbasis, &pbasis, &conf, &pgrid, false);

  TEST_CHECK( eqn->num_equations == 2 );
  TEST_CHECK( eqn->vol_term != NULL );
  TEST_CHECK( eqn->surf_term != NULL );
  TEST_CHECK( eqn->boundary_surf_term != NULL );

  gkyl_dg_eqn_release(eqn);
}

void
test_lbo_pkpm_diff_1x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order);

  struct gkyl_range conf, conf_ext;
  struct gkyl_rect_grid pgrid;
  mk_setup(1, &conf, &conf_ext, &pgrid);

  struct gkyl_dg_eqn *eqn =
    gkyl_dg_lbo_pkpm_diff_new(&cbasis, &pbasis, &conf, &pgrid, false);

  TEST_CHECK( eqn->num_equations == 2 );
  TEST_CHECK( eqn->vol_term != NULL );
  TEST_CHECK( eqn->surf_term != NULL );

  gkyl_dg_eqn_release(eqn);
}

void
test_lbo_pkpm_drag_1x1v_p2()
{
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order);

  struct gkyl_range conf, conf_ext;
  struct gkyl_rect_grid pgrid;
  mk_setup(1, &conf, &conf_ext, &pgrid);

  struct gkyl_dg_eqn *drag =
    gkyl_dg_lbo_pkpm_drag_new(&cbasis, &pbasis, &conf, &pgrid, false);
  struct gkyl_dg_eqn *diff =
    gkyl_dg_lbo_pkpm_diff_new(&cbasis, &pbasis, &conf, &pgrid, false);
  TEST_CHECK( drag->num_equations == 2 );
  TEST_CHECK( diff->num_equations == 2 );
  gkyl_dg_eqn_release(drag);
  gkyl_dg_eqn_release(diff);
}

void
test_lbo_pkpm_2x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 3, poly_order);

  struct gkyl_range conf, conf_ext;
  struct gkyl_rect_grid pgrid;
  mk_setup(2, &conf, &conf_ext, &pgrid);

  struct gkyl_dg_eqn *drag =
    gkyl_dg_lbo_pkpm_drag_new(&cbasis, &pbasis, &conf, &pgrid, false);
  struct gkyl_dg_eqn *diff =
    gkyl_dg_lbo_pkpm_diff_new(&cbasis, &pbasis, &conf, &pgrid, false);
  TEST_CHECK( drag->num_equations == 2 );
  TEST_CHECK( diff->num_equations == 2 );
  TEST_CHECK( drag->surf_term != NULL );
  TEST_CHECK( diff->surf_term != NULL );
  gkyl_dg_eqn_release(drag);
  gkyl_dg_eqn_release(diff);
}

TEST_LIST = {
  { "lbo_pkpm_drag_1x1v_p1", test_lbo_pkpm_drag_1x1v_p1 },
  { "lbo_pkpm_diff_1x1v_p1", test_lbo_pkpm_diff_1x1v_p1 },
  { "lbo_pkpm_drag_1x1v_p2", test_lbo_pkpm_drag_1x1v_p2 },
  { "lbo_pkpm_2x1v_p1", test_lbo_pkpm_2x1v_p1 },
  { NULL, NULL },
};
