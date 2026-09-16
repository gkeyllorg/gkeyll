// Test construction of the PKPM Vlasov DG equation object.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_vlasov_pkpm.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static void
mk_ranges(int cdim, struct gkyl_range *conf, struct gkyl_range *conf_ext,
  struct gkyl_range *phase, struct gkyl_range *phase_ext)
{
  int pdim = cdim+1;
  double clower[GKYL_MAX_DIM], cupper[GKYL_MAX_DIM];
  int ccells[GKYL_MAX_DIM], cghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) { clower[d] = -1.0; cupper[d] = 1.0; ccells[d] = 4; cghost[d] = 1; }
  struct gkyl_rect_grid cgrid;
  gkyl_rect_grid_init(&cgrid, cdim, clower, cupper, ccells);
  gkyl_create_grid_ranges(&cgrid, cghost, conf_ext, conf);

  double plower[GKYL_MAX_DIM], pupper[GKYL_MAX_DIM];
  int pcells[GKYL_MAX_DIM], pghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) { plower[d] = -1.0; pupper[d] = 1.0; pcells[d] = 4; pghost[d] = 1; }
  plower[cdim] = -6.0; pupper[cdim] = 6.0; pcells[cdim] = 8; pghost[cdim] = 0;
  struct gkyl_rect_grid pgrid;
  gkyl_rect_grid_init(&pgrid, pdim, plower, pupper, pcells);
  gkyl_create_grid_ranges(&pgrid, pghost, phase_ext, phase);
}

void
test_vlasov_pkpm_1x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order);

  struct gkyl_range conf, conf_ext, phase, phase_ext;
  mk_ranges(1, &conf, &conf_ext, &phase, &phase_ext);

  struct gkyl_dg_eqn *eqn =
    gkyl_dg_vlasov_pkpm_new(&cbasis, &pbasis, &conf, &phase, false);

  // PKPM Vlasov evolves 2 distributions (F_0 and G).
  TEST_CHECK( eqn->num_equations == 2 );
  TEST_CHECK( eqn->vol_term != NULL );
  TEST_CHECK( eqn->surf_term != NULL );
  TEST_CHECK( eqn->boundary_surf_term != NULL );

  gkyl_dg_eqn_release(eqn);
}

void
test_vlasov_pkpm_1x1v_p2()
{
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order);

  struct gkyl_range conf, conf_ext, phase, phase_ext;
  mk_ranges(1, &conf, &conf_ext, &phase, &phase_ext);

  struct gkyl_dg_eqn *eqn =
    gkyl_dg_vlasov_pkpm_new(&cbasis, &pbasis, &conf, &phase, false);
  TEST_CHECK( eqn->num_equations == 2 );
  TEST_CHECK( eqn->vol_term != NULL );

  // Acquire/release reference-counting smoke test.
  struct gkyl_dg_eqn *eqn2 = gkyl_dg_eqn_acquire(eqn);
  TEST_CHECK( eqn2 == eqn );
  gkyl_dg_eqn_release(eqn2);
  gkyl_dg_eqn_release(eqn);
}

void
test_vlasov_pkpm_2x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 3, poly_order);

  struct gkyl_range conf, conf_ext, phase, phase_ext;
  mk_ranges(2, &conf, &conf_ext, &phase, &phase_ext);

  struct gkyl_dg_eqn *eqn =
    gkyl_dg_vlasov_pkpm_new(&cbasis, &pbasis, &conf, &phase, false);
  TEST_CHECK( eqn->num_equations == 2 );
  TEST_CHECK( eqn->surf_term != NULL );
  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "vlasov_pkpm_1x1v_p1", test_vlasov_pkpm_1x1v_p1 },
  { "vlasov_pkpm_1x1v_p2", test_vlasov_pkpm_1x1v_p2 },
  { "vlasov_pkpm_2x1v_p1", test_vlasov_pkpm_2x1v_p1 },
  { NULL, NULL },
};
