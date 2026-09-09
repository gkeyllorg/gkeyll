// Test construction of the PKPM Euler (fluid) DG equation object.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_euler_pkpm.h>
#include <gkyl_range.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_rect_grid.h>

static void
mk_conf_range(int cdim, struct gkyl_range *local, struct gkyl_range *local_ext)
{
  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  int cells[GKYL_MAX_DIM], ghost[GKYL_MAX_DIM];
  for (int d=0; d<cdim; ++d) { lower[d] = -1.0; upper[d] = 1.0; cells[d] = 6; ghost[d] = 1; }
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, cdim, lower, upper, cells);
  gkyl_create_grid_ranges(&grid, ghost, local_ext, local);
}

void
test_euler_pkpm_1x_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);

  struct gkyl_range local, local_ext;
  mk_conf_range(1, &local, &local_ext);

  struct gkyl_dg_eqn *eqn = gkyl_dg_euler_pkpm_new(&cbasis, &local, false);

  // PKPM Euler solves the 3 momentum components (rho ux, rho uy, rho uz).
  TEST_CHECK( eqn->num_equations == 3 );
  TEST_CHECK( eqn->vol_term != NULL );
  TEST_CHECK( eqn->surf_term != NULL );

  gkyl_dg_eqn_release(eqn);
}

void
test_euler_pkpm_1x_p2()
{
  int poly_order = 2;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);

  struct gkyl_range local, local_ext;
  mk_conf_range(1, &local, &local_ext);

  struct gkyl_dg_eqn *eqn = gkyl_dg_euler_pkpm_new(&cbasis, &local, false);
  TEST_CHECK( eqn->num_equations == 3 );

  struct gkyl_dg_eqn *eqn2 = gkyl_dg_eqn_acquire(eqn);
  TEST_CHECK( eqn2 == eqn );
  gkyl_dg_eqn_release(eqn2);
  gkyl_dg_eqn_release(eqn);
}

void
test_euler_pkpm_2x_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, poly_order);

  struct gkyl_range local, local_ext;
  mk_conf_range(2, &local, &local_ext);

  struct gkyl_dg_eqn *eqn = gkyl_dg_euler_pkpm_new(&cbasis, &local, false);
  TEST_CHECK( eqn->num_equations == 3 );
  TEST_CHECK( eqn->surf_term != NULL );
  gkyl_dg_eqn_release(eqn);
}

void
test_euler_pkpm_3x_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, poly_order);

  struct gkyl_range local, local_ext;
  mk_conf_range(3, &local, &local_ext);

  struct gkyl_dg_eqn *eqn = gkyl_dg_euler_pkpm_new(&cbasis, &local, false);
  TEST_CHECK( eqn->num_equations == 3 );
  TEST_CHECK( eqn->vol_term != NULL );
  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "euler_pkpm_1x_p1", test_euler_pkpm_1x_p1 },
  { "euler_pkpm_1x_p2", test_euler_pkpm_1x_p2 },
  { "euler_pkpm_2x_p1", test_euler_pkpm_2x_p1 },
  { "euler_pkpm_3x_p1", test_euler_pkpm_3x_p1 },
  { NULL, NULL },
};
