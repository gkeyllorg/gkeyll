// Tests for the Euler (fluid) DG equation object constructor.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_euler.h>
#include <gkyl_dg_euler_priv.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wave_geom.h>

static struct gkyl_wave_geom*
mk_geom(struct gkyl_rect_grid *grid, struct gkyl_range *range, int ndim,
  const double *lower, const double *upper, const int *cells)
{
  gkyl_rect_grid_init(grid, ndim, lower, upper, cells);
  gkyl_range_init_from_shape(range, ndim, cells);
  // NULL mapc2p -> identity (nomapc2p)
  return gkyl_wave_geom_new(grid, range, NULL, NULL, false);
}

void
test_euler_1x_p1()
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);

  struct gkyl_rect_grid grid;
  struct gkyl_range crange;
  struct gkyl_wave_geom *geom = mk_geom(&grid, &crange, 1,
    (double[]) { 0.0 }, (double[]) { 1.0 }, (int[]) { 16 });

  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *wv = gkyl_wv_euler_new(gas_gamma, false);

  struct gkyl_dg_eqn *eqn = gkyl_dg_euler_new(&cbasis, &crange, wv, geom, false);

  // Euler has 5 conserved variables (rho, momentum x3, energy)
  TEST_CHECK( eqn->num_equations == 5 );
  TEST_CHECK( eqn->num_equations == wv->num_equations );

  TEST_CHECK( eqn->vol_term != 0 );
  TEST_CHECK( eqn->surf_term != 0 );
  TEST_CHECK( eqn->boundary_surf_term != 0 );
  TEST_CHECK( eqn->on_dev == eqn );

  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  TEST_CHECK( euler->conf_range.volume == 16 );
  TEST_CHECK( euler->conf_range.ndim == 1 );
  TEST_CHECK( gkyl_compare(euler->gas_gamma, gas_gamma, 1e-15) );
  TEST_CHECK( euler->surf[0] != 0 );
  TEST_CHECK( euler->wv_eqn == wv );
  TEST_CHECK( euler->geom == geom );

  // aux fields start NULL
  TEST_CHECK( euler->auxfields.u == 0 );
  TEST_CHECK( euler->auxfields.p == 0 );

  gkyl_dg_eqn_release(eqn);
  gkyl_wv_eqn_release(wv);
  gkyl_wave_geom_release(geom);
}

void
test_euler_2x_p1()
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, 1);

  struct gkyl_rect_grid grid;
  struct gkyl_range crange;
  struct gkyl_wave_geom *geom = mk_geom(&grid, &crange, 2,
    (double[]) { 0.0, 0.0 }, (double[]) { 1.0, 1.0 }, (int[]) { 8, 4 });

  double gas_gamma = 5.0/3.0;
  struct gkyl_wv_eqn *wv = gkyl_wv_euler_new(gas_gamma, false);

  struct gkyl_dg_eqn *eqn = gkyl_dg_euler_new(&cbasis, &crange, wv, geom, false);

  TEST_CHECK( eqn->num_equations == 5 );
  TEST_CHECK( eqn->vol_term != 0 );

  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  TEST_CHECK( euler->conf_range.volume == 32 );
  TEST_CHECK( euler->conf_range.ndim == 2 );
  TEST_CHECK( gkyl_compare(euler->gas_gamma, gas_gamma, 1e-15) );
  TEST_CHECK( euler->surf[0] != 0 );
  TEST_CHECK( euler->surf[1] != 0 );

  gkyl_dg_eqn_release(eqn);
  gkyl_wv_eqn_release(wv);
  gkyl_wave_geom_release(geom);
}

void
test_euler_3x_p1()
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, 1);

  struct gkyl_rect_grid grid;
  struct gkyl_range crange;
  struct gkyl_wave_geom *geom = mk_geom(&grid, &crange, 3,
    (double[]) { 0.0, 0.0, 0.0 }, (double[]) { 1.0, 1.0, 1.0 }, (int[]) { 4, 3, 2 });

  struct gkyl_wv_eqn *wv = gkyl_wv_euler_new(1.4, false);

  struct gkyl_dg_eqn *eqn = gkyl_dg_euler_new(&cbasis, &crange, wv, geom, false);

  TEST_CHECK( eqn->num_equations == 5 );

  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  TEST_CHECK( euler->conf_range.volume == 24 );
  TEST_CHECK( euler->conf_range.ndim == 3 );
  TEST_CHECK( euler->surf[0] != 0 );
  TEST_CHECK( euler->surf[1] != 0 );
  TEST_CHECK( euler->surf[2] != 0 );

  gkyl_dg_eqn_release(eqn);
  gkyl_wv_eqn_release(wv);
  gkyl_wave_geom_release(geom);
}

TEST_LIST = {
  { "euler_1x_p1", test_euler_1x_p1 },
  { "euler_2x_p1", test_euler_2x_p1 },
  { "euler_3x_p1", test_euler_3x_p1 },
  { NULL, NULL },
};
