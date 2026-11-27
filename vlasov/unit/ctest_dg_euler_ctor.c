// Tests for the Euler (fluid) DG equation object constructor.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_euler.h>
#include <gkyl_dg_euler_priv.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wave_geom.h>

static struct gkyl_wave_geom *
mk_geom(
  struct gkyl_rect_grid *grid, struct gkyl_range *range, int ndim, const double *lower,
  const double *upper, const int *cells
)
{
  gkyl_rect_grid_init(grid, ndim, lower, upper, cells);
  gkyl_range_init_from_shape(range, ndim, cells);
  // NULL mapc2p -> identity (nomapc2p)
  return gkyl_wave_geom_new(grid, range, NULL, NULL, false);
}

// An Euler equation object on the grid with the identity position map and its own u and p.
struct euler_eqn_test {
  struct gkyl_vlasov_position_map *pos_map;
  struct gkyl_array *u, *p;
  struct gkyl_dg_eqn *eqn;
};

static struct euler_eqn_test
mk_eqn(
  const struct gkyl_basis *cbasis, const struct gkyl_rect_grid *grid,
  const struct gkyl_range *crange, const struct gkyl_wv_eqn *wv, const struct gkyl_wave_geom *geom
)
{
  struct euler_eqn_test et;
  struct gkyl_vlasov_position_map_inp inp_pmap[GKYL_MAX_CDIM] = {0};
  et.pos_map = gkyl_vlasov_position_map_new(grid, crange, crange, cbasis, inp_pmap, false);
  et.u = gkyl_array_new(GKYL_DOUBLE, 3 * cbasis->num_basis, crange->volume);
  et.p = gkyl_array_new(GKYL_DOUBLE, cbasis->num_basis, crange->volume);
  et.eqn = gkyl_dg_euler_inew(&(struct gkyl_dg_euler_inp){
    .cbasis = cbasis,
    .conf_range = crange,
    .wv_eqn = wv,
    .geom = geom,
    .pos_map = et.pos_map,
    .u = et.u,
    .p = et.p,
    .use_gpu = false,
  });
  return et;
}

static void
release_eqn(struct euler_eqn_test *et)
{
  gkyl_dg_eqn_release(et->eqn);
  gkyl_array_release(et->u);
  gkyl_array_release(et->p);
  gkyl_vlasov_position_map_release(et->pos_map);
}

void
test_euler_1x_p1()
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);

  struct gkyl_rect_grid grid;
  struct gkyl_range crange;
  struct gkyl_wave_geom *geom =
    mk_geom(&grid, &crange, 1, (double[]){0.0}, (double[]){1.0}, (int[]){16});

  double gas_gamma = 1.4;
  struct gkyl_wv_eqn *wv = gkyl_wv_euler_new(gas_gamma, false);

  struct euler_eqn_test et = mk_eqn(&cbasis, &grid, &crange, wv, geom);
  struct gkyl_dg_eqn *eqn = et.eqn;

  // Euler has 5 conserved variables (rho, momentum x3, energy)
  TEST_CHECK(eqn->num_equations == 5);
  TEST_CHECK(eqn->num_equations == wv->num_equations);

  TEST_CHECK(eqn->vol_term != 0);
  TEST_CHECK(eqn->surf_term != 0);
  TEST_CHECK(eqn->boundary_surf_term != 0);
  TEST_CHECK(eqn->on_dev == eqn);

  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  TEST_CHECK(euler->conf_range.volume == 16);
  TEST_CHECK(euler->conf_range.ndim == 1);
  TEST_CHECK(gkyl_compare(euler->gas_gamma, gas_gamma, 1e-15));
  TEST_CHECK(euler->surf[0] != 0);
  TEST_CHECK(euler->wv_eqn == wv);
  TEST_CHECK(euler->geom == geom);

  // the arrays of the input and the kernels of the volume term and of the absorbing wall
  TEST_CHECK(euler->u == et.u);
  TEST_CHECK(euler->p == et.p);
  TEST_CHECK(euler->wall_flux[0] == euler_wall_fluxx_1x_ser_p1);
  TEST_CHECK(euler->boundary_surf[0] == euler_boundary_surfx_1x_ser_p1);
  TEST_CHECK(euler->vol == euler_vol_1x_ser_p1);

  release_eqn(&et);
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
  struct gkyl_wave_geom *geom =
    mk_geom(&grid, &crange, 2, (double[]){0.0, 0.0}, (double[]){1.0, 1.0}, (int[]){8, 4});

  double gas_gamma = 5.0 / 3.0;
  struct gkyl_wv_eqn *wv = gkyl_wv_euler_new(gas_gamma, false);

  struct euler_eqn_test et = mk_eqn(&cbasis, &grid, &crange, wv, geom);
  struct gkyl_dg_eqn *eqn = et.eqn;

  TEST_CHECK(eqn->num_equations == 5);
  TEST_CHECK(eqn->vol_term != 0);

  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  TEST_CHECK(euler->conf_range.volume == 32);
  TEST_CHECK(euler->conf_range.ndim == 2);
  TEST_CHECK(gkyl_compare(euler->gas_gamma, gas_gamma, 1e-15));
  TEST_CHECK(euler->surf[0] != 0);
  TEST_CHECK(euler->surf[1] != 0);

  release_eqn(&et);
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
  struct gkyl_wave_geom *geom = mk_geom(
    &grid, &crange, 3, (double[]){0.0, 0.0, 0.0}, (double[]){1.0, 1.0, 1.0}, (int[]){4, 3, 2}
  );

  struct gkyl_wv_eqn *wv = gkyl_wv_euler_new(1.4, false);

  struct euler_eqn_test et = mk_eqn(&cbasis, &grid, &crange, wv, geom);
  struct gkyl_dg_eqn *eqn = et.eqn;

  TEST_CHECK(eqn->num_equations == 5);

  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  TEST_CHECK(euler->conf_range.volume == 24);
  TEST_CHECK(euler->conf_range.ndim == 3);
  TEST_CHECK(euler->surf[0] != 0);
  TEST_CHECK(euler->surf[1] != 0);
  TEST_CHECK(euler->surf[2] != 0);

  release_eqn(&et);
  gkyl_wv_eqn_release(wv);
  gkyl_wave_geom_release(geom);
}

// Isothermal Euler: four conserved variables (rho, momentum x3) and the kernels of every direction.
static void
test_iso_euler(int cdim)
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, 1);
  struct gkyl_rect_grid grid;
  struct gkyl_range crange;
  struct gkyl_wave_geom *geom = mk_geom(
    &grid, &crange, cdim, (double[]){0.0, 0.0, 0.0}, (double[]){1.0, 1.0, 1.0}, (int[]){4, 3, 2}
  );

  double vt = 0.7;
  struct gkyl_wv_eqn *wv = gkyl_wv_iso_euler_new(vt, false);
  struct euler_eqn_test et = mk_eqn(&cbasis, &grid, &crange, wv, geom);
  struct gkyl_dg_eqn *eqn = et.eqn;

  TEST_CHECK(eqn->num_equations == 4);
  TEST_CHECK(eqn->num_equations == wv->num_equations);
  TEST_CHECK(eqn->vol_term != 0);
  TEST_CHECK(eqn->surf_term != 0);
  TEST_CHECK(eqn->on_dev == eqn);

  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  TEST_CHECK(euler->eqn_type == GKYL_EQN_ISO_EULER);
  TEST_CHECK(gkyl_compare(euler->vt, vt, 1e-15));
  for (int d = 0; d < cdim; ++d) {
    TEST_CHECK(euler->surf[d] != 0);
  }
  TEST_CHECK(euler->wv_eqn == wv);

  release_eqn(&et);
  gkyl_wv_eqn_release(wv);
  gkyl_wave_geom_release(geom);
}

void
test_iso_euler_1x_p1()
{
  test_iso_euler(1);
}

void
test_iso_euler_2x_p1()
{
  test_iso_euler(2);
}

void
test_iso_euler_3x_p1()
{
  test_iso_euler(3);
}

TEST_LIST = {
  {"euler_1x_p1", test_euler_1x_p1},
  {"euler_2x_p1", test_euler_2x_p1},
  {"euler_3x_p1", test_euler_3x_p1},
  {"iso_euler_1x_p1", test_iso_euler_1x_p1},
  {"iso_euler_2x_p1", test_iso_euler_2x_p1},
  {"iso_euler_3x_p1", test_iso_euler_3x_p1},
  {NULL, NULL}
};
