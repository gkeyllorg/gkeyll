// Tests for the advection DG equation object constructor.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_advection_priv.h>

void
test_advection_1x()
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 2);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, 1, (int[]) { 16 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_advection_new(&cbasis, &crange, false);

  // advection is a scalar equation
  TEST_CHECK( eqn->num_equations == 1 );

  // volume/surface terms must be set
  TEST_CHECK( eqn->vol_term != 0 );
  TEST_CHECK( eqn->surf_term != 0 );
  TEST_CHECK( eqn->boundary_surf_term != 0 );

  // CPU eqn obj should point to itself
  TEST_CHECK( eqn->on_dev == eqn );

  // inspect internals (testing only)
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  TEST_CHECK( advection->conf_range.volume == 16 );
  TEST_CHECK( advection->conf_range.ndim == 1 );
  // aux fields start NULL
  TEST_CHECK( advection->auxfields.u_i == 0 );

  gkyl_dg_eqn_release(eqn);
}

void
test_advection_2x()
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, 1);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, 2, (int[]) { 8, 4 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_advection_new(&cbasis, &crange, false);

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != 0 );

  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  TEST_CHECK( advection->conf_range.volume == 32 );
  TEST_CHECK( advection->conf_range.ndim == 2 );

  gkyl_dg_eqn_release(eqn);
}

void
test_advection_3x()
{
  struct gkyl_basis cbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, 1);

  struct gkyl_range crange;
  gkyl_range_init_from_shape(&crange, 3, (int[]) { 4, 3, 2 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_advection_new(&cbasis, &crange, false);

  TEST_CHECK( eqn->num_equations == 1 );

  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  TEST_CHECK( advection->conf_range.volume == 24 );
  TEST_CHECK( advection->conf_range.ndim == 3 );

  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "advection_1x", test_advection_1x },
  { "advection_2x", test_advection_2x },
  { "advection_3x", test_advection_3x },
  { NULL, NULL },
};
