// Tests for the fluid diffusion DG equation object constructor.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_dg_diffusion_fluid.h>
#include <gkyl_dg_diffusion_fluid_priv.h>

void
test_diffusion_fluid_1x_scalar()
{
  int cdim = 1;
  int poly_order = 2;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  bool diff_in_dir[] = { true };
  int diff_order = 2;
  int num_eq = 1;

  struct gkyl_range drange;
  gkyl_range_init_from_shape(&drange, cdim, (int[]) { 10 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_diffusion_fluid_new(&basis, true, num_eq,
    diff_in_dir, diff_order, &drange, false);

  TEST_CHECK( eqn->num_equations == num_eq );
  TEST_CHECK( eqn->vol_term != 0 );
  TEST_CHECK( eqn->surf_term != 0 );
  TEST_CHECK( eqn->boundary_surf_term != 0 );
  TEST_CHECK( eqn->on_dev == eqn );

  struct dg_diffusion_fluid *diff = container_of(eqn, struct dg_diffusion_fluid, eqn);
  TEST_CHECK( diff->const_coeff == true );
  TEST_CHECK( diff->num_equations == num_eq );
  TEST_CHECK( diff->num_basis == basis.num_basis );
  TEST_CHECK( diff->diff_in_dir[0] == true );
  TEST_CHECK( diff->diff_range.volume == 10 );
  TEST_CHECK( diff->auxfields.D == 0 );

  gkyl_dg_eqn_release(eqn);
}

void
test_diffusion_fluid_2x_system()
{
  int cdim = 2;
  int poly_order = 1;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  bool diff_in_dir[] = { true, true };
  int diff_order = 2;
  int num_eq = 5; // e.g. 5-moment fluid system

  struct gkyl_range drange;
  gkyl_range_init_from_shape(&drange, cdim, (int[]) { 8, 6 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_diffusion_fluid_new(&basis, true, num_eq,
    diff_in_dir, diff_order, &drange, false);

  TEST_CHECK( eqn->num_equations == num_eq );

  struct dg_diffusion_fluid *diff = container_of(eqn, struct dg_diffusion_fluid, eqn);
  TEST_CHECK( diff->num_equations == 5 );
  TEST_CHECK( diff->diff_in_dir[0] == true );
  TEST_CHECK( diff->diff_in_dir[1] == true );
  TEST_CHECK( diff->diff_range.volume == 48 );

  gkyl_dg_eqn_release(eqn);
}

void
test_diffusion_fluid_1x_varcoeff()
{
  int cdim = 1;
  int poly_order = 2;
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, cdim, poly_order);

  bool diff_in_dir[] = { true };
  int diff_order = 2;
  int num_eq = 1;

  struct gkyl_range drange;
  gkyl_range_init_from_shape(&drange, cdim, (int[]) { 7 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_diffusion_fluid_new(&basis, false, num_eq,
    diff_in_dir, diff_order, &drange, false);

  struct dg_diffusion_fluid *diff = container_of(eqn, struct dg_diffusion_fluid, eqn);
  TEST_CHECK( diff->const_coeff == false );
  TEST_CHECK( diff->diff_range.volume == 7 );

  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "diffusion_fluid_1x_scalar", test_diffusion_fluid_1x_scalar },
  { "diffusion_fluid_2x_system", test_diffusion_fluid_2x_system },
  { "diffusion_fluid_1x_varcoeff", test_diffusion_fluid_1x_varcoeff },
  { NULL, NULL },
};
