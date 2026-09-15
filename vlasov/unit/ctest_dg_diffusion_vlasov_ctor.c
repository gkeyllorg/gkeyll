// Tests for the Vlasov diffusion DG equation object constructor.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_dg_diffusion_vlasov.h>
#include <gkyl_dg_diffusion_vlasov_priv.h>

void
test_diffusion_vlasov_1x1v_const()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  bool diff_in_dir[] = { true };
  int diff_order = 2;

  struct gkyl_range drange;
  gkyl_range_init_from_shape(&drange, cdim, (int[]) { 12 });

  struct gkyl_dg_eqn *eqn = gkyl_dg_diffusion_vlasov_new(&pbasis, &cbasis,
    true, diff_in_dir, diff_order, &drange, false);

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != 0 );
  TEST_CHECK( eqn->surf_term != 0 );
  TEST_CHECK( eqn->boundary_surf_term != 0 );
  TEST_CHECK( eqn->on_dev == eqn );

  struct dg_diffusion_vlasov *diff = container_of(eqn, struct dg_diffusion_vlasov, eqn);
  TEST_CHECK( diff->const_coeff == true );
  TEST_CHECK( diff->num_basis == pbasis.num_basis );
  TEST_CHECK( diff->diff_in_dir[0] == true );
  TEST_CHECK( diff->diff_range.volume == 12 );
  TEST_CHECK( diff->auxfields.D == 0 );

  gkyl_dg_eqn_release(eqn);
}

void
test_diffusion_vlasov_1x1v_varcoeff()
{
  int cdim = 1, vdim = 1, pdim = cdim+vdim;
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  bool diff_in_dir[] = { true };
  int diff_order = 2;

  struct gkyl_range drange;
  gkyl_range_init_from_shape(&drange, cdim, (int[]) { 5 });

  // non-constant diffusion coefficient path
  struct gkyl_dg_eqn *eqn = gkyl_dg_diffusion_vlasov_new(&pbasis, &cbasis,
    false, diff_in_dir, diff_order, &drange, false);

  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->vol_term != 0 );

  struct dg_diffusion_vlasov *diff = container_of(eqn, struct dg_diffusion_vlasov, eqn);
  TEST_CHECK( diff->const_coeff == false );
  TEST_CHECK( diff->diff_range.volume == 5 );

  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "diffusion_vlasov_1x1v_const", test_diffusion_vlasov_1x1v_const },
  { "diffusion_vlasov_1x1v_varcoeff", test_diffusion_vlasov_1x1v_varcoeff },
  { NULL, NULL },
};
