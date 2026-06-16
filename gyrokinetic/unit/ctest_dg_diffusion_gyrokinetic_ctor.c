// Test construction of the gyrokinetic diffusion DG equation object.
// Verifies that gkyl_dg_diffusion_gyrokinetic_new returns an object with the
// expected base-class field (num_equations), correctly recorded const-coeff
// flag, per-direction diffusion flags and number of phase-space basis
// functions, and that reference counting via acquire/release works.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_dg_diffusion_gyrokinetic.h>
#include <gkyl_dg_diffusion_gyrokinetic_priv.h>
#include <gkyl_range.h>

static void
check_diffusion(int cdim, int vdim, int poly_order, bool const_coeff, const bool *diff_in_dir)
{
  int pdim = cdim + vdim;

  struct gkyl_basis basis, cbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&basis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&basis, pdim, poly_order);

  // Conf-space range used to index the diffusion coefficient.
  struct gkyl_range diff_range;
  int lower[GKYL_MAX_CDIM], upper[GKYL_MAX_CDIM];
  for (int d=0; d<cdim; d++) { lower[d] = 0; upper[d] = 3; }
  gkyl_range_init(&diff_range, cdim, lower, upper);

  int diff_order = 2;
  struct gkyl_dg_eqn *eqn = gkyl_dg_diffusion_gyrokinetic_new(&basis, &cbasis,
    const_coeff, diff_in_dir, diff_order, &diff_range, false);

  TEST_CHECK( eqn != NULL );
  TEST_CHECK( eqn->num_equations == 1 );

  // Inspect the derived object (test-only access).
  struct dg_diffusion_gyrokinetic *diffusion =
    container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);

  TEST_CHECK( diffusion->const_coeff == const_coeff );
  TEST_CHECK( diffusion->num_basis == basis.num_basis );
  for (int d=0; d<cdim; d++)
    TEST_CHECK( diffusion->diff_in_dir[d] == diff_in_dir[d] );

  // Volume and surface kernels must be wired up.
  TEST_CHECK( eqn->vol_term != NULL );
  TEST_CHECK( eqn->surf_term != NULL );
  TEST_CHECK( eqn->boundary_surf_term != NULL );
  TEST_CHECK( diffusion->surf[0] != NULL );

  gkyl_dg_eqn_release(eqn);
}

void test_diff_1x1v_p1_const()
{
  bool dir[] = { true };
  check_diffusion(1, 1, 1, true, dir);
}

void test_diff_1x2v_p1_var()
{
  bool dir[] = { true };
  check_diffusion(1, 2, 1, false, dir);
}

void test_diff_2x2v_p1_xy()
{
  bool dir[] = { true, true };
  check_diffusion(2, 2, 1, true, dir);
}

void test_diff_2x2v_p1_x_only()
{
  bool dir[] = { true, false };
  check_diffusion(2, 2, 1, true, dir);
}

void test_diff_3x2v_p1()
{
  bool dir[] = { true, true, true };
  check_diffusion(3, 2, 1, true, dir);
}

void test_diff_1x1v_p2_const()
{
  bool dir[] = { true };
  check_diffusion(1, 1, 2, true, dir);
}

void
test_diff_acquire()
{
  struct gkyl_basis basis, cbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);
  gkyl_cart_modal_gkhybrid(&basis, 1, 1);

  struct gkyl_range diff_range;
  gkyl_range_init(&diff_range, 1, (int[]) { 0 }, (int[]) { 3 });

  bool dir[] = { true };
  struct gkyl_dg_eqn *eqn = gkyl_dg_diffusion_gyrokinetic_new(&basis, &cbasis,
    true, dir, 2, &diff_range, false);

  struct gkyl_dg_eqn *eqn2 = gkyl_dg_eqn_acquire(eqn);
  TEST_CHECK( eqn2 == eqn );
  TEST_CHECK( eqn2->num_equations == 1 );

  // Drop first ref; object must survive via the second reference.
  gkyl_dg_eqn_release(eqn);
  TEST_CHECK( eqn2->num_equations == 1 );

  gkyl_dg_eqn_release(eqn2);
}

TEST_LIST = {
  { "diff_1x1v_p1_const", test_diff_1x1v_p1_const },
  { "diff_1x2v_p1_var", test_diff_1x2v_p1_var },
  { "diff_2x2v_p1_xy", test_diff_2x2v_p1_xy },
  { "diff_2x2v_p1_x_only", test_diff_2x2v_p1_x_only },
  { "diff_3x2v_p1", test_diff_3x2v_p1 },
  { "diff_1x1v_p2_const", test_diff_1x1v_p2_const },
  { "diff_acquire", test_diff_acquire },
  { NULL, NULL },
};
