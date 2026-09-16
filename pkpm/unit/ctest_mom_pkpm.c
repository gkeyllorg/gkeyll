// Test construction of PKPM moment-type objects.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_mom_pkpm.h>
#include <gkyl_mom_type.h>

void
test_mom_pkpm_1x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order); // 1X
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order); // 1X1V

  // Non-diagnostic moments (rho, p_par, p_perp, M1).
  struct gkyl_mom_type *m = gkyl_mom_pkpm_new(&cbasis, &pbasis, 1.0, false, false);
  TEST_CHECK( m->cdim == 1 );
  TEST_CHECK( m->pdim == 2 );
  TEST_CHECK( m->poly_order == 1 );
  TEST_CHECK( m->num_config == cbasis.num_basis );
  TEST_CHECK( m->num_phase == pbasis.num_basis );
  TEST_CHECK( m->num_mom == 4 );
  TEST_CHECK( m->kernel != NULL );

  // Diagnostic moments has 8 components.
  struct gkyl_mom_type *md = gkyl_mom_pkpm_new(&cbasis, &pbasis, 2.0, true, false);
  TEST_CHECK( md->num_mom == 8 );
  TEST_CHECK( md->cdim == 1 );
  TEST_CHECK( md->kernel != NULL );

  gkyl_mom_type_release(m);
  gkyl_mom_type_release(md);
}

void
test_mom_pkpm_1x1v_p2()
{
  int poly_order = 2;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, poly_order);
  gkyl_cart_modal_serendip(&pbasis, 2, poly_order);

  struct gkyl_mom_type *m = gkyl_mom_pkpm_new(&cbasis, &pbasis, 1.0, false, false);
  TEST_CHECK( m->poly_order == 2 );
  TEST_CHECK( m->num_mom == 4 );
  TEST_CHECK( m->num_config == cbasis.num_basis );
  gkyl_mom_type_release(m);
}

void
test_mom_pkpm_2x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 2, poly_order); // 2X
  gkyl_cart_modal_serendip(&pbasis, 3, poly_order); // 2X1V

  struct gkyl_mom_type *m = gkyl_mom_pkpm_new(&cbasis, &pbasis, 1.0, false, false);
  TEST_CHECK( m->cdim == 2 );
  TEST_CHECK( m->pdim == 3 );
  TEST_CHECK( m->num_mom == 4 );
  struct gkyl_mom_type *md = gkyl_mom_pkpm_new(&cbasis, &pbasis, 1.0, true, false);
  TEST_CHECK( md->num_mom == 8 );
  gkyl_mom_type_release(m);
  gkyl_mom_type_release(md);
}

void
test_mom_pkpm_3x1v_p1()
{
  int poly_order = 1;
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 3, poly_order); // 3X
  gkyl_cart_modal_serendip(&pbasis, 4, poly_order); // 3X1V

  struct gkyl_mom_type *m = gkyl_mom_pkpm_new(&cbasis, &pbasis, 1.0, false, false);
  TEST_CHECK( m->cdim == 3 );
  TEST_CHECK( m->pdim == 4 );
  TEST_CHECK( m->num_mom == 4 );
  TEST_CHECK( m->kernel != NULL );
  gkyl_mom_type_release(m);
}

TEST_LIST = {
  { "mom_pkpm_1x1v_p1", test_mom_pkpm_1x1v_p1 },
  { "mom_pkpm_1x1v_p2", test_mom_pkpm_1x1v_p2 },
  { "mom_pkpm_2x1v_p1", test_mom_pkpm_2x1v_p1 },
  { "mom_pkpm_3x1v_p1", test_mom_pkpm_3x1v_p1 },
  { NULL, NULL },
};
