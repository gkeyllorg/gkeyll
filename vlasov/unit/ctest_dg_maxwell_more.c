// Additional tests for the Maxwell DG equation object constructor,
// verifying the error-speed-factor scaling and multiple dimensions/orders.
//
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_maxwell_priv.h>

void
test_maxwell_scaling()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 2);

  // gkyl_dg_maxwell_new(basis, lightSpeed, elcErrorSpeedFactor, mgnErrorSpeedFactor, use_gpu)
  double c = 2.0, elcFac = 1.5, mgnFac = 0.5;
  struct gkyl_dg_eqn *eqn = gkyl_dg_maxwell_new(&basis, c, elcFac, mgnFac, false);

  TEST_CHECK( eqn->num_equations == 8 );
  TEST_CHECK( eqn->vol_term != 0 );
  TEST_CHECK( eqn->surf_term != 0 );
  TEST_CHECK( eqn->on_dev == eqn );

  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  TEST_CHECK( maxwell->maxwell_data.c == c );
  // chi = c * elcErrorSpeedFactor, gamma = c * mgnErrorSpeedFactor
  TEST_CHECK( maxwell->maxwell_data.chi == c*elcFac );    // 3.0
  TEST_CHECK( maxwell->maxwell_data.gamma == c*mgnFac );  // 1.0

  gkyl_dg_eqn_release(eqn);
}

void
test_maxwell_2x()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 2, 1);

  struct gkyl_dg_eqn *eqn = gkyl_dg_maxwell_new(&basis, 1.0, 1.0, 1.0, false);

  TEST_CHECK( eqn->num_equations == 8 );

  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  TEST_CHECK( maxwell->maxwell_data.c == 1.0 );
  TEST_CHECK( maxwell->maxwell_data.chi == 1.0 );
  TEST_CHECK( maxwell->maxwell_data.gamma == 1.0 );

  gkyl_dg_eqn_release(eqn);
}

void
test_maxwell_3x()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 3, 1);

  struct gkyl_dg_eqn *eqn = gkyl_dg_maxwell_new(&basis, 3.0, 0.0, 2.0, false);

  TEST_CHECK( eqn->num_equations == 8 );

  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);
  TEST_CHECK( maxwell->maxwell_data.c == 3.0 );
  TEST_CHECK( maxwell->maxwell_data.chi == 0.0 );   // 3.0*0.0
  TEST_CHECK( maxwell->maxwell_data.gamma == 6.0 ); // 3.0*2.0

  gkyl_dg_eqn_release(eqn);
}

TEST_LIST = {
  { "maxwell_scaling", test_maxwell_scaling },
  { "maxwell_2x", test_maxwell_2x },
  { "maxwell_3x", test_maxwell_3x },
  { NULL, NULL },
};
