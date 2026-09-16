#include <acutest.h>
#include <math.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_wv_canonical_pb_fluid.h>

// Incompressible Euler: scalar (vorticity) equation, type tag and equation count.
void
test_can_pb_incompress_euler_basic()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_can_pb_incompress_euler_new();

  TEST_CHECK( eqn->type == GKYL_EQN_CAN_PB_INCOMPRESS_EULER );
  TEST_CHECK( eqn->num_equations == 1 );

  // on_dev points back to the object itself for the CPU build.
  TEST_CHECK( eqn->on_dev == eqn );

  gkyl_wv_eqn_release(eqn);
}

// Reference-count round-trip: acquire bumps and release returns the same object.
void
test_can_pb_incompress_euler_refcount()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_can_pb_incompress_euler_new();

  struct gkyl_wv_eqn *eqn2 = gkyl_wv_eqn_acquire(eqn);
  TEST_CHECK( eqn2 == eqn );

  // Release the extra reference; object stays alive.
  gkyl_wv_eqn_release(eqn2);
  TEST_CHECK( eqn->num_equations == 1 );

  gkyl_wv_eqn_release(eqn);
}

// Hasegawa-Mima: scalar equation with distinct type tag.
void
test_can_pb_hasegawa_mima_basic()
{
  struct gkyl_wv_eqn *eqn = gkyl_wv_can_pb_hasegawa_mima_new();

  TEST_CHECK( eqn->type == GKYL_EQN_CAN_PB_HASEGAWA_MIMA );
  TEST_CHECK( eqn->num_equations == 1 );
  TEST_CHECK( eqn->on_dev == eqn );

  gkyl_wv_eqn_release(eqn);
}

// Hasegawa-Wakatani: two-component system (vorticity + density); alpha and
// is_modified accessors return exactly what was passed at construction.
void
test_can_pb_hasegawa_wakatani_basic()
{
  double alpha = 2.5;
  struct gkyl_wv_eqn *eqn = gkyl_wv_can_pb_hasegawa_wakatani_new(alpha, false);

  TEST_CHECK( eqn->type == GKYL_EQN_CAN_PB_HASEGAWA_WAKATANI );
  TEST_CHECK( eqn->num_equations == 2 );
  TEST_CHECK( eqn->on_dev == eqn );

  TEST_CHECK( gkyl_compare(gkyl_wv_can_pb_hasegawa_wakatani_alpha(eqn), alpha, 1e-15) );
  TEST_CHECK( gkyl_wv_can_pb_hasegawa_wakatani_is_modified(eqn) == false );

  gkyl_wv_eqn_release(eqn);
}

// Modified Hasegawa-Wakatani with a different adiabaticity parameter.
void
test_can_pb_hasegawa_wakatani_modified()
{
  double alpha = 0.375;
  struct gkyl_wv_eqn *eqn = gkyl_wv_can_pb_hasegawa_wakatani_new(alpha, true);

  TEST_CHECK( eqn->num_equations == 2 );
  TEST_CHECK( gkyl_compare(gkyl_wv_can_pb_hasegawa_wakatani_alpha(eqn), alpha, 1e-15) );
  TEST_CHECK( gkyl_wv_can_pb_hasegawa_wakatani_is_modified(eqn) == true );

  gkyl_wv_eqn_release(eqn);
}

// Distinct Hasegawa-Wakatani objects keep independent parameters.
void
test_can_pb_hasegawa_wakatani_independent()
{
  struct gkyl_wv_eqn *e1 = gkyl_wv_can_pb_hasegawa_wakatani_new(1.0, true);
  struct gkyl_wv_eqn *e2 = gkyl_wv_can_pb_hasegawa_wakatani_new(7.0, false);

  TEST_CHECK( e1 != e2 );
  TEST_CHECK( gkyl_compare(gkyl_wv_can_pb_hasegawa_wakatani_alpha(e1), 1.0, 1e-15) );
  TEST_CHECK( gkyl_compare(gkyl_wv_can_pb_hasegawa_wakatani_alpha(e2), 7.0, 1e-15) );
  TEST_CHECK( gkyl_wv_can_pb_hasegawa_wakatani_is_modified(e1) == true );
  TEST_CHECK( gkyl_wv_can_pb_hasegawa_wakatani_is_modified(e2) == false );

  gkyl_wv_eqn_release(e1);
  gkyl_wv_eqn_release(e2);
}

TEST_LIST = {
  { "can_pb_incompress_euler_basic", test_can_pb_incompress_euler_basic },
  { "can_pb_incompress_euler_refcount", test_can_pb_incompress_euler_refcount },
  { "can_pb_hasegawa_mima_basic", test_can_pb_hasegawa_mima_basic },
  { "can_pb_hasegawa_wakatani_basic", test_can_pb_hasegawa_wakatani_basic },
  { "can_pb_hasegawa_wakatani_modified", test_can_pb_hasegawa_wakatani_modified },
  { "can_pb_hasegawa_wakatani_independent", test_can_pb_hasegawa_wakatani_independent },
  { NULL, NULL },
};
