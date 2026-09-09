// Test construction of the gyrokinetic LBO primitive-moment type object.
// Verifies that the returned gkyl_prim_lbo_type has the expected dimension,
// basis-count and udim fields for a variety of (cdim, vdim, poly_order)
// configurations, and that acquire/release reference counting works.
#include <acutest.h>

#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_prim_lbo_gyrokinetic.h>

static void
check_prim(int cdim, int vdim, int poly_order)
{
  int pdim = cdim + vdim;

  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, cdim, poly_order);
  // Gyrokinetic phase space uses the gkhybrid basis at p=1, serendip otherwise.
  if (poly_order == 1)
    gkyl_cart_modal_gkhybrid(&pbasis, cdim, vdim);
  else
    gkyl_cart_modal_serendip(&pbasis, pdim, poly_order);

  struct gkyl_prim_lbo_type *prim = gkyl_prim_lbo_gyrokinetic_new(&cbasis, &pbasis, false);

  TEST_CHECK( prim != NULL );
  TEST_CHECK( prim->cdim == cdim );
  TEST_CHECK( prim->pdim == pdim );
  TEST_CHECK( prim->poly_order == poly_order );
  TEST_CHECK( prim->num_config == cbasis.num_basis );
  TEST_CHECK( prim->num_phase == pbasis.num_basis );
  // Gyrokinetic flow has a single (parallel) velocity component.
  TEST_CHECK( prim->udim == 1 );
  // Kernels must be wired up.
  TEST_CHECK( prim->self_prim != NULL );
  TEST_CHECK( prim->cross_prim != NULL );

  gkyl_prim_lbo_type_release(prim);
}

void test_prim_ctor_1x1v_p1() { check_prim(1, 1, 1); }
void test_prim_ctor_1x2v_p1() { check_prim(1, 2, 1); }
void test_prim_ctor_2x2v_p1() { check_prim(2, 2, 1); }
void test_prim_ctor_3x2v_p1() { check_prim(3, 2, 1); }
void test_prim_ctor_1x1v_p2() { check_prim(1, 1, 2); }
void test_prim_ctor_2x2v_p2() { check_prim(2, 2, 2); }

void
test_prim_ctor_acquire()
{
  struct gkyl_basis cbasis, pbasis;
  gkyl_cart_modal_serendip(&cbasis, 1, 1);
  gkyl_cart_modal_gkhybrid(&pbasis, 1, 2);

  struct gkyl_prim_lbo_type *prim = gkyl_prim_lbo_gyrokinetic_new(&cbasis, &pbasis, false);

  // Acquire a second reference; the underlying object must survive the first
  // release, and report identical fields.
  struct gkyl_prim_lbo_type *prim2 = gkyl_prim_lbo_type_acquire(prim);
  TEST_CHECK( prim2 == prim );

  gkyl_prim_lbo_type_release(prim);

  TEST_CHECK( prim2->cdim == 1 );
  TEST_CHECK( prim2->pdim == 3 );
  TEST_CHECK( prim2->udim == 1 );

  gkyl_prim_lbo_type_release(prim2);
}

TEST_LIST = {
  { "prim_ctor_1x1v_p1", test_prim_ctor_1x1v_p1 },
  { "prim_ctor_1x2v_p1", test_prim_ctor_1x2v_p1 },
  { "prim_ctor_2x2v_p1", test_prim_ctor_2x2v_p1 },
  { "prim_ctor_3x2v_p1", test_prim_ctor_3x2v_p1 },
  { "prim_ctor_1x1v_p2", test_prim_ctor_1x1v_p2 },
  { "prim_ctor_2x2v_p2", test_prim_ctor_2x2v_p2 },
  { "prim_ctor_acquire", test_prim_ctor_acquire },
  { NULL, NULL },
};
