#include <acutest.h>

#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_maxwell_priv.h>

void
test_dg_max()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);

  struct gkyl_dg_eqn* eqn = gkyl_dg_maxwell_new(&basis, 1.0, 0.5, 0.25, false);

  TEST_CHECK( eqn->num_equations == 8 );

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);

  TEST_CHECK( maxwell->maxwell_data.c == 1.0 );
  TEST_CHECK( maxwell->maxwell_data.chi == 0.5 );
  TEST_CHECK( maxwell->maxwell_data.gamma == 0.25 );

  struct dg_maxwell *m_on_dev = container_of(eqn->on_dev,
    struct dg_maxwell, eqn);

  TEST_CHECK( m_on_dev->maxwell_data.c == 1.0 );
  TEST_CHECK( m_on_dev->maxwell_data.chi == 0.5 );
  TEST_CHECK( m_on_dev->maxwell_data.gamma == 0.25 );

  gkyl_dg_eqn_release(eqn);
}

#ifdef GKYL_HAVE_CUDA

int cu_maxwell_test(const struct gkyl_dg_eqn *eqn);

void
test_cu_dg_max()
{
  struct gkyl_basis basis;
  gkyl_cart_modal_serendip(&basis, 1, 1);

  struct gkyl_dg_eqn* eqn = gkyl_dg_maxwell_cu_dev_new(&basis, 1.0, 0.5, 0.25);

  // this is not possible from user code and should NOT be done. This
  // is for testing only
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);

  TEST_CHECK( maxwell->maxwell_data.c == 1.0 );
  TEST_CHECK( maxwell->maxwell_data.chi == 0.5 );
  TEST_CHECK( maxwell->maxwell_data.gamma == 0.25 );
  
  // call CUDA test
  /* int nfail = cu_maxwell_test(eqn->on_dev); */

  /* TEST_CHECK( nfail == 0 ); */

  gkyl_dg_eqn_release(eqn);
}

#endif

TEST_LIST = {
  { "dg_max", test_dg_max },
#ifdef GKYL_HAVE_CUDA
  { "cu_dg_max", test_cu_dg_max },
#endif  
  { NULL, NULL },
};
