// GPU sparse linear solvers (gkyl_culinsolver_ops): cuDSS when available, cuSOLVER otherwise. The tests
// themselves are device code in ctest_cudss.cu and ctest_cusolver.cu; the host SuperLU solver is tested in
// moments/unit/ctest_linsolvers.c.
#include <acutest.h>

#ifdef GKYL_HAVE_CUDA
#ifdef GKYL_HAVE_CUDSS
void test_cudss_simple_dev();
void test_cudss_ops_dev();
void test_cudss_ops_update_amat_dev();
void test_cudss_ops_multiple_rhs_dev();
void test_cudss_ops_update_amat_multiple_rhs_dev();
#else
void test_cusolver_qr_dev();
void test_cusolver_rf_dev();
void test_cusolver_ops_dev();
void test_cusolver_ops_multiple_rhs_dev();
void test_cusolver_ops_multiple_prob_dev();
#endif
#endif

TEST_LIST = {
#ifdef GKYL_HAVE_CUDA
#ifdef GKYL_HAVE_CUDSS
  {"cudss_simple_dev", test_cudss_simple_dev},
  {"cudss_ops_dev", test_cudss_ops_dev},
  {"cudss_ops_update_amat_dev", test_cudss_ops_update_amat_dev},
  {"cudss_ops_multiple_rhs_dev", test_cudss_ops_multiple_rhs_dev},
  {"cudss_ops_update_amat_multiple_rhs_dev", test_cudss_ops_update_amat_multiple_rhs_dev},
#else
  {"cusolver_qr_dev", test_cusolver_qr_dev},
  {"cusolver_rf_dev", test_cusolver_rf_dev},
  {"cusolver_ops_dev", test_cusolver_ops_dev},
  {"cusolver_ops_multiple_rhs_dev", test_cusolver_ops_multiple_rhs_dev},
  {"cusolver_ops_multiple_prob_dev", test_cusolver_ops_multiple_prob_dev},
#endif
#endif
  {NULL, NULL}
};
