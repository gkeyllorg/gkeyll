#include <acutest.h>

#ifdef GKYL_HAVE_CUDA
#ifdef GKYL_HAVE_CUDSS
void test_cudss_simple_dev();
void test_cudss_ops_dev();
void test_cudss_ops_update_amat_dev();
void test_cudss_ops_multiple_rhs_dev();
#else
void test_cusolver_qr();
void test_cusolver_rf();
void test_cusolver_ops();
void test_cusolver_ops_multiple_rhs();
void test_cusolver_ops_multiple_prob();
#endif
#endif

TEST_LIST = {
#ifdef GKYL_HAVE_CUDA
#ifdef GKYL_HAVE_CUDSS
  {"cudss_simple_dev", test_cudss_simple_dev},
  {"cudss_ops_dev", test_cudss_ops_dev},
  {"cudss_ops_update_amat_dev", test_cudss_ops_update_amat_dev},
  {"cudss_ops_multiple_rhs_dev", test_cudss_ops_multiple_rhs_dev},
#else
  {"cusolver_qr_dev", test_cusolver_qr},
  {"cusolver_rf_dev", test_cusolver_rf},
  {"cusolver_ops_dev", test_cusolver_ops},
  {"cusolver_ops_multiple_rhs_dev", test_cusolver_ops_multiple_rhs},
  {"cusolver_ops_multiple_prob_dev", test_cusolver_ops_multiple_prob},
#endif
#endif
  {NULL, NULL}
};
