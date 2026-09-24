// Gkyl ------------------------------------------------------------------------
//
// C back-end for use in multimomlinear Tool. This code computes
// eigenvalues and eigenvectors of the complex matrix supplied to it.
// 
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <gkyl_multi_mom_linear.h>
#include <complex.h>

// BLAS and LAPACKE includes
#ifdef GKYL_USING_FRAMEWORK_ACCELERATE
# include <Accelerate/Accelerate.h>
#else
// On non-Darwin platforms use OpenBLAS
# include <cblas.h>
# include <lapacke.h>
#endif

bool
gkyl_multi_mom_eigensolve(double _Complex *A, double _Complex *x, double _Complex *vl, double _Complex *vr, int N, int eig_vec)
{
  char calc_vec = 'N';
  calc_vec = (eig_vec == 1) ? 'V' : 'N';
#ifdef GKYL_USING_FRAMEWORK_ACCELERATE  
  int info = 1;
#else  
  int info = LAPACKE_zgeev(LAPACK_COL_MAJOR, calc_vec, calc_vec, N, A,
    N, x, vl, N, vr, N);
#endif
  return info == 0 ? true : false;
}
