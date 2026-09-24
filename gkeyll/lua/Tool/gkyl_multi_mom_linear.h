// Gkyl ------------------------------------------------------------------------
//
// C back-end for use in multimomlinear Tool. This code computes
// eigenvalues and eigenvectors of the complex matrix supplied to it.
// 
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

// std includes
#include <complex.h>
#include <stdbool.h>

bool gkyl_multi_mom_eigensolve(double _Complex *A, double _Complex *x, double _Complex *vl, double _Complex *vr, int N, int eig_vec);
