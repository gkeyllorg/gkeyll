#pragma once

#include <math.h>
#include <stdlib.h>

// Maximum order ordinate/weight data
enum { gkyl_gauss_max = 8 };

// Ordinates
static const double gkyl_gauss_ordinates_1[] =
{ 0.0 };
static const double gkyl_gauss_ordinates_2[] =
{ -0.5773502691896257645091, 0.5773502691896257645091 };
static const double gkyl_gauss_ordinates_3[] =
{ -0.7745966692414833770359, 0, 0.7745966692414833770359 };
static const double gkyl_gauss_ordinates_4[] =
{ -0.8611363115940525752239, -0.3399810435848562648027, 0.3399810435848562648027, 0.8611363115940525752239 };
static const double gkyl_gauss_ordinates_5[] =
{ -0.9061798459386639927976, -0.5384693101056830910363, 0, 0.5384693101056830910363, 0.9061798459386639927976 };
static const double gkyl_gauss_ordinates_6[] =
{ -0.9324695142031520278123, -0.6612093864662645136614, -0.2386191860831969086305, 0.2386191860831969086305, 0.6612093864662645136614, 0.9324695142031520278123 };
static const double gkyl_gauss_ordinates_7[] =
{ -0.9491079123427585245262, -0.7415311855993944398639, -0.4058451513773971669066, 0, 0.4058451513773971669066, 0.7415311855993944398639, 0.9491079123427585245262 };
static const double gkyl_gauss_ordinates_8[] =
{ -0.960289856497536231684, -0.7966664774136267395916, -0.5255324099163289858177, -0.1834346424956498049395, 0.1834346424956498049395, 0.525532409916328985818, 0.796666477413626739592, 0.9602898564975362316836 };

// Weights
static const double gkyl_gauss_weights_1[] =
{ 2.0 };
static const double gkyl_gauss_weights_2[] =
{ 1.0, 1.0 };
static const double gkyl_gauss_weights_3[] =
{ 0.5555555555555555555556, 0.888888888888888888889, 0.555555555555555555556 };
static const double gkyl_gauss_weights_4[] =
{ 0.3478548451374538573731, 0.6521451548625461426269, 0.652145154862546142627, 0.3478548451374538573731 };
static const double gkyl_gauss_weights_5[] =
{ 0.2369268850561890875143, 0.4786286704993664680413, 0.568888888888888888889, 0.478628670499366468041, 0.236926885056189087514 };
static const double gkyl_gauss_weights_6[] =
{ 0.1713244923791703450403, 0.36076157304813860757, 0.46791393457269104739, 0.46791393457269104739, 0.36076157304813860757, 0.1713244923791703450403 };
static const double gkyl_gauss_weights_7[] =
{ 0.1294849661688696932706, 0.279705391489276667901, 0.38183005050511894495, 0.4179591836734693877552, 0.38183005050511894495, 0.279705391489276667901, 0.1294849661688696932706 };
static const double gkyl_gauss_weights_8[] =
{ 0.1012285362903762591525, 0.222381034453374470544, 0.313706645877887287338, 0.3626837833783619829652, 0.3626837833783619829652, 0.31370664587788728734, 0.222381034453374470544, 0.1012285362903762591525 };

// gkyl_gauss_ordinates[N] are ordinates for N-point Guassian
// integration
static const double* gkyl_gauss_ordinates[] = {
  0, // N=0 makes no sense,
  gkyl_gauss_ordinates_1,
  gkyl_gauss_ordinates_2,
  gkyl_gauss_ordinates_3,
  gkyl_gauss_ordinates_4,
  gkyl_gauss_ordinates_5,
  gkyl_gauss_ordinates_6,
  gkyl_gauss_ordinates_7,
  gkyl_gauss_ordinates_8
};

// gkyl_gauss_weights[N] are weights for N-point Guassian integration
static const double *gkyl_gauss_weights[] = {
  0, // N=0 makes no sense,
  gkyl_gauss_weights_1,
  gkyl_gauss_weights_2,
  gkyl_gauss_weights_3,
  gkyl_gauss_weights_4,
  gkyl_gauss_weights_5,
  gkyl_gauss_weights_6,
  gkyl_gauss_weights_7,
  gkyl_gauss_weights_8
};


// Lobatto quadrature

// Ordinates
static const double gkyl_gauss_lobatto_ordinates_2[] = 
{ -1.0,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_3[] = 
{ -1.0,0.0,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_4[] = 
{ -1.0,-0.4472135954999579,0.4472135954999579,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_5[] = 
{ -1.0,-0.6546536707079771,0.0,0.6546536707079771,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_6[] = 
{ -1.0,-0.7650553239294646,-0.285231516480645,0.285231516480645,0.7650553239294646,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_7[] = 
{ -1.0,-0.8302238962785669,-0.4688487934707142,0.0,0.4688487934707142,0.8302238962785669,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_8[] = 
{ -1.0,-0.8717401485096066,-0.5917001814331423,-0.2092992179024789,0.2092992179024789,0.5917001814331423,0.8717401485096066,1.0 };

// Weights
static const double gkyl_gauss_lobatto_weights_2[] = 
{ 1.0,1.0 }; 
static const double gkyl_gauss_lobatto_weights_3[] = 
{ 0.3333333333333333,1.333333333333333,0.3333333333333333 }; 
static const double gkyl_gauss_lobatto_weights_4[] = 
{ 0.1666666666666667,0.8333333333333334,0.8333333333333334,0.1666666666666667 }; 
static const double gkyl_gauss_lobatto_weights_5[] = 
{ 0.1,0.5444444444444444,0.7111111111111111,0.5444444444444444,0.1 }; 
static const double gkyl_gauss_lobatto_weights_6[] = 
{ 0.06666666666666667,0.378474956297847,0.5548583770354863,0.5548583770354863,0.378474956297847,0.06666666666666667 }; 
static const double gkyl_gauss_lobatto_weights_7[] = 
{ 0.04761904761904762,0.276826047361566,0.4317453812098626,0.4876190476190476,0.4317453812098626,0.276826047361566,0.04761904761904762 }; 
static const double gkyl_gauss_lobatto_weights_8[] = 
{ 0.03571428571428571,0.210704227143506,0.3411226924835044,0.4124587946587039,0.4124587946587039,0.3411226924835044,0.210704227143506,0.03571428571428571 };

// gkyl_gauss_lobatto_ordinates[N] are ordinates for N-point
// Guass-Lobatto integration
static const double* gkyl_gauss_lobatto_ordinates[] = {
  0, // N=0 makes no sense,
  0, // N=1 makes no sense,
  gkyl_gauss_lobatto_ordinates_2,
  gkyl_gauss_lobatto_ordinates_3,
  gkyl_gauss_lobatto_ordinates_4,
  gkyl_gauss_lobatto_ordinates_5,
  gkyl_gauss_lobatto_ordinates_6,
  gkyl_gauss_lobatto_ordinates_7,
  gkyl_gauss_lobatto_ordinates_8
};

// gkyl_gauss_lobatto_weights[N] are weights for N-point Guass-Lobatto
// integration
static const double *gkyl_gauss_lobatto_weights[] = {
  0, // N=0 makes no sense,
  0, // N=1 makes no sense,  
  gkyl_gauss_lobatto_weights_2,
  gkyl_gauss_lobatto_weights_3,
  gkyl_gauss_lobatto_weights_4,
  gkyl_gauss_lobatto_weights_5,
  gkyl_gauss_lobatto_weights_6,
  gkyl_gauss_lobatto_weights_7,
  gkyl_gauss_lobatto_weights_8
};

/**
 * Compute ordinates and weights for use in Gaussian quadrature.
 *
 * @param x1 Left coordinate of domain.
 * @param x2 Right coordinate of domain.
 * @param x On output, ordinates.
 * @param w On output, weights.
 * @param n Order of the quadrature.
 */
void gkyl_gauleg(double x1, double x2,  double x[], double w[], int n);

/**
 * Compute ordinates and weights for use in n-dimensiona Gaussian
 * quadrature of order q. The array "x" needs to be at least
 * ndim*pow(nq,ndim) long and the array "w" needs to be at least
 * pow(nq,ndim) long.
 *
 * @param x On output, ordinates.
 * @param w On output, weights.
 * @param nq Order of the quadrature.
 */
void gkyl_ndim_ordinates_weights(int ndim, double *x, double *w, int nq);
