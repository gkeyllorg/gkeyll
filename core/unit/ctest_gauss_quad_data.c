// Unit tests for Gauss-Legendre quadrature data in gkyl_gauss_quad_data.h
#include <acutest.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>
#include <math.h>

// Integrate x^p over [-1,1] using N-point Gauss-Legendre quadrature.
static double
quad_int_monomial(int N, int p)
{
  const double *x = gkyl_gauss_ordinates[N];
  const double *w = gkyl_gauss_weights[N];
  double sum = 0.0;
  for (int i=0; i<N; ++i)
    sum += w[i]*pow(x[i], p);
  return sum;
}

void
test_gauss_weights_sum()
{
  // Sum of weights = integral of 1 over [-1,1] = 2, for every N.
  for (int N=1; N<=gkyl_gauss_max; ++N) {
    const double *w = gkyl_gauss_weights[N];
    double sum = 0.0;
    for (int i=0; i<N; ++i) sum += w[i];
    TEST_CHECK( gkyl_compare_double(sum, 2.0, 1e-13) );
    TEST_MSG("N=%d sum=%g", N, sum);
  }
}

void
test_gauss_ordinates_symmetric()
{
  // Ordinates are symmetric about 0; reversed list = negated list.
  for (int N=1; N<=gkyl_gauss_max; ++N) {
    const double *x = gkyl_gauss_ordinates[N];
    for (int i=0; i<N; ++i)
      TEST_CHECK( gkyl_compare_double(x[i], -x[N-1-i], 1e-13) );
  }
}

void
test_gauss_ordinates_in_range()
{
  // All ordinates lie strictly inside (-1,1).
  for (int N=1; N<=gkyl_gauss_max; ++N) {
    const double *x = gkyl_gauss_ordinates[N];
    for (int i=0; i<N; ++i)
      TEST_CHECK( x[i] > -1.0 && x[i] < 1.0 );
  }
}

void
test_gauss_exactness()
{
  // N-point Gauss-Legendre integrates polynomials up to degree 2N-1 exactly.
  // Integral of x^p over [-1,1] = 0 (odd p) or 2/(p+1) (even p).
  for (int N=1; N<=gkyl_gauss_max; ++N) {
    int maxdeg = 2*N - 1;
    for (int p=0; p<=maxdeg; ++p) {
      double exact = (p%2==1) ? 0.0 : 2.0/(p+1);
      double approx = quad_int_monomial(N, p);
      TEST_CHECK( gkyl_compare_double(approx, exact, 1e-11) );
      TEST_MSG("N=%d p=%d exact=%g approx=%g", N, p, exact, approx);
    }
  }
}

void
test_gauss_specific_values()
{
  // 2-point Gauss: ordinates +/- 1/sqrt(3), weights 1,1.
  const double *x2 = gkyl_gauss_ordinates[2];
  const double *w2 = gkyl_gauss_weights[2];
  TEST_CHECK( gkyl_compare_double(fabs(x2[0]), 1.0/sqrt(3.0), 1e-12) );
  TEST_CHECK( gkyl_compare_double(w2[0], 1.0, 1e-12) );
  TEST_CHECK( gkyl_compare_double(w2[1], 1.0, 1e-12) );

  // 1-point Gauss: ordinate 0, weight 2.
  TEST_CHECK( gkyl_compare_double(gkyl_gauss_ordinates[1][0], 0.0, 1e-13) );
  TEST_CHECK( gkyl_compare_double(gkyl_gauss_weights[1][0], 2.0, 1e-13) );

  // 3-point Gauss: middle ordinate 0, weights 5/9, 8/9, 5/9.
  const double *x3 = gkyl_gauss_ordinates[3];
  const double *w3 = gkyl_gauss_weights[3];
  TEST_CHECK( gkyl_compare_double(x3[1], 0.0, 1e-12) );
  TEST_CHECK( gkyl_compare_double(w3[1], 8.0/9.0, 1e-12) );
  TEST_CHECK( gkyl_compare_double(w3[0], 5.0/9.0, 1e-12) );
}

TEST_LIST = {
  { "gauss_weights_sum", test_gauss_weights_sum },
  { "gauss_ordinates_symmetric", test_gauss_ordinates_symmetric },
  { "gauss_ordinates_in_range", test_gauss_ordinates_in_range },
  { "gauss_exactness", test_gauss_exactness },
  { "gauss_specific_values", test_gauss_specific_values },
  { NULL, NULL },
};
