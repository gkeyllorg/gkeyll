#include <gkyl_const.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_range.h>

#include <string.h>

#define GUASS_QUAD_EPS 3.0e-15

// This is based on an implementation in Numerical Recipes in C book
static void
priv_gkyl_gauleg( double x1, double x2,  double x[], double w[], int n)
{
  double z1, xm, xl, pp, p3, p2, p1;
  int m = (n+1)/2;
  
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);
  for (int i = 1; i <= m; i++) {
    double z = cos(GKYL_PI*(i-0.25)/(n+0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (int j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp = n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z = z1-p1/pp;
    } while( fabs(z-z1) > GUASS_QUAD_EPS );
    x[i] = xm-xl*z;
    x[n+1-i] = xm+xl*z;
    w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i] = w[i];
  }
}

void
gkyl_gauleg(double x1, double x2,  double x[], double w[], int n)
{
  priv_gkyl_gauleg(x1, x2, x-1, w-1, n); // actual routine assumes 1-offset arrays
}

void
gkyl_ndim_ordinates_weights(int ndim, double *x, double *w, int nq)
{
  double ordinates1[nq], weights1[nq];

  if (nq <= gkyl_gauss_max) {
    // use pre-computed values if possible (these are more accurate
    // than computing them on the fly)
    memcpy(ordinates1, gkyl_gauss_ordinates[nq], sizeof(double[nq]));
    memcpy(weights1, gkyl_gauss_weights[nq], sizeof(double[nq]));
  }
  else {
    gkyl_gauleg(-1, 1, ordinates1, weights1, nq);
  }
  
  int shape[ndim];
  for (int d=0; d<ndim; ++d) shape[d] = nq;
  
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, ndim, shape);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&qrange, iter.idx);

    for (int d=0; d<ndim; ++d)
      x[lidx*ndim+d] = ordinates1[iter.idx[d]];

    w[lidx] = 1.0;
    for (int d=0; d<ndim; ++d)
      w[lidx] *= weights1[iter.idx[d]];
  }
}
