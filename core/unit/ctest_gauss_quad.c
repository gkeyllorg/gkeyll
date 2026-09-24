#include <acutest.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_util.h>
#include <gkyl_range.h>

static void
test_g1()
{
  int gauss_max = gkyl_gauss_max;
  double w[gauss_max], x[gauss_max];
  
  for (unsigned n=1; n<=gkyl_gauss_max; ++n) {
    // use the generic routine to get the ordinates and weights
    gkyl_gauleg(-1, 1, &x[0], &w[0], n);

    // compare with pre-computed values
    const double *xp = gkyl_gauss_ordinates[n];
    const double *wp = gkyl_gauss_weights[n];
    
    for (unsigned i=0; i<n; ++i) {
      TEST_CHECK( gkyl_compare(x[i], xp[i], 1e-12) );
      TEST_CHECK( gkyl_compare(w[i], wp[i], 1e-12) );
    }
  }
}

static void
test_ndim()
{
  int nquad = 4;
  double *x = gkyl_malloc(sizeof(double)*nquad*nquad*2);
  double *w = gkyl_malloc(sizeof(double)*nquad*nquad);

  gkyl_ndim_ordinates_weights(2, x, w, nquad);

  const double *xp = gkyl_gauss_ordinates[nquad];
  const double *wp = gkyl_gauss_weights[nquad];

  int shape[2] = { nquad, nquad };
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, 2, shape);

  double area = 0.0;
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);
  while (gkyl_range_iter_next(&iter)) {
    long lidx = gkyl_range_idx(&qrange, iter.idx);

    for (int d=0; d<2; ++d)
      TEST_CHECK( gkyl_compare_double(xp[iter.idx[d]], x[2*lidx+d], 1e-14) );
    
    area += w[lidx];
  }

  TEST_CHECK( gkyl_compare_double(4.0, area, 1e-14) );
  
  gkyl_free(x);
  gkyl_free(w);
}

TEST_LIST = {
  { "basic", test_g1 },
  { "ndim", test_ndim },
  { NULL, NULL },
};
