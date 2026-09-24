/* -*- c -*- */

#include <stdio.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
  int cu_range_test(const struct gkyl_range rng);
}

__global__
void ker_cu_range_test(const struct gkyl_range rng, int *nfail)
{
  *nfail = 0;

  int lower[] = {0, 0}, upper[] = {24, 49};

  GKYL_CU_CHECK( rng.ndim == 2, nfail );
  GKYL_CU_CHECK( rng.volume == 25*50, nfail );

  for (unsigned i=0; i<2; ++i) {
    GKYL_CU_CHECK( rng.lower[i] == lower[i], nfail );
    GKYL_CU_CHECK( rng.upper[i] == upper[i], nfail );
  }  
}

int cu_range_test(const struct gkyl_range rng)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));  
  ker_cu_range_test<<<1,1>>>(rng, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;  
}


