/* -*- c -*- */

#include <stdio.h>

extern "C" {
#include <gkyl_util.h>
#include <gkyl_alloc.h>
  int dev_cu_malloc_array(double **arr, int narr, int nelem);
}

__global__ void
ker_dev_cu_malloc_array(double **arr, int narr, int nelem, int *nfail)
{
  *nfail = 0;

  for (int k=0; k<narr; k++) {
    double *arr_d = arr[k];
    for (int i=0; i<nelem; i++)
      GKYL_CU_CHECK( arr_d[i] == (double)(nelem*k+i), nfail);
  }
}

int dev_cu_malloc_array(double **arr, int narr, int nelem)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));
  ker_dev_cu_malloc_array<<<1,1>>>(arr, narr, nelem, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;
}
