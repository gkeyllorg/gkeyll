/* -*- c -*- */

#include <stdio.h>

extern "C" {
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_alloc.h>
  int cu_rect_grid_test(const struct gkyl_rect_grid grid);
}

__global__
void ker_cu_rect_grid_test(const struct gkyl_rect_grid grid, int *nfail)
{
  *nfail = 0;

  double lower[] = {1.0, 1.0}, upper[] = {2.5, 5.0};
  int cells[] = {20, 20};

  GKYL_CU_CHECK( grid.ndim == 2, nfail );
  for (int i=0; i<grid.ndim; ++i) {
    GKYL_CU_CHECK( grid.lower[i] == lower[i], nfail );
    GKYL_CU_CHECK( grid.upper[i] == upper[i], nfail );
    GKYL_CU_CHECK( grid.cells[i] == cells[i], nfail );
    GKYL_CU_CHECK( grid.dx[i] == (upper[i]-lower[i])/cells[i], nfail );
  }
  GKYL_CU_CHECK( grid.cellVolume == 0.075*0.2, nfail );
  
}

int cu_rect_grid_test(const struct gkyl_rect_grid grid)
{
  int *nfail_dev = (int *) gkyl_cu_malloc(sizeof(int));  
  ker_cu_rect_grid_test<<<1,1>>>(grid, nfail_dev);

  int nfail;
  gkyl_cu_memcpy(&nfail, nfail_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  gkyl_cu_free(nfail_dev);

  return nfail;  
}
