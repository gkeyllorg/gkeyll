/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_priv.h>
#include <gkyl_mom_type.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_mom_calc_advance_cu_ker(const gkyl_mom_calc* mcalc,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_array* GKYL_RESTRICT fin, struct gkyl_array* GKYL_RESTRICT mout)
{
  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    gkyl_rect_grid_cell_center(&mcalc->grid, pidx, xc);

    long lincP = gkyl_range_idx(&phase_range, pidx);
    const double* fptr = (const double*) gkyl_array_cfetch(fin, lincP);
    double momLocal[96]; // hard-coded to 3 * max confBasis.num_basis (3x p=3 Ser) for now.
    for (unsigned int k=0; k<96; ++k)
      momLocal[k] = 0.0;

    // reduce local f to local mom
    mcalc->momt->kernel(mcalc->momt, xc, mcalc->grid.dx, pidx, fptr, &momLocal[0], 0);

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_range.ndim; k++)
      cidx[k] = pidx[k];
    long lincC = gkyl_range_idx(&conf_range, cidx);

    double* mptr = (double*) gkyl_array_fetch(mout, lincC);
    for (unsigned int k = 0; k < mout->ncomp; ++k) {
       atomicAdd(&mptr[k], momLocal[k]);
    }
  }
}

void
gkyl_mom_calc_advance_cu(const gkyl_mom_calc* mcalc,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT mout)
{
  int nblocks = phase_range->nblocks, nthreads = phase_range->nthreads;
  gkyl_array_clear_range(mout, 0.0, conf_range);
  gkyl_mom_calc_advance_cu_ker<<<nblocks, nthreads>>>
    (mcalc->on_dev, *phase_range, *conf_range, fin->on_dev, mout->on_dev);
}

gkyl_mom_calc*
gkyl_mom_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt)
{
  gkyl_mom_calc *up = (gkyl_mom_calc*) gkyl_malloc(sizeof(gkyl_mom_calc));
  up->grid = *grid;

  struct gkyl_mom_type *mt = gkyl_mom_type_acquire(momt);
  up->momt = mt->on_dev; // so memcpy below gets dev copy

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  gkyl_mom_calc *up_cu = (gkyl_mom_calc*) gkyl_cu_malloc(sizeof(gkyl_mom_calc));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_mom_calc), GKYL_CU_MEMCPY_H2D);

  up->momt = mt; // host portion of struct should have host copy
  up->on_dev = up_cu; // host pointer

  return up;
}
