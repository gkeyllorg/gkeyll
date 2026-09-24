/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_interpolate.h>
#include <gkyl_dg_interpolate_priv.h>
#include <gkyl_util.h>
}

// CUDA kernel to set device pointer to interpolating kernel.
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_interp_choose_kernel_ptrs_cu(struct gkyl_dg_interpolate_kernels *kernels, int cdim, struct gkyl_basis basis, int dir, double dxRat)
{
  kernels->interp = dg_interp_choose_gk_interp_kernel(cdim, basis, dir);

  // Map from grid to stencil index in each direction.
  if (dxRat > 1)
    kernels->grid2stencil = dg_interp_index_stencil_map_refine;
  else
    kernels->grid2stencil = dg_interp_index_stencil_map_coarsen;
}

void
dg_interp_choose_kernel_cu(struct gkyl_dg_interpolate_kernels *kernels, int cdim, struct gkyl_basis basis, int dir, double dxRat)
{
  dg_interp_choose_kernel_ptrs_cu<<<1,1>>>(kernels, cdim, basis, dir, dxRat); 
}

__global__ static void
gkyl_dg_interpolate_advance_1x_cu_ker(struct gkyl_dg_interpolate_kernels *kernels,
  int dir, double dxRat, int *offset_upper, 
  struct gkyl_rect_grid grid_do, struct gkyl_rect_grid grid_tar,
  struct gkyl_range range_do, struct gkyl_range range_tar,
  const struct gkyl_array* GKYL_RESTRICT fdo, struct gkyl_array* GKYL_RESTRICT ftar)
{
  int idx_do[GKYL_MAX_DIM] = {-1};
  int idx_tar[GKYL_MAX_DIM] = {-1};
  int idx_tar_lo;
  double xc_do[GKYL_MAX_DIM];
  double xc_tar[GKYL_MAX_DIM];

  int ndim = range_do.ndim;

  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < range_do.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&range_do, tid, idx_do);

    gkyl_rect_grid_cell_center(&grid_do, idx_do, xc_do);

    long linidx_do = gkyl_range_idx(&range_do, idx_do);
    const double *fdo_c = (const double *) gkyl_array_cfetch(fdo, linidx_do);

    // Compute the index of the lower cell this cell contributes to.
    double eveOI = dxRat*(idx_do[dir]-1);
    idx_tar_lo = ceil(eveOI)+((int) ceil(eveOI-floor(eveOI))+1) % 2;

    // Get the index to the stencil for this donor cell.
    int idx_sten = kernels->grid2stencil(idx_do[dir], grid_do.cells[dir], dxRat);

    for (int d=0; d<ndim; d++)
      idx_tar[d] = idx_do[d];

    // Loop over the target-grid cells this donor cell contributes to.
    for (int off=0; off<offset_upper[idx_sten]; off++) {

      idx_tar[dir] = idx_tar_lo + off;

      gkyl_rect_grid_cell_center(&grid_tar, idx_tar, xc_tar);

      long linidx_tar = gkyl_range_idx(&range_tar, idx_tar);
      double *ftar_c = (double *) gkyl_array_fetch(ftar, linidx_tar);

      kernels->interp(xc_do, xc_tar, grid_do.dx, grid_tar.dx, fdo_c, ftar_c);

    }
  }
}


void
gkyl_dg_interpolate_advance_1x_cu(gkyl_dg_interpolate* up,
  const struct gkyl_range *range_do, const struct gkyl_range *range_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar)
{
  gkyl_array_clear_range(ftar, 0.0, range_tar);

  int nblocks = range_do->nblocks, nthreads = range_do->nthreads;

  gkyl_dg_interpolate_advance_1x_cu_ker<<<nblocks, nthreads>>>
    (up->kernels, up->dir, up->dxRat, up->offset_upper, up->grid_do, up->grid_tar,
     *range_do, *range_tar, fdo->on_dev, ftar->on_dev);
}
