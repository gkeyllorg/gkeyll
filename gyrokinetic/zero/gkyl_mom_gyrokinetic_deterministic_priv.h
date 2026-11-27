#pragma once

// One block owns each configuration cell. Its threads visit velocity cells in
// fixed strides and reduce each moment coefficient with a fixed CUB tree.

#include <cub/block/block_reduce.cuh>

enum { MOM_GK_BLOCK_SIZE = 128 };

__global__ static void
mom_gk_advance_deterministic_cu_ker(
  double mass, double charge, struct gkyl_mom_gyrokinetic_kernels *kernels,
  struct gkyl_rect_grid phase_grid, struct gkyl_range conf_rng, struct gkyl_range vel_rng,
  struct gkyl_range phase_rng, const struct gkyl_array *bmag, const struct gkyl_array *vmap,
  const struct gkyl_array *phi, const struct gkyl_array *fin, struct gkyl_array *mout, long num_vel
)
{
  typedef cub::BlockReduce<double, MOM_GK_BLOCK_SIZE, cub::BLOCK_REDUCE_RAKING> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  int cidx[GKYL_MAX_CDIM], pidx[GKYL_MAX_DIM], vidx[GKYL_MAX_VDIM];
  gkyl_sub_range_inv_idx(&conf_rng, blockIdx.x, cidx);
  for (int d = 0; d < conf_rng.ndim; ++d) {
    pidx[d] = cidx[d];
  }
  long clin = gkyl_range_idx(&conf_rng, cidx);
  const double *bmag_c = (const double *)gkyl_array_cfetch(bmag, clin);
  const double *phi_c = phi ? (const double *)gkyl_array_cfetch(phi, clin) : NULL;
  double local[48] = {0.0};
  // Each thread accumulates its fixed stride of velocity cells.
  for (long v = threadIdx.x; v < num_vel; v += MOM_GK_BLOCK_SIZE) {
    long remaining = v;
    for (int d = conf_rng.ndim; d < phase_rng.ndim; ++d) {
      int extent = phase_rng.upper[d] - phase_rng.lower[d] + 1;
      pidx[d] = phase_rng.lower[d] + remaining % extent;
      vidx[d - conf_rng.ndim] = pidx[d];
      remaining /= extent;
    }
    double moment[48] = {0.0};
    mom_gk_call(
      kernels->wmom, phase_grid.dx,
      (const double *)gkyl_array_cfetch(vmap, gkyl_range_idx(&vel_rng, vidx)), mass, charge, bmag_c,
      phi_c, (const double *)gkyl_array_cfetch(fin, gkyl_range_idx(&phase_rng, pidx)), moment
    );
    for (unsigned k = 0; k < mout->ncomp; ++k) {
      local[k] += moment[k];
    }
  }
  // Combine thread sums and write this configuration cell once.
  double *out = (double *)gkyl_array_fetch(mout, clin);
  for (unsigned k = 0; k < mout->ncomp; ++k) {
    double sum = BlockReduceT(temp).Sum(local[k]);
    if (threadIdx.x == 0) {
      out[k] = sum;
    }
    __syncthreads();
  }
}
