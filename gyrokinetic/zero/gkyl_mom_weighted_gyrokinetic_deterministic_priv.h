#pragma once

// One block owns each configuration cell and reduces weighted velocity-cell
// contributions in a fixed order for every moment coefficient.

#include <cub/block/block_reduce.cuh>

enum { MOM_WEIGHTED_GK_BLOCK_SIZE = 128 };

__global__ static void
mom_weighted_gk_advance_deterministic_cu_ker(
  double mass, double charge, struct gkyl_mom_weighted_gyrokinetic_kernels *kernels,
  struct gkyl_rect_grid phase_grid, struct gkyl_range conf_rng, struct gkyl_range vel_rng,
  struct gkyl_range phase_rng, struct gkyl_range wgt_rng, enum gkyl_mom_weight_type wgt_type,
  const struct gkyl_array *bmag, const struct gkyl_array *vmap, const struct gkyl_array *wgt,
  const struct gkyl_array *phi, const struct gkyl_array *fin, struct gkyl_array *mout, long num_vel
)
{
  typedef cub::BlockReduce<double, MOM_WEIGHTED_GK_BLOCK_SIZE, cub::BLOCK_REDUCE_RAKING>
    BlockReduceT;
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
  // Each thread accumulates weighted moments from a fixed velocity stride.
  for (long v = threadIdx.x; v < num_vel; v += MOM_WEIGHTED_GK_BLOCK_SIZE) {
    long remaining = v;
    for (int d = conf_rng.ndim; d < phase_rng.ndim; ++d) {
      int extent = phase_rng.upper[d] - phase_rng.lower[d] + 1;
      pidx[d] = phase_rng.lower[d] + remaining % extent;
      vidx[d - conf_rng.ndim] = pidx[d];
      remaining /= extent;
    }
    const double *wgt_c =
      wgt ? (const double *)gkyl_array_cfetch(
              wgt, gkyl_range_idx(&wgt_rng, wgt_type == GKYL_F_MOMENT_WEIGHT_CONF ? cidx : vidx)
            ) :
            NULL;
    double moment[48] = {0.0};
    mom_weighted_gk_call(
      kernels->wmom, phase_grid.dx,
      (const double *)gkyl_array_cfetch(vmap, gkyl_range_idx(&vel_rng, vidx)), mass, charge, bmag_c,
      phi_c, (const double *)gkyl_array_cfetch(fin, gkyl_range_idx(&phase_rng, pidx)), wgt_c, moment
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
