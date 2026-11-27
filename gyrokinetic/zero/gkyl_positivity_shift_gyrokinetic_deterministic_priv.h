#pragma once

// One block shifts the velocity cells of a configuration cell and reduces its
// old and new densities in a fixed order. The previous floor stays read-only.

#include <cub/block/block_reduce.cuh>

enum { POS_SHIFT_GK_BLOCK_SIZE = 128 };

struct pos_shift_gk_max {
  __device__ double operator()(double a, double b) const
  {
    return fmax(a, b);
  }
};

__global__ static void
pos_shift_gk_deterministic_cu_ker(
  struct gkyl_positivity_shift_gyrokinetic_kernels *kers, struct gkyl_rect_grid grid,
  struct gkyl_range conf_range, struct gkyl_range vel_range, struct gkyl_range phase_range,
  const double *ffloor, double mass, const struct gkyl_array *bmag,
  const struct gkyl_array *jacobtot, const struct gkyl_array *jacobtot_inv,
  const struct gkyl_array *vmap, const struct gkyl_array *jacobvel, struct gkyl_array *shiftedf,
  struct gkyl_array *distf, struct gkyl_array *m0, struct gkyl_array *delta_m0,
  double *max_buffer_dev, long num_vel
)
{
  typedef cub::BlockReduce<double, POS_SHIFT_GK_BLOCK_SIZE, cub::BLOCK_REDUCE_RAKING> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  int cidx[GKYL_MAX_CDIM], pidx[GKYL_MAX_DIM], vidx[GKYL_MAX_VDIM];
  gkyl_sub_range_inv_idx(&conf_range, blockIdx.x, cidx);
  for (int d = 0; d < conf_range.ndim; ++d) {
    pidx[d] = cidx[d];
  }
  long clin = gkyl_range_idx(&conf_range, cidx);
  const double *bmag_c = (const double *)gkyl_array_cfetch(bmag, clin);
  const double *jacobtot_c = (const double *)gkyl_array_cfetch(jacobtot, clin);
  const double *jacobtot_inv_c = (const double *)gkyl_array_cfetch(jacobtot_inv, clin);
  double old_local[20] = {0.0}, new_local[20] = {0.0};
  double local_max = -DBL_MAX, shifted_local = 0.0;
  // Process each thread's fixed velocity stride while the old floor is read-only.
  for (long v = threadIdx.x; v < num_vel; v += POS_SHIFT_GK_BLOCK_SIZE) {
    long remaining = v;
    for (int d = conf_range.ndim; d < phase_range.ndim; ++d) {
      int extent = phase_range.upper[d] - phase_range.lower[d] + 1;
      pidx[d] = phase_range.lower[d] + remaining % extent;
      vidx[d - conf_range.ndim] = pidx[d];
      remaining /= extent;
    }
    long plin = gkyl_range_idx(&phase_range, pidx);
    const double *vmap_c =
      (const double *)gkyl_array_cfetch(vmap, gkyl_range_idx(&vel_range, vidx));
    const double *jacobvel_c = (const double *)gkyl_array_cfetch(jacobvel, plin);
    double *distf_c = (double *)gkyl_array_fetch(distf, plin);
    double old_m0[20] = {0.0};
    kers->m0(grid.dx, vmap_c, mass, 0.0, bmag_c, NULL, distf_c, old_m0);
    for (unsigned k = 0; k < delta_m0->ncomp; ++k) {
      old_local[k] += old_m0[k];
    }
    // Shift the physical distribution, then restore its geometric factors.
    kers->conf_phase_mul_op(jacobtot_inv_c, distf_c, distf_c);
    for (unsigned k = 0; k < distf->ncomp; ++k) {
      distf_c[k] /= jacobvel_c[0];
    }
    bool shifted = kers->shift(ffloor[0], distf_c);
    kers->conf_phase_mul_op(jacobtot_c, distf_c, distf_c);
    for (unsigned k = 0; k < distf->ncomp; ++k) {
      distf_c[k] *= jacobvel_c[0];
    }
    if (shifted) {
      double new_m0[20] = {0.0};
      kers->m0(grid.dx, vmap_c, mass, 0.0, bmag_c, NULL, distf_c, new_m0);
      // Preserve the cell-averaged density when both local densities are positive.
      if (old_m0[0] > 0.0 && new_m0[0] > 0.0) {
        double ratio = old_m0[0] / new_m0[0];
        for (unsigned k = 0; k < distf->ncomp; ++k) {
          distf_c[k] *= ratio;
        }
        for (unsigned k = 0; k < m0->ncomp; ++k) {
          new_local[k] += ratio * new_m0[k];
        }
      } else {
        for (unsigned k = 0; k < m0->ncomp; ++k) {
          new_local[k] += new_m0[k];
        }
        shifted_local = 1.0;
      }
    } else {
      for (unsigned k = 0; k < m0->ncomp; ++k) {
        new_local[k] += old_m0[k];
      }
    }
    local_max = fmax(local_max, distf_c[0]);
  }
  // Reduce old and new density coefficients separately for this configuration cell.
  double *m0_c = (double *)gkyl_array_fetch(m0, clin);
  double *delta_c = (double *)gkyl_array_fetch(delta_m0, clin);
  for (unsigned k = 0; k < m0->ncomp; ++k) {
    double sum = BlockReduceT(temp).Sum(new_local[k]);
    if (threadIdx.x == 0) {
      m0_c[k] = sum;
    }
    __syncthreads();
    sum = BlockReduceT(temp).Sum(old_local[k]);
    if (threadIdx.x == 0) {
      delta_c[k] = sum;
    }
    __syncthreads();
  }
  // Record whether a non-rescalable shift occurred and the largest updated value.
  double shifted_sum = BlockReduceT(temp).Sum(shifted_local);
  if (threadIdx.x == 0) {
    *(int *)gkyl_array_fetch(shiftedf, clin) = shifted_sum > 0.0;
  }
  __syncthreads();
  double max_value = BlockReduceT(temp).Reduce(local_max, pos_shift_gk_max());
  if (threadIdx.x == 0) {
    max_buffer_dev[blockIdx.x] = max_value;
  }
}

// Update the floor only after all phase cells have finished using its old value.
__global__ static void
pos_shift_gk_update_floor_cu_ker(
  double *ffloor, const double *max_buffer_dev, long count, double ffloor_fac, double cellav_fac
)
{
  if (threadIdx.x == 0) {
    double maximum = -DBL_MAX;
    for (long i = 0; i < count; ++i) {
      maximum = fmax(maximum, max_buffer_dev[i]);
    }
    ffloor[0] = ffloor_fac * maximum * cellav_fac;
  }
}
