#pragma once

// Sum each SR pressure component over ordered velocity tiles, then reduce
// the tile results without atomic additions.

#include <cub/cub.cuh>

enum { SR_PRESSURE_TILE_SIZE = 128 };

__global__ static void
dg_calc_sr_vars_pressure_partial_cu_kernel(
  struct gkyl_dg_calc_sr_vars *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  bool use_vmap, const struct gkyl_array *vmap, const struct gkyl_array *gamma,
  const struct gkyl_array *gamma_inv, const struct gkyl_array *u_i, const struct gkyl_array *u_i_sq,
  const struct gkyl_array *GammaV, const struct gkyl_array *GammaV_sq, const struct gkyl_array *f,
  double *partial_dev, long conf_start, long num_tiles, long num_vel, unsigned int ncomp
)
{
  typedef cub::BlockReduce<double, SR_PRESSURE_TILE_SIZE, cub::BLOCK_REDUCE_RAKING> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  long conf_local = conf_start + blockIdx.x / num_tiles;
  long tile = blockIdx.x % num_tiles;
  long vel_local = tile * SR_PRESSURE_TILE_SIZE + threadIdx.x;
  int idx[GKYL_MAX_DIM];
  // Pair this velocity tile with its configuration cell.
  gkyl_sub_range_inv_idx(&conf_range, conf_local, idx);

  double mom_local[32] = {0.0};
  if (vel_local < num_vel) {
    long remaining = vel_local;
    for (int d = conf_range.ndim; d < phase_range.ndim; ++d) {
      int extent = phase_range.upper[d] - phase_range.lower[d] + 1;
      idx[d] = phase_range.lower[d] + remaining % extent;
      remaining /= extent;
    }
    dg_calc_sr_vars_pressure_cell(
      up, conf_range, phase_range, use_vmap, vmap, gamma, gamma_inv, u_i, u_i_sq, GammaV, GammaV_sq,
      f, idx, mom_local
    );
  }

  // Reduce each pressure component within the tile.
  for (unsigned int k = 0; k < ncomp; ++k) {
    double sum = BlockReduceT(temp).Sum(mom_local[k]);
    if (threadIdx.x == 0) {
      long segment = (blockIdx.x / num_tiles) * ncomp + k;
      partial_dev[segment * num_tiles + tile] = sum;
    }
    __syncthreads();
  }
}

__global__ static void
dg_calc_sr_vars_pressure_finish_cu_kernel(
  struct gkyl_range conf_range, struct gkyl_array *sr_pressure, const double *partial_dev,
  long conf_start, long num_tiles, unsigned int ncomp
)
{
  typedef cub::BlockReduce<double, SR_PRESSURE_TILE_SIZE, cub::BLOCK_REDUCE_RAKING> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  long segment = blockIdx.x;
  // Combine tile partials for one configuration cell and pressure component.
  double local_sum = 0.0;
  for (long tile = threadIdx.x; tile < num_tiles; tile += blockDim.x) {
    local_sum += partial_dev[segment * num_tiles + tile];
  }
  double sum = BlockReduceT(temp).Sum(local_sum);
  if (threadIdx.x == 0) {
    long conf_local = conf_start + segment / ncomp;
    int idx[GKYL_MAX_DIM];
    gkyl_sub_range_inv_idx(&conf_range, conf_local, idx);
    double *pressure = (double *)gkyl_array_fetch(sr_pressure, gkyl_range_idx(&conf_range, idx));
    pressure[segment % ncomp] = sum;
  }
}
