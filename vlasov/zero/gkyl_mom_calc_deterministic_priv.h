#pragma once

// Reduce velocity cells in fixed-size tiles, then combine tile sums for each
// configuration cell and moment coefficient in a fixed order.

#include <cub/cub.cuh>

enum { MOM_CALC_TILE_SIZE = 128 };

__global__ static void
mom_calc_partial_cu_ker(
  const gkyl_mom_calc *mcalc, struct gkyl_range phase_range, struct gkyl_range conf_range,
  const struct gkyl_array *fin, double *partial_dev, long conf_start, long num_tiles, long num_vel,
  unsigned int ncomp
)
{
  typedef cub::BlockReduce<double, MOM_CALC_TILE_SIZE, cub::BLOCK_REDUCE_RAKING> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  long conf_local = conf_start + blockIdx.x / num_tiles;
  long tile = blockIdx.x % num_tiles;
  long vel_local = tile * MOM_CALC_TILE_SIZE + threadIdx.x;
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];
  // Map this block and thread to one configuration cell and velocity tile.
  gkyl_sub_range_inv_idx(&conf_range, conf_local, cidx);
  for (int d = 0; d < conf_range.ndim; ++d) {
    pidx[d] = cidx[d];
  }

  double mom_local[96] = {0.0};
  if (vel_local < num_vel) {
    long remaining = vel_local;
    for (int d = conf_range.ndim; d < phase_range.ndim; ++d) {
      int extent = phase_range.upper[d] - phase_range.lower[d] + 1;
      pidx[d] = phase_range.lower[d] + remaining % extent;
      remaining /= extent;
    }
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&mcalc->grid, pidx, xc);
    long phase_index = gkyl_range_idx(&phase_range, pidx);
    const double *fptr = (const double *)gkyl_array_cfetch(fin, phase_index);
    mcalc->momt->kernel(mcalc->momt, xc, mcalc->grid.dx, pidx, fptr, mom_local, GKYL_VEL_EDGE_NONE);
  }

  // Write one partial per tile and moment component.
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
mom_calc_finish_cu_ker(
  struct gkyl_range conf_range, struct gkyl_array *mout, const double *partial_dev, long conf_start,
  long num_tiles, unsigned int ncomp
)
{
  typedef cub::BlockReduce<double, MOM_CALC_TILE_SIZE, cub::BLOCK_REDUCE_RAKING> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  long segment = blockIdx.x;
  // Reduce the ordered tile partials for one output component.
  double local_sum = 0.0;
  for (long tile = threadIdx.x; tile < num_tiles; tile += blockDim.x) {
    local_sum += partial_dev[segment * num_tiles + tile];
  }
  double sum = BlockReduceT(temp).Sum(local_sum);
  if (threadIdx.x == 0) {
    long conf_local = conf_start + segment / ncomp;
    int cidx[GKYL_MAX_CDIM];
    gkyl_sub_range_inv_idx(&conf_range, conf_local, cidx);
    double *mptr = (double *)gkyl_array_fetch(mout, gkyl_range_idx(&conf_range, cidx));
    mptr[segment % ncomp] = sum;
  }
}

template <typename T>
static void
mom_calc_reserve(T **buffer_dev, size_t *capacity, size_t required)
{
  // Grow the reusable device buffer only when a larger batch requires it.
  if (required > *capacity) {
    if (*buffer_dev) {
      gkyl_cu_free(*buffer_dev);
    }
    *buffer_dev = (T *)gkyl_cu_malloc(required);
    *capacity = required;
  }
}
