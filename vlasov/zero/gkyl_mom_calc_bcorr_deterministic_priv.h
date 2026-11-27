#pragma once

// Form a fixed-order boundary-correction sum for each velocity-space edge.
// Partial tiles are combined before the next edge is added to the output.

#include <cub/cub.cuh>

enum { MOM_CALC_BCORR_TILE_SIZE = 128 };

__global__ static void
mom_calc_bcorr_partial_cu_ker(
  const struct gkyl_mom_calc_bcorr *bcorr, struct gkyl_range conf_rng, struct gkyl_range vel_rng,
  enum gkyl_vel_edge edge, const struct gkyl_array *fin, double *partial_dev, long conf_start,
  long num_tiles, long num_vel, unsigned int ncomp
)
{
  typedef cub::BlockReduce<double, MOM_CALC_BCORR_TILE_SIZE, cub::BLOCK_REDUCE_RAKING> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  long conf_local = conf_start + blockIdx.x / num_tiles;
  long tile = blockIdx.x % num_tiles;
  long vel_local = tile * MOM_CALC_BCORR_TILE_SIZE + threadIdx.x;
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];
  // Reconstruct the phase-space cell on this velocity-space edge.
  gkyl_sub_range_inv_idx(&conf_rng, conf_local, cidx);
  for (int d = 0; d < conf_rng.ndim; ++d) {
    pidx[d] = cidx[d];
  }

  double mom_local[96] = {0.0};
  if (vel_local < num_vel) {
    long remaining = vel_local;
    for (int d = conf_rng.ndim; d < vel_rng.ndim; ++d) {
      int extent = vel_rng.upper[d] - vel_rng.lower[d] + 1;
      pidx[d] = vel_rng.lower[d] + remaining % extent;
      remaining /= extent;
    }
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&bcorr->grid, pidx, xc);
    long phase_index = gkyl_range_idx(&vel_rng, pidx);
    const double *fptr = (const double *)gkyl_array_cfetch(fin, phase_index);
    bcorr->momt->kernel(bcorr->momt, xc, bcorr->grid.dx, pidx, fptr, mom_local, edge);
  }

  // Store a partial edge contribution for each moment component.
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
mom_calc_bcorr_finish_cu_ker(
  struct gkyl_range conf_rng, struct gkyl_array *out, const double *partial_dev, long conf_start,
  long num_tiles, unsigned int ncomp
)
{
  typedef cub::BlockReduce<double, MOM_CALC_BCORR_TILE_SIZE, cub::BLOCK_REDUCE_RAKING> BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  long segment = blockIdx.x;
  // Combine tiles in a fixed order, then add this edge to the output.
  double local_sum = 0.0;
  for (long tile = threadIdx.x; tile < num_tiles; tile += blockDim.x) {
    local_sum += partial_dev[segment * num_tiles + tile];
  }
  double sum = BlockReduceT(temp).Sum(local_sum);
  if (threadIdx.x == 0) {
    long conf_local = conf_start + segment / ncomp;
    int cidx[GKYL_MAX_CDIM];
    gkyl_sub_range_inv_idx(&conf_rng, conf_local, cidx);
    double *mptr = (double *)gkyl_array_fetch(out, gkyl_range_idx(&conf_rng, cidx));
    mptr[segment % ncomp] += sum;
  }
}

static void
mom_calc_bcorr_advance_edge_cu(
  const struct gkyl_mom_calc_bcorr *bcorr, const struct gkyl_range *conf_rng,
  const struct gkyl_range *vel_rng, enum gkyl_vel_edge edge, const struct gkyl_array *fin,
  struct gkyl_array *out
)
{
  long num_vel = 1;
  for (int d = conf_rng->ndim; d < vel_rng->ndim; ++d) {
    num_vel *= vel_rng->upper[d] - vel_rng->lower[d] + 1;
  }
  if (conf_rng->volume == 0 || num_vel == 0) {
    return;
  }
  long num_tiles = (num_vel + MOM_CALC_BCORR_TILE_SIZE - 1) / MOM_CALC_BCORR_TILE_SIZE;
  unsigned int ncomp = out->ncomp;
  // Cap the partial buffer by processing configuration cells in batches.
  const long max_partial_values = (64L * 1024 * 1024) / sizeof(double);
  long batch_conf = max_partial_values / (num_tiles * ncomp);
  if (batch_conf < 1) {
    batch_conf = 1;
  }
  if (batch_conf > conf_rng->volume) {
    batch_conf = conf_rng->volume;
  }
  gkyl_mom_calc_bcorr_buffer *buffer = bcorr->buffer;
  size_t required_bytes = batch_conf * num_tiles * ncomp * sizeof(double);
  if (required_bytes > buffer->partial_bytes) {
    if (buffer->partial_dev) {
      gkyl_cu_free(buffer->partial_dev);
    }
    buffer->partial_dev = (double *)gkyl_cu_malloc(required_bytes);
    buffer->partial_bytes = required_bytes;
  }
  // Complete both reduction passes for each batch before reusing its buffer.
  for (long start = 0; start < conf_rng->volume; start += batch_conf) {
    long count = conf_rng->volume - start < batch_conf ? conf_rng->volume - start : batch_conf;
    long segments = count * ncomp;
    mom_calc_bcorr_partial_cu_ker<<<count * num_tiles, MOM_CALC_BCORR_TILE_SIZE>>>(
      bcorr->on_dev, *conf_rng, *vel_rng, edge, fin->on_dev, buffer->partial_dev, start, num_tiles,
      num_vel, ncomp
    );
    mom_calc_bcorr_finish_cu_ker<<<segments, MOM_CALC_BCORR_TILE_SIZE>>>(
      *conf_rng, out->on_dev, buffer->partial_dev, start, num_tiles, ncomp
    );
  }
}
