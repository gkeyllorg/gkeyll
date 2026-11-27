#pragma once

// Sum the weighted-emission numerator and denominator over ordered impact
// velocity tiles before normalizing the emission spectrum.

#include <cub/cub.cuh>

enum { BC_EMISSION_WEIGHT_TILE_SIZE = 128 };

__global__ static void
bc_emission_spectrum_weight_partial_cu_ker(
  int cdim, int dir, enum gkyl_edge_loc edge, const struct gkyl_array *bflux,
  struct gkyl_rect_grid grid, const struct gkyl_array *yield, struct gkyl_range impact_buff_r,
  struct gkyl_range impact_cbuff_r, double *partial_dev, long conf_start, long num_tiles,
  long num_vel
)
{
  typedef cub::BlockReduce<double, BC_EMISSION_WEIGHT_TILE_SIZE, cub::BLOCK_REDUCE_RAKING>
    BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  long conf_local = conf_start + blockIdx.x / num_tiles;
  long tile = blockIdx.x % num_tiles;
  long vel_local = tile * BC_EMISSION_WEIGHT_TILE_SIZE + threadIdx.x;
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];
  // Map this thread to an impact velocity cell at the selected boundary.
  gkyl_sub_range_inv_idx(&impact_cbuff_r, conf_local, cidx);
  for (int d = 0; d < impact_cbuff_r.ndim; ++d) {
    pidx[d] = cidx[d];
  }

  double weight_local[2] = {0.0};
  if (vel_local < num_vel) {
    long remaining = vel_local;
    for (int d = impact_cbuff_r.ndim; d < impact_buff_r.ndim; ++d) {
      int extent = impact_buff_r.upper[d] - impact_buff_r.lower[d] + 1;
      pidx[d] = impact_buff_r.lower[d] + remaining % extent;
      remaining /= extent;
    }
    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&grid, pidx, xc);
    long phase_index = gkyl_range_idx(&impact_buff_r, pidx);
    const double *inp = (const double *)gkyl_array_cfetch(bflux, phase_index);
    const double *gain = (const double *)gkyl_array_cfetch(yield, phase_index);
    bc_weighted_delta(inp, cdim, dir, edge, xc, gain, weight_local);
  }

  // Store separate partial sums for the weighted numerator and denominator.
  for (unsigned int k = 0; k < 2; ++k) {
    double sum = BlockReduceT(temp).Sum(weight_local[k]);
    if (threadIdx.x == 0) {
      long segment = (blockIdx.x / num_tiles) * 2 + k;
      partial_dev[segment * num_tiles + tile] = sum;
    }
    __syncthreads();
  }
}

__global__ static void
bc_emission_spectrum_weight_finish_cu_ker(
  struct gkyl_range impact_cbuff_r, struct gkyl_array *weight, const double *partial_dev,
  long conf_start, long num_tiles
)
{
  typedef cub::BlockReduce<double, BC_EMISSION_WEIGHT_TILE_SIZE, cub::BLOCK_REDUCE_RAKING>
    BlockReduceT;
  __shared__ typename BlockReduceT::TempStorage temp;
  long segment = blockIdx.x;
  // Reduce the tiles for one normalization component.
  double local_sum = 0.0;
  for (long tile = threadIdx.x; tile < num_tiles; tile += blockDim.x) {
    local_sum += partial_dev[segment * num_tiles + tile];
  }
  double sum = BlockReduceT(temp).Sum(local_sum);
  if (threadIdx.x == 0) {
    long conf_local = conf_start + segment / 2;
    int cidx[GKYL_MAX_CDIM];
    gkyl_sub_range_inv_idx(&impact_cbuff_r, conf_local, cidx);
    double *wptr = (double *)gkyl_array_fetch(weight, gkyl_range_idx(&impact_cbuff_r, cidx));
    wptr[segment % 2] = sum;
  }
}

static void
bc_emission_spectrum_advance_weight_cu(
  const struct gkyl_bc_emission_spectrum *up, const struct gkyl_range *impact_buff_r,
  const struct gkyl_range *impact_cbuff_r, const struct gkyl_array *bflux,
  const struct gkyl_array *yield, struct gkyl_array *weight
)
{
  long num_vel = 1;
  for (int d = impact_cbuff_r->ndim; d < impact_buff_r->ndim; ++d) {
    num_vel *= impact_buff_r->upper[d] - impact_buff_r->lower[d] + 1;
  }
  if (impact_cbuff_r->volume == 0 || num_vel == 0) {
    return;
  }
  assert(weight->ncomp == 2);
  long num_tiles = (num_vel + BC_EMISSION_WEIGHT_TILE_SIZE - 1) / BC_EMISSION_WEIGHT_TILE_SIZE;
  // Bound temporary storage by batching configuration cells.
  const long max_partial_values = (64L * 1024 * 1024) / sizeof(double);
  long batch_conf = max_partial_values / (num_tiles * 2);
  if (batch_conf < 1) {
    batch_conf = 1;
  }
  if (batch_conf > impact_cbuff_r->volume) {
    batch_conf = impact_cbuff_r->volume;
  }
  gkyl_bc_emission_spectrum_buffer *buffer = up->buffer;
  size_t required_bytes = batch_conf * num_tiles * 2 * sizeof(double);
  if (required_bytes > buffer->partial_bytes) {
    if (buffer->partial_dev) {
      gkyl_cu_free(buffer->partial_dev);
    }
    buffer->partial_dev = (double *)gkyl_cu_malloc(required_bytes);
    buffer->partial_bytes = required_bytes;
  }
  // Finish each batch before reusing its partial buffer.
  for (long start = 0; start < impact_cbuff_r->volume; start += batch_conf) {
    long count = impact_cbuff_r->volume - start < batch_conf ? impact_cbuff_r->volume - start :
                                                               batch_conf;
    bc_emission_spectrum_weight_partial_cu_ker<<<count * num_tiles, BC_EMISSION_WEIGHT_TILE_SIZE>>>(
      up->cdim, up->dir, up->edge, bflux->on_dev, *up->grid, yield->on_dev, *impact_buff_r,
      *impact_cbuff_r, buffer->partial_dev, start, num_tiles, num_vel
    );
    bc_emission_spectrum_weight_finish_cu_ker<<<count * 2, BC_EMISSION_WEIGHT_TILE_SIZE>>>(
      *impact_cbuff_r, weight->on_dev, buffer->partial_dev, start, num_tiles
    );
  }
}
