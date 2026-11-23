/* -*- c++ -*- */
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
#include <cub/cub.cuh>
#endif

extern "C" {
#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_bc_emission_spectrum_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
}

// start ID for use in various loops
#define START_ID (threadIdx.x + blockIdx.x * blockDim.x)

__global__ static void
gkyl_bc_emission_spectrum_set_exterm_params_cu_ker(
  struct gkyl_emission_spectrum_model *spectrum_model,
  struct gkyl_emission_yield_model *yield_model, int cdim, int vdim, double mass_in, double mass_out
)
{
  spectrum_model->cdim = cdim;
  spectrum_model->vdim = vdim;
  spectrum_model->mass = mass_out;

  yield_model->cdim = cdim;
  yield_model->vdim = vdim;
  yield_model->mass = mass_in;
}

#ifndef GKYL_DETERMINISTIC_REDUCTIONS_GPU
__global__ static void
gkyl_bc_emission_spectrum_advance_cu_weight_ker(
  int cdim, int dir, enum gkyl_edge_loc edge, const struct gkyl_array *bflux,
  struct gkyl_array *weight, struct gkyl_rect_grid grid, struct gkyl_array *yield,
  const struct gkyl_range impact_buff_r, const struct gkyl_range impact_cbuff_r
)
{
  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < impact_buff_r.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&impact_buff_r, tid, pidx);

    gkyl_rect_grid_cell_center(&grid, pidx, xc);

    long lincP = gkyl_range_idx(&impact_buff_r, pidx);

    const double *inp = (const double *)gkyl_array_cfetch(bflux, lincP);
    const double *gain = (const double *)gkyl_array_cfetch(yield, lincP);
    double wLocal[2];
    for (unsigned int k = 0; k < weight->ncomp; ++k) {
      wLocal[k] = 0.0;
    }

    bc_weighted_delta(inp, cdim, dir, edge, xc, gain, &wLocal[0]);

    // get conf-space linear index.
    for (unsigned int i = 0; i < impact_cbuff_r.ndim; i++) {
      cidx[i] = pidx[i];
    }
    long lincC = gkyl_range_idx(&impact_cbuff_r, cidx);

    double *wptr = (double *)gkyl_array_fetch(weight, lincC);
    for (unsigned int k = 0; k < weight->ncomp; ++k) {
      atomicAdd(&wptr[k], wLocal[k]);
    }
  }
}
#else
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
#endif

__global__ static void
gkyl_bc_emission_spectrum_advance_cu_normalize_ker(
  const struct gkyl_array *weight, struct gkyl_array *k, const struct gkyl_array *flux,
  const struct gkyl_range impact_cbuff_r, struct gkyl_emission_spectrum_model *spectrum_model
)
{
  int cidx[GKYL_MAX_CDIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < impact_cbuff_r.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&impact_cbuff_r, tid, cidx);
    long lincC = gkyl_range_idx(&impact_cbuff_r, cidx);

    const double *w = (const double *)gkyl_array_cfetch(weight, lincC);
    const double *boundary_flux = (const double *)gkyl_array_cfetch(flux, lincC);
    double *fac = (double *)gkyl_array_fetch(k, lincC);
    // With no incident flux there is nothing to emit (and the mean yield is 0/0).
    double effective_delta = w[1] != 0.0 ? w[0] / w[1] : 0.0;
    spectrum_model->normalization(fac, spectrum_model, boundary_flux, effective_delta);
  }
}

__global__ static void
gkyl_bc_emission_spectrum_advance_cu_accumulate_ker(
  const struct gkyl_array *spectrum, struct gkyl_array *f_emit, const struct gkyl_array *k,
  const struct gkyl_range emit_buff_r, const struct gkyl_range impact_cbuff_r
)
{
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < emit_buff_r.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&emit_buff_r, tid, pidx);

    // get conf-space linear index.
    for (unsigned int i = 0; i < impact_cbuff_r.ndim; i++) {
      cidx[i] = pidx[i];
    }
    long lincC = gkyl_range_idx(&impact_cbuff_r, cidx);

    long lincP = gkyl_range_idx(&emit_buff_r, pidx);

    const double *inp = (const double *)gkyl_array_cfetch(spectrum, lincP);
    const double *fac = (const double *)gkyl_array_cfetch(k, lincC);
    double *out = (double *)gkyl_array_fetch(f_emit, lincP);

    for (int c = 0; c < spectrum->ncomp; ++c) {
      out[c] += fac[0] * inp[c];
    }
  }
}

void
gkyl_bc_emission_spectrum_set_extern_params_cu(
  const struct gkyl_bc_emission_spectrum *up, int cdim, int vdim, double mass_in, double mass_out
)
{
  gkyl_bc_emission_spectrum_set_exterm_params_cu_ker<<<1, 1>>>(
    up->spectrum_model->on_dev, up->yield_model->on_dev, cdim, vdim, mass_in, mass_out
  );
}

void
gkyl_bc_emission_spectrum_advance_cu(
  const struct gkyl_bc_emission_spectrum *up, struct gkyl_range *impact_buff_r,
  struct gkyl_range *impact_cbuff_r, struct gkyl_range *emit_buff_r, struct gkyl_array *bflux,
  struct gkyl_array *f_emit, struct gkyl_array *yield, struct gkyl_array *spectrum,
  struct gkyl_array *weight, struct gkyl_array *flux, struct gkyl_array *k
)
{
#ifndef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  int nblocks = impact_buff_r->nblocks, nthreads = impact_buff_r->nthreads;
#endif

  // The weighted-mean numerator and denominator are accumulated with atomics, so they
  // must start from zero on every call (as on the host); otherwise the effective yield
  // becomes a running average over all previous stages.
  gkyl_array_clear_range(weight, 0.0, impact_cbuff_r);

  // Calculate weighted mean numerator and denominator.
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  bc_emission_spectrum_advance_weight_cu(up, impact_buff_r, impact_cbuff_r, bflux, yield, weight);
#else
  gkyl_bc_emission_spectrum_advance_cu_weight_ker<<<nblocks, nthreads>>>(
    up->cdim, up->dir, up->edge, bflux->on_dev, weight->on_dev, *up->grid, yield->on_dev,
    *impact_buff_r, *impact_cbuff_r
  );
#endif

#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  int nblocks = impact_cbuff_r->nblocks, nthreads = impact_cbuff_r->nthreads;
#else
  nblocks = impact_cbuff_r->nblocks;
  nthreads = impact_cbuff_r->nthreads;
#endif
  gkyl_bc_emission_spectrum_advance_cu_normalize_ker<<<nblocks, nthreads>>>(
    weight->on_dev, k->on_dev, flux->on_dev, *impact_cbuff_r, up->spectrum_model->on_dev
  );

  nblocks = emit_buff_r->nblocks;
  nthreads = emit_buff_r->nthreads;

  // Finish weighted mean calculation and accumulate to buffer
  gkyl_bc_emission_spectrum_advance_cu_accumulate_ker<<<nblocks, nthreads>>>(
    spectrum->on_dev, f_emit->on_dev, k->on_dev, *emit_buff_r, *impact_cbuff_r
  );
}
