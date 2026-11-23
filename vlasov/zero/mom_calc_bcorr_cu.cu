/* -*- c++ -*- */

#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
#include <cub/cub.cuh>
#endif

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_calc_bcorr_priv.h>
#include <gkyl_mom_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
}

#ifndef GKYL_DETERMINISTIC_REDUCTIONS_GPU
__global__ static void
gkyl_mom_calc_bcorr_advance_cu_ker(
  const struct gkyl_mom_calc_bcorr *bcorr, const struct gkyl_range conf_rng,
  struct gkyl_range vel_rng, enum gkyl_vel_edge edge, const struct gkyl_array *fin,
  struct gkyl_array *out
)
{
  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < vel_rng.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&vel_rng, tid, pidx);
    gkyl_rect_grid_cell_center(&bcorr->grid, pidx, xc);

    long lincP = gkyl_range_idx(&vel_rng, pidx);
    const double *fptr = (const double *)gkyl_array_cfetch(fin, lincP);
    double momLocal[96]; // hard-coded to max confBasis.num_basis (3x p=3 Ser) for now.
    for (unsigned int k = 0; k < 96; ++k) {
      momLocal[k] = 0.0;
    }

    // reduce local f to local mom
    bcorr->momt->kernel(bcorr->momt, xc, bcorr->grid.dx, pidx, fptr, &momLocal[0], edge);

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_rng.ndim; k++) {
      cidx[k] = pidx[k];
    }
    long lincC = gkyl_range_idx(&conf_rng, cidx);

    double *mptr = (double *)gkyl_array_fetch(out, lincC);
    for (unsigned int k = 0; k < out->ncomp; ++k) {
      atomicAdd(&mptr[k], momLocal[k]);
    }
  }
}
#else
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
#endif

void
gkyl_mom_calc_bcorr_advance_cu(
  const struct gkyl_mom_calc_bcorr *bcorr, const struct gkyl_range *phase_rng,
  const struct gkyl_range *conf_rng, const struct gkyl_array *GKYL_RESTRICT fin,
  struct gkyl_array *GKYL_RESTRICT out
)
{
  struct gkyl_range vel_rng;
#ifndef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  int nblocks, nthreads;
#endif
  int vlower_idx[GKYL_MAX_DIM], vupper_idx[GKYL_MAX_DIM] = {0};
  for (int dim = 0; dim < phase_rng->ndim; ++dim) {
    vlower_idx[dim] = phase_rng->lower[dim];
    vupper_idx[dim] = phase_rng->upper[dim];
  }
  enum gkyl_vel_edge edge;

  gkyl_array_clear_range(out, 0.0, conf_rng);

  for (int d = 0; d < phase_rng->ndim - conf_rng->ndim; ++d) {
    edge = gkyl_vel_edge(d + GKYL_MAX_CDIM);
    vlower_idx[conf_rng->ndim + d] = phase_rng->upper[conf_rng->ndim + d];
    vupper_idx[conf_rng->ndim + d] = phase_rng->upper[conf_rng->ndim + d];
    gkyl_sub_range_init(&vel_rng, phase_rng, vlower_idx, vupper_idx);
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
    mom_calc_bcorr_advance_edge_cu(bcorr, conf_rng, &vel_rng, edge, fin, out);
#else
    nblocks = vel_rng.nblocks;
    nthreads = vel_rng.nthreads;

    gkyl_mom_calc_bcorr_advance_cu_ker<<<nblocks, nthreads>>>(
      bcorr->on_dev, *conf_rng, vel_rng, edge, fin->on_dev, out->on_dev
    );
#endif

    edge = gkyl_vel_edge(d);
    vlower_idx[conf_rng->ndim + d] = phase_rng->lower[conf_rng->ndim + d];
    vupper_idx[conf_rng->ndim + d] = phase_rng->lower[conf_rng->ndim + d];
    gkyl_sub_range_init(&vel_rng, phase_rng, vlower_idx, vupper_idx);
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
    mom_calc_bcorr_advance_edge_cu(bcorr, conf_rng, &vel_rng, edge, fin, out);
#else
    nblocks = vel_rng.nblocks;
    nthreads = vel_rng.nthreads;

    gkyl_mom_calc_bcorr_advance_cu_ker<<<nblocks, nthreads>>>(
      bcorr->on_dev, *conf_rng, vel_rng, edge, fin->on_dev, out->on_dev
    );
#endif

    // Reset indices for loop over each velocity dimension
    vlower_idx[conf_rng->ndim + d] = phase_rng->lower[conf_rng->ndim + d];
    vupper_idx[conf_rng->ndim + d] = phase_rng->upper[conf_rng->ndim + d];
  }
}

gkyl_mom_calc_bcorr *
gkyl_mom_calc_bcorr_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_mom_type *momt)
{
  gkyl_mom_calc_bcorr *up = (gkyl_mom_calc_bcorr *)gkyl_malloc(sizeof(gkyl_mom_calc_bcorr));
  up->grid = *grid;

  struct gkyl_mom_type *mt = gkyl_mom_type_acquire(momt);
  up->momt = mt->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  up->buffer = (gkyl_mom_calc_bcorr_buffer *)gkyl_malloc(sizeof(gkyl_mom_calc_bcorr_buffer));
  *up->buffer = {};
#else
  up->buffer = 0;
#endif

  gkyl_mom_calc_bcorr *up_cu = (gkyl_mom_calc_bcorr *)gkyl_cu_malloc(sizeof(gkyl_mom_calc_bcorr));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_mom_calc_bcorr), GKYL_CU_MEMCPY_H2D);

  up->momt = mt;
  up->on_dev = up_cu;

  return up;
}
