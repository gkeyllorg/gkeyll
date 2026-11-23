/* -*- c++ -*- */

#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
#include <cub/cub.cuh>
#endif

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_calc_priv.h>
#include <gkyl_mom_type.h>
#include <gkyl_util.h>
}

#ifndef GKYL_DETERMINISTIC_REDUCTIONS_GPU
__global__ static void
gkyl_mom_calc_advance_cu_ker(
  const gkyl_mom_calc *mcalc, const struct gkyl_range phase_range,
  const struct gkyl_range conf_range, const struct gkyl_array *GKYL_RESTRICT fin,
  struct gkyl_array *GKYL_RESTRICT mout
)
{
  double xc[GKYL_MAX_DIM];
  int pidx[GKYL_MAX_DIM], cidx[GKYL_MAX_CDIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x * blockDim.x; tid < phase_range.volume;
       tid += blockDim.x * gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    gkyl_rect_grid_cell_center(&mcalc->grid, pidx, xc);

    long lincP = gkyl_range_idx(&phase_range, pidx);
    const double *fptr = (const double *)gkyl_array_cfetch(fin, lincP);
    double momLocal[96]; // hard-coded to 3 * max confBasis.num_basis (3x p=3 Ser) for now.
    for (unsigned int k = 0; k < 96; ++k) {
      momLocal[k] = 0.0;
    }

    // reduce local f to local mom
    mcalc->momt->kernel(
      mcalc->momt, xc, mcalc->grid.dx, pidx, fptr, &momLocal[0], GKYL_VEL_EDGE_NONE
    );

    // get conf-space linear index.
    for (unsigned int k = 0; k < conf_range.ndim; k++) {
      cidx[k] = pidx[k];
    }
    long lincC = gkyl_range_idx(&conf_range, cidx);

    double *mptr = (double *)gkyl_array_fetch(mout, lincC);
    for (unsigned int k = 0; k < mout->ncomp; ++k) {
      atomicAdd(&mptr[k], momLocal[k]);
    }
  }
}
#else
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
  if (required > *capacity) {
    if (*buffer_dev) {
      gkyl_cu_free(*buffer_dev);
    }
    *buffer_dev = (T *)gkyl_cu_malloc(required);
    *capacity = required;
  }
}
#endif

void
gkyl_mom_calc_advance_cu(
  const gkyl_mom_calc *mcalc, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_array *GKYL_RESTRICT fin,
  struct gkyl_array *GKYL_RESTRICT mout
)
{
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  long num_vel = 1;
  for (int d = conf_range->ndim; d < phase_range->ndim; ++d) {
    num_vel *= phase_range->upper[d] - phase_range->lower[d] + 1;
  }
  long num_tiles = (num_vel + MOM_CALC_TILE_SIZE - 1) / MOM_CALC_TILE_SIZE;
  unsigned int ncomp = mout->ncomp;
  if (conf_range->volume == 0 || num_vel == 0) {
    return;
  }
  const long max_partial_values = (64L * 1024 * 1024) / sizeof(double);
  long batch_conf = max_partial_values / (num_tiles * ncomp);
  if (batch_conf < 1) {
    batch_conf = 1;
  }
  if (batch_conf > conf_range->volume) {
    batch_conf = conf_range->volume;
  }
  gkyl_mom_calc_buffer *buffer = mcalc->buffer;
  mom_calc_reserve(
    &buffer->partial_dev, &buffer->partial_bytes, batch_conf * num_tiles * ncomp * sizeof(double)
  );
  gkyl_array_clear_range(mout, 0.0, conf_range);
  for (long start = 0; start < conf_range->volume; start += batch_conf) {
    long count = conf_range->volume - start < batch_conf ? conf_range->volume - start : batch_conf;
    long segments = count * ncomp;
    mom_calc_partial_cu_ker<<<count * num_tiles, MOM_CALC_TILE_SIZE>>>(
      mcalc->on_dev, *phase_range, *conf_range, fin->on_dev, buffer->partial_dev, start, num_tiles,
      num_vel, ncomp
    );
    mom_calc_finish_cu_ker<<<segments, MOM_CALC_TILE_SIZE>>>(
      *conf_range, mout->on_dev, buffer->partial_dev, start, num_tiles, ncomp
    );
  }
#else
  int nblocks = phase_range->nblocks, nthreads = phase_range->nthreads;
  gkyl_array_clear_range(mout, 0.0, conf_range);
  gkyl_mom_calc_advance_cu_ker<<<nblocks, nthreads>>>(
    mcalc->on_dev, *phase_range, *conf_range, fin->on_dev, mout->on_dev
  );
#endif
}

gkyl_mom_calc *
gkyl_mom_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_mom_type *momt)
{
  gkyl_mom_calc *up = (gkyl_mom_calc *)gkyl_malloc(sizeof(gkyl_mom_calc));
  up->grid = *grid;

  struct gkyl_mom_type *mt = gkyl_mom_type_acquire(momt);
  up->momt = mt->on_dev; // so memcpy below gets dev copy

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  up->buffer = (gkyl_mom_calc_buffer *)gkyl_malloc(sizeof(gkyl_mom_calc_buffer));
  *up->buffer = {};
#endif

  gkyl_mom_calc *up_cu = (gkyl_mom_calc *)gkyl_cu_malloc(sizeof(gkyl_mom_calc));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_mom_calc), GKYL_CU_MEMCPY_H2D);

  up->momt = mt; // host portion of struct should have host copy
  up->on_dev = up_cu; // host pointer
  up->use_gpu = true;

  return up;
}
