/* -*- c++ -*- */

#include <gkyl_array_reduce_cu_priv.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_reduce_priv.h>
#include <gkyl_basis.h>
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayMax_blockRedAtomic_cub(const struct gkyl_array *inp, double *out)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double *)inp->data;

  for (size_t k = 0; k < nComp; ++k) {
    double f = -DBL_MAX;
    if (linc < nCells) {
      f = inp_d[linc * nComp + k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_MAX>(temp, f, &out[k]);
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayMax_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_range range, double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    const double *fptr = NULL;
    if (linc < nCells) {
      gkyl_sub_range_inv_idx(&range, linc, idx);
      long start = gkyl_range_idx(&range, idx);
      fptr = (const double *)gkyl_array_cfetch(inp, start);
    }
    double f = -DBL_MAX;
    if (linc < nCells) {
      f = fptr[k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_MAX>(temp, f, &out[k]);
  }
}

void
gkyl_array_reduce_max_cu(double *out_d, const struct gkyl_array *inp)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, -DBL_MAX
  );
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  if (inp->size > 0) {
    arrayMax_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_range_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_range *range
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, -DBL_MAX
  );
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  if (range->volume > 0) {
    arrayMax_range_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, *range, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayMin_blockRedAtomic_cub(const struct gkyl_array *inp, double *out)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double *)inp->data;

  for (size_t k = 0; k < nComp; ++k) {
    double f = DBL_MAX;
    if (linc < nCells) {
      f = inp_d[linc * nComp + k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_MIN>(temp, f, &out[k]);
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayMin_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_range range, double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    const double *fptr = NULL;
    if (linc < nCells) {
      gkyl_sub_range_inv_idx(&range, linc, idx);
      long start = gkyl_range_idx(&range, idx);
      fptr = (const double *)gkyl_array_cfetch(inp, start);
    }
    double f = DBL_MAX;
    if (linc < nCells) {
      f = fptr[k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_MIN>(temp, f, &out[k]);
  }
}

void
gkyl_array_reduce_min_cu(double *out_d, const struct gkyl_array *inp)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, DBL_MAX
  );
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  if (inp->size > 0) {
    arrayMin_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_range_min_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_range *range
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, DBL_MAX
  );
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  if (range->volume > 0) {
    arrayMin_range_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, *range, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

#ifndef GKYL_DETERMINISTIC_REDUCTIONS_GPU
template <unsigned int BLOCKSIZE>
__global__ void
arraySum_blockRedAtomic_cub(const struct gkyl_array *inp, double *out)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double *)inp->data;

  for (size_t k = 0; k < nComp; ++k) {
    double f = 0;
    if (linc < nCells) {
      f = inp_d[linc * nComp + k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_SUM>(temp, f, &out[k]);
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arraySum_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_range range, double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    const double *fptr = NULL;
    if (linc < nCells) {
      gkyl_sub_range_inv_idx(&range, linc, idx);
      long start = gkyl_range_idx(&range, idx);
      fptr = (const double *)gkyl_array_cfetch(inp, start);
    }
    double f = 0;
    if (linc < nCells) {
      f = fptr[k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_SUM>(temp, f, &out[k]);
  }
}
#endif

#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
#include <gkyl_array_reduce_deterministic_priv.h>
#endif

void
gkyl_array_reduce_sum_cu(double *out_d, const struct gkyl_array *inp)
{
  if (inp->ncomp == 0) {
    return;
  }
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  array_sum_deterministic_cu<false>(out_d, inp, nullptr);
#else
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  if (inp->size > 0) {
    arraySum_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
#endif
}

void
gkyl_array_reduce_range_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_range *range
)
{
  if (inp->ncomp == 0) {
    return;
  }
#ifdef GKYL_DETERMINISTIC_REDUCTIONS_GPU
  array_sum_deterministic_cu<true>(out_d, inp, range);
#else
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  if (range->volume > 0) {
    arraySum_range_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, *range, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
#endif
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayAbsMax_blockRedAtomic_cub(const struct gkyl_array *inp, double *out)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double *)inp->data;

  for (size_t k = 0; k < nComp; ++k) {
    double f = 0;
    if (linc < nCells) {
      f = fabs(inp_d[linc * nComp + k]);
    }
    array_reduce_block<BLOCKSIZE, GKYL_MAX>(temp, f, &out[k]);
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayAbsMax_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_range range, double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    const double *fptr = NULL;
    if (linc < nCells) {
      gkyl_sub_range_inv_idx(&range, linc, idx);
      long start = gkyl_range_idx(&range, idx);
      fptr = (const double *)gkyl_array_cfetch(inp, start);
    }
    double f = 0;
    if (linc < nCells) {
      f = fabs(fptr[k]);
    }
    array_reduce_block<BLOCKSIZE, GKYL_MAX>(temp, f, &out[k]);
  }
}

void
gkyl_array_reduce_abs_max_cu(double *out_d, const struct gkyl_array *inp)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  if (inp->size > 0) {
    arrayAbsMax_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_range_abs_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_range *range
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  if (range->volume > 0) {
    arrayAbsMax_range_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, *range, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE>
__global__ void
arraySqSum_blockRedAtomic_cub(const struct gkyl_array *inp, double *out)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double *)inp->data;

  for (size_t k = 0; k < nComp; ++k) {
    double f = 0;
    if (linc < nCells) {
      f = pow(inp_d[linc * nComp + k], 2);
    }
    array_reduce_block<BLOCKSIZE, GKYL_SUM>(temp, f, &out[k]);
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arraySqSum_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_range range, double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    const double *fptr = NULL;
    if (linc < nCells) {
      gkyl_sub_range_inv_idx(&range, linc, idx);
      long start = gkyl_range_idx(&range, idx);
      fptr = (const double *)gkyl_array_cfetch(inp, start);
    }
    double f = 0;
    if (linc < nCells) {
      f = pow(fptr[k], 2);
    }
    array_reduce_block<BLOCKSIZE, GKYL_SUM>(temp, f, &out[k]);
  }
}

void
gkyl_array_reduce_sq_sum_cu(double *out_d, const struct gkyl_array *inp)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  if (inp->size > 0) {
    arraySqSum_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_range_sq_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_range *range
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  if (range->volume > 0) {
    arraySqSum_range_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, *range, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayMax_weighted_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double *)inp->data;
  const double *wgt_d = (const double *)wgt->data;

  for (size_t k = 0; k < nComp; ++k) {
    double f = -DBL_MAX;
    if (linc < nCells) {
      f = wgt_d[linc * nComp + k] * inp_d[linc * nComp + k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_MAX>(temp, f, &out[k]);
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayMax_weighted_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, const struct gkyl_range range,
  double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    const double *fptr = NULL;
    const double *wptr = NULL;
    if (linc < nCells) {
      gkyl_sub_range_inv_idx(&range, linc, idx);
      long start = gkyl_range_idx(&range, idx);
      fptr = (const double *)gkyl_array_cfetch(inp, start);
      wptr = (const double *)gkyl_array_cfetch(wgt, start);
    }
    double f = -DBL_MAX;
    if (linc < nCells) {
      f = wptr[k] * fptr[k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_MAX>(temp, f, &out[k]);
  }
}

void
gkyl_array_reduce_weighted_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, -DBL_MAX
  );
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  if (inp->size > 0) {
    arrayMax_weighted_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_weighted_range_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, -DBL_MAX
  );
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  if (range->volume > 0) {
    arrayMax_weighted_range_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, *range, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayMin_weighted_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double *)inp->data;
  const double *wgt_d = (const double *)wgt->data;

  for (size_t k = 0; k < nComp; ++k) {
    double f = DBL_MAX;
    if (linc < nCells) {
      f = wgt_d[linc * nComp + k] * inp_d[linc * nComp + k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_MIN>(temp, f, &out[k]);
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayMin_weighted_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, const struct gkyl_range range,
  double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    const double *fptr = NULL;
    const double *wptr = NULL;
    if (linc < nCells) {
      gkyl_sub_range_inv_idx(&range, linc, idx);
      long start = gkyl_range_idx(&range, idx);
      fptr = (const double *)gkyl_array_cfetch(inp, start);
      wptr = (const double *)gkyl_array_cfetch(wgt, start);
    }
    double f = DBL_MAX;
    if (linc < nCells) {
      f = wptr[k] * fptr[k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_MIN>(temp, f, &out[k]);
  }
}

void
gkyl_array_reduce_weighted_min_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, DBL_MAX
  );
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  if (inp->size > 0) {
    arrayMin_weighted_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_weighted_range_min_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, DBL_MAX
  );
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  if (range->volume > 0) {
    arrayMin_weighted_range_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, *range, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE>
__global__ void
arraySum_weighted_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double *)inp->data;
  const double *wgt_d = (const double *)wgt->data;

  for (size_t k = 0; k < nComp; ++k) {
    double f = 0;
    if (linc < nCells) {
      f = wgt_d[linc * nComp + k] * inp_d[linc * nComp + k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_SUM>(temp, f, &out[k]);
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arraySum_weighted_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, const struct gkyl_range range,
  double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    const double *fptr = NULL;
    const double *wptr = NULL;
    if (linc < nCells) {
      gkyl_sub_range_inv_idx(&range, linc, idx);
      long start = gkyl_range_idx(&range, idx);
      fptr = (const double *)gkyl_array_cfetch(inp, start);
      wptr = (const double *)gkyl_array_cfetch(wgt, start);
    }
    double f = 0;
    if (linc < nCells) {
      f = wptr[k] * fptr[k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_SUM>(temp, f, &out[k]);
  }
}

void
gkyl_array_reduce_weighted_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  if (inp->size > 0) {
    arraySum_weighted_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_weighted_range_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  if (range->volume > 0) {
    arraySum_weighted_range_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, *range, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayAbsMax_weighted_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double *)inp->data;
  const double *wgt_d = (const double *)wgt->data;

  for (size_t k = 0; k < nComp; ++k) {
    double f = 0;
    if (linc < nCells) {
      f = fabs(wgt_d[linc * nComp + k] * inp_d[linc * nComp + k]);
    }
    array_reduce_block<BLOCKSIZE, GKYL_MAX>(temp, f, &out[k]);
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayAbsMax_weighted_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, const struct gkyl_range range,
  double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    const double *fptr = NULL;
    const double *wptr = NULL;
    if (linc < nCells) {
      gkyl_sub_range_inv_idx(&range, linc, idx);
      long start = gkyl_range_idx(&range, idx);
      fptr = (const double *)gkyl_array_cfetch(inp, start);
      wptr = (const double *)gkyl_array_cfetch(wgt, start);
    }
    double f = 0;
    if (linc < nCells) {
      f = fabs(wptr[k] * fptr[k]);
    }
    array_reduce_block<BLOCKSIZE, GKYL_MAX>(temp, f, &out[k]);
  }
}

void
gkyl_array_reduce_weighted_abs_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  if (inp->size > 0) {
    arrayAbsMax_weighted_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_weighted_range_abs_max_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  if (range->volume > 0) {
    arrayAbsMax_weighted_range_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, *range, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE>
__global__ void
arraySqSum_weighted_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double *)inp->data;
  const double *wgt_d = (const double *)wgt->data;

  for (size_t k = 0; k < nComp; ++k) {
    double f = 0;
    if (linc < nCells) {
      double value = inp_d[linc * nComp + k];
      f = (wgt_d[linc * nComp + k] * value) * value;
    }
    array_reduce_block<BLOCKSIZE, GKYL_SUM>(temp, f, &out[k]);
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arraySqSum_weighted_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, const struct gkyl_range range,
  double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    const double *fptr = NULL;
    const double *wptr = NULL;
    if (linc < nCells) {
      gkyl_sub_range_inv_idx(&range, linc, idx);
      long start = gkyl_range_idx(&range, idx);
      fptr = (const double *)gkyl_array_cfetch(inp, start);
      wptr = (const double *)gkyl_array_cfetch(wgt, start);
    }
    double f = 0;
    if (linc < nCells) {
      f = (wptr[k] * fptr[k]) * fptr[k];
    }
    array_reduce_block<BLOCKSIZE, GKYL_SUM>(temp, f, &out[k]);
  }
}

void
gkyl_array_reduce_weighted_sq_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  if (inp->size > 0) {
    arraySqSum_weighted_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_weighted_range_sq_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  if (range->volume > 0) {
    arraySqSum_weighted_range_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, *range, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayCellMeanSqSum_weighted_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double *)inp->data;
  const double *wgt_d = (const double *)wgt->data;

  double mean_sq = 0.0;
  if (linc < nCells) {
    for (size_t k = 0; k < nComp; ++k) {
      mean_sq += pow(inp_d[linc * nComp + k], 2);
    }
    mean_sq /= nComp;
  }

  for (size_t k = 0; k < nComp; ++k) {
    double f = 0;
    if (linc < nCells) {
      f = wgt_d[linc * nComp + k] * mean_sq;
    }
    array_reduce_block<BLOCKSIZE, GKYL_SUM>(temp, f, &out[k]);
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayCellMeanSqSum_weighted_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, const struct gkyl_array *wgt, const struct gkyl_range range,
  double *out
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];
  const double *fptr = NULL;
  const double *wptr = NULL;
  if (linc < nCells) {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    fptr = (const double *)gkyl_array_cfetch(inp, start);
    wptr = (const double *)gkyl_array_cfetch(wgt, start);
  }

  double mean_sq = 0;
  if (linc < nCells) {
    for (size_t k = 0; k < nComp; ++k) {
      mean_sq += pow(fptr[k], 2);
    }
    mean_sq /= nComp;
  }

  for (size_t k = 0; k < nComp; ++k) {
    double f = 0;
    if (linc < nCells) {
      f = wptr[k] * mean_sq;
    }
    array_reduce_block<BLOCKSIZE, GKYL_SUM>(temp, f, &out[k]);
  }
}

void
gkyl_array_reduce_weighted_cell_mean_sq_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  if (inp->size > 0) {
    arrayCellMeanSqSum_weighted_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_weighted_range_cell_mean_sq_sum_cu(
  double *out_d, const struct gkyl_array *inp, const struct gkyl_array *wgt,
  const struct gkyl_range *range
)
{
  if (inp->ncomp == 0) {
    return;
  }
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<gkyl_int_div_up(inp->ncomp, nthreads), nthreads>>>(
    out_d, inp->ncomp, 0.0
  );
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  if (range->volume > 0) {
    arrayCellMeanSqSum_weighted_range_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, wgt->on_dev, *range, out_d);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}
