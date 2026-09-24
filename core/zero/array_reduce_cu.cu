/* -*- c++ -*- */

// CUB for reductions.
#include <cub/cub.cuh>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_array_reduce.h>
#include <gkyl_array_reduce_priv.h>
#include <gkyl_basis.h>
}

// CUDA does not natively support atomics for MAX and MIN on doubles (addition/sum is fine).
// These functions utilize the atomicCAS (compare and swap) to thread-by-thread find if the
// input value is greater than (for max) or less than (for min) the output and swap the values
// if the condition is satisfied, along the way doing the double_as_longlong and longlong_as_double
// conversions needed to determine if the double (as a long long) is indeed greater than or less than
// the output. Note that because this operation is done thread-by-thread we still use CUB to perform
// the reduction over CUDA blocks, but then the threads are compared thread-by-thread. 
// These particular functions are adapted from (adapted by JJ on 03/14/24): 
// https://github.com/treecode/Bonsai/blob/master/runtime/profiling/derived_atomic_functions.h
__device__ static __forceinline__ double 
atomicMax_double(double *address, double val)
{
  unsigned long long int ret = __double_as_longlong(*address);
  while(val > __longlong_as_double(ret))
  {
    unsigned long long int old = ret;
    if((ret = atomicCAS((unsigned long long int*)address, old, __double_as_longlong(val))) == old)
      break;
  }
  return __longlong_as_double(ret);
}

__device__ static __forceinline__ double 
atomicMin_double(double *address, double val)
{
  unsigned long long int ret = __double_as_longlong(*address);
  while(val < __longlong_as_double(ret))
  {
    unsigned long long int old = ret;
    if((ret = atomicCAS((unsigned long long int*)address, old, __double_as_longlong(val))) == old)
      break;
  }
  return __longlong_as_double(ret);
}

template <unsigned int BLOCKSIZE> 
__global__ void
arrayMax_blockRedAtomic_cub(const struct gkyl_array* inp, double* out)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double*) inp->data;

  for (size_t k = 0; k < nComp; ++k) {
    out[k] = -DBL_MAX;
    double f = -DBL_MAX;
    if (linc < nCells) f = inp_d[linc*nComp+k];
    double bResult = 0;
    bResult = BlockReduceT(temp).Reduce(f,
#if CUDART_VERSION > 12090
    ::cuda::maximum()
#else
    cub::Max()
#endif
    );
    if (threadIdx.x < BLOCKSIZE) {
      atomicMax_double(&out[k], bResult);
    }
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayMax_range_blockRedAtomic_cub(const struct gkyl_array* inp, const struct gkyl_range range, double* out)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    out[k] = -DBL_MAX;
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    const double* fptr = (const double*) gkyl_array_cfetch(inp, start);
    double f = -DBL_MAX;
    if (linc < nCells) f = fptr[k];
    double bResult = 0;
    bResult = BlockReduceT(temp).Reduce(f, 
#if CUDART_VERSION > 12090
    ::cuda::maximum()
#else
    cub::Max()
#endif
    );
    if (threadIdx.x < BLOCKSIZE) {
      atomicMax_double(&out[k], bResult);
    }
  }
}

void
gkyl_array_reduce_max_cu(double *out_d, const struct gkyl_array* inp)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;  
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  arrayMax_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_range_max_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  arrayMax_range_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, *range, out_d);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE> 
__global__ void
arrayMin_blockRedAtomic_cub(const struct gkyl_array* inp, double* out)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double*) inp->data;

  for (size_t k = 0; k < nComp; ++k) {
    out[k] = DBL_MAX;
    double f = DBL_MAX;
    if (linc < nCells) f = inp_d[linc*nComp+k];
    double bResult = 0;
    bResult = BlockReduceT(temp).Reduce(f, 
#if CUDART_VERSION > 12090
    ::cuda::minimum()
#else
    cub::Min()
#endif
    );
    if (threadIdx.x < BLOCKSIZE) {
      atomicMin_double(&out[k], bResult);
    }
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayMin_range_blockRedAtomic_cub(const struct gkyl_array* inp, const struct gkyl_range range, double* out)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    out[k] = DBL_MAX;
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    const double* fptr = (const double*) gkyl_array_cfetch(inp, start);
    double f = DBL_MAX;
    if (linc < nCells) f = fptr[k];
    double bResult = 0;
    bResult = BlockReduceT(temp).Reduce(f, 
#if CUDART_VERSION > 12090
    ::cuda::minimum()
#else
    cub::Min()
#endif
    );
    if (threadIdx.x < BLOCKSIZE) { 
      atomicMin_double(&out[k], bResult);
    }
  }
}

void
gkyl_array_reduce_min_cu(double *out_d, const struct gkyl_array* inp)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;  
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  arrayMin_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_range_min_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  arrayMin_range_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, *range, out_d);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE> 
__global__ void
arraySum_blockRedAtomic_cub(const struct gkyl_array* inp, double* out)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;

  const double *inp_d = (const double*) inp->data;

  for (size_t k = 0; k < nComp; ++k) {
    double f = 0;
    if (linc < nCells) f = inp_d[linc*nComp+k];
    double bResult = 0;
    bResult = BlockReduceT(temp).Reduce(f, 
#if CUDART_VERSION > 12090
    ::cuda::std::plus()
#else
    cub::Sum()
#endif
    );
    if (threadIdx.x == 0) {
      atomicAdd(&out[k], bResult);
    }
  }
}

template <unsigned int BLOCKSIZE>
__global__ void
arraySum_range_blockRedAtomic_cub(const struct gkyl_array* inp, const struct gkyl_range range, double* out)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  size_t nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (size_t k = 0; k < nComp; ++k) {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    const double* fptr = (const double*) gkyl_array_cfetch(inp, start);
    double f = 0;
    if (linc < nCells) f = fptr[k];
    double bResult = 0;
    bResult = BlockReduceT(temp).Reduce(f, 
#if CUDART_VERSION > 12090
    ::cuda::std::plus()
#else
    cub::Sum()
#endif
    );
    if (threadIdx.x == 0) {
      atomicAdd(&out[k], bResult);
    }
  }
}

void
gkyl_array_reduce_sum_cu(double *out_d, const struct gkyl_array* inp)
{
  gkyl_cu_memset(out_d, 0, inp->ncomp*sizeof(double));
  
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;  
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  arraySum_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_range_sum_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range *range)
{
  gkyl_cu_memset(out_d, 0, inp->ncomp*sizeof(double));
  
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  arraySum_range_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, *range, out_d);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

