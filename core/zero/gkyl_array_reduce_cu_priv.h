#pragma once

#include <cub/cub.cuh>
#include <gkyl_elem_type.h>

// CAS-based extrema for doubles, used once per block after the CUB reduction.
// Adapted from https://github.com/treecode/Bonsai/blob/master/runtime/profiling/derived_atomic_functions.h.
__device__ static __forceinline__ double
array_reduce_atomic_max(double *address, double val)
{
  unsigned long long int ret = __double_as_longlong(*address);
  while (val > __longlong_as_double(ret)) {
    unsigned long long int old = ret;
    if ((ret = atomicCAS((unsigned long long int *)address, old, __double_as_longlong(val))) ==
        old) {
      break;
    }
  }
  return __longlong_as_double(ret);
}

__device__ static __forceinline__ double
array_reduce_atomic_min(double *address, double val)
{
  unsigned long long int ret = __double_as_longlong(*address);
  while (val < __longlong_as_double(ret)) {
    unsigned long long int old = ret;
    if ((ret = atomicCAS((unsigned long long int *)address, old, __double_as_longlong(val))) ==
        old) {
      break;
    }
  }
  return __longlong_as_double(ret);
}

// Initialize outputs in a separate launch before any block contributes.
__global__ static void
array_reduce_init_cu_kernel(double *out, size_t ncomp, double identity)
{
  size_t comp = blockIdx.x * blockDim.x + threadIdx.x;
  if (comp < ncomp) {
    out[comp] = identity;
  }
}

/**
 * Reduce a block and atomically combine its result with an initialized output.
 * All block threads must call this function. With REUSE_STORAGE=true, a final
 * barrier allows temp to be reused on return; otherwise the kernel must finish
 * without reusing temp.
 *
 * @param temp Shared CUB storage belonging to the calling block.
 * @param value Thread contribution, or the reduction identity for inactive cells.
 * @param out Device-accessible output scalar, initialized before the kernel launch.
 */
template <unsigned int BLOCKSIZE, enum gkyl_array_op OP, bool REUSE_STORAGE = true>
__device__ static __forceinline__ void
array_reduce_block(
  typename cub::BlockReduce<double, BLOCKSIZE>::TempStorage &temp, double value, double *out
)
{
  static_assert(OP == GKYL_MIN || OP == GKYL_MAX || OP == GKYL_SUM, "Unsupported block reduction");
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;
  double result;
  if (OP == GKYL_MIN) {
    result = BlockReduceT(temp).Reduce(
      value,
#if CUDART_VERSION > 12090
      ::cuda::minimum()
#else
      cub::Min()
#endif
    );
    if (threadIdx.x == 0) {
      array_reduce_atomic_min(out, result);
    }
  } else if (OP == GKYL_MAX) {
    result = BlockReduceT(temp).Reduce(
      value,
#if CUDART_VERSION > 12090
      ::cuda::maximum()
#else
      cub::Max()
#endif
    );
    if (threadIdx.x == 0) {
      array_reduce_atomic_max(out, result);
    }
  } else {
    result = BlockReduceT(temp).Reduce(
      value,
#if CUDART_VERSION > 12090
      ::cuda::std::plus()
#else
      cub::Sum()
#endif
    );
    if (threadIdx.x == 0) {
      atomicAdd(out, result);
    }
  }
  if (REUSE_STORAGE) {
    __syncthreads();
  }
}
