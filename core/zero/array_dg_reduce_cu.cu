/* -*- c++ -*- */

// CUB for reductions.
#include <cub/cub.cuh>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_array_dg_reduce.h>
#include <gkyl_array_dg_reduce_priv.h>
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
dg_arrayMax_blockRedAtomic_cub(const struct gkyl_array* inp, double* out, int comp, const struct gkyl_basis *basis)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  const double *inp_d = (const double*) inp->data;

  out[0] = -DBL_MAX;
  double f = -DBL_MAX;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k=0; k<num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[linc*nComp+comp*basis->num_basis], arr_nodal, k);
      f = fmax(f, arr_nodal[k]);
    }
  }
  double bResult = 0;
  bResult = BlockReduceT(temp).Reduce(f, 
#if CUDART_VERSION > 12090
    ::cuda::maximum()
#else
    cub::Max()
#endif
    );
  if (threadIdx.x < BLOCKSIZE) {
    atomicMax_double(&out[0], bResult);
  }
}

template <unsigned int BLOCKSIZE> 
__global__ void
dg_arrayMax_range_blockRedAtomic_cub(const struct gkyl_array* inp, double* out,
  int comp, const struct gkyl_basis *basis, struct gkyl_range range)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  int idx[GKYL_MAX_DIM];
  gkyl_sub_range_inv_idx(&range, linc, idx);
  long start = gkyl_range_idx(&range, idx);
  const double* inp_d = (const double*) gkyl_array_cfetch(inp, start);

  out[0] = -DBL_MAX;
  double f = -DBL_MAX;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k=0; k<num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[comp*basis->num_basis], arr_nodal, k);
      f = fmax(f, arr_nodal[k]);
    }
  }
  double bResult = 0;
  bResult = BlockReduceT(temp).Reduce(f, 
#if CUDART_VERSION > 12090
    ::cuda::maximum()
#else
    cub::Max()
#endif
    );
  if (threadIdx.x < BLOCKSIZE) {
    atomicMax_double(&out[0], bResult);
  }
}

void
gkyl_array_dg_reducec_max_cu(double *out_d, const struct gkyl_array* inp, int comp, const struct gkyl_basis *basis)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;  
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  dg_arrayMax_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_dg_reducec_range_max_cu(double *out_d, const struct gkyl_array* inp,
  int comp, const struct gkyl_basis *basis, const struct gkyl_range *range)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;  
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  dg_arrayMax_range_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis, *range);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE> 
__global__ void
dg_arrayMin_blockRedAtomic_cub(const struct gkyl_array* inp, double* out, int comp, const struct gkyl_basis *basis)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  const double *inp_d = (const double*) inp->data;

  out[0] = DBL_MAX;
  double f = DBL_MAX;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k=0; k<num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[linc*nComp+comp*basis->num_basis], arr_nodal, k);
      f = fmin(f, arr_nodal[k]);
    }
  }
  double bResult = 0;
  bResult = BlockReduceT(temp).Reduce(f, 
#if CUDART_VERSION > 12090
    ::cuda::minimum()
#else
    cub::Min()
#endif
    );
  if (threadIdx.x < BLOCKSIZE) {
    atomicMin_double(&out[0], bResult);
  }
}

template <unsigned int BLOCKSIZE> 
__global__ void
dg_arrayMin_range_blockRedAtomic_cub(const struct gkyl_array* inp, double* out,
  int comp, const struct gkyl_basis *basis, struct gkyl_range range)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  int idx[GKYL_MAX_DIM];
  gkyl_sub_range_inv_idx(&range, linc, idx);
  long start = gkyl_range_idx(&range, idx);
  const double* inp_d = (const double*) gkyl_array_cfetch(inp, start);

  out[0] = DBL_MAX;
  double f = DBL_MAX;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k=0; k<num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[comp*basis->num_basis], arr_nodal, k);
      f = fmin(f, arr_nodal[k]);
    }
  }
  double bResult = 0;
  bResult = BlockReduceT(temp).Reduce(f, 
#if CUDART_VERSION > 12090
    ::cuda::minimum()
#else
    cub::Min()
#endif
    );
  if (threadIdx.x < BLOCKSIZE) {
    atomicMin_double(&out[0], bResult);
  }
}

void
gkyl_array_dg_reducec_min_cu(double *out_d, const struct gkyl_array* inp, int comp, const struct gkyl_basis *basis)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;  
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  dg_arrayMin_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_dg_reducec_range_min_cu(double *out_d, const struct gkyl_array* inp,
  int comp, const struct gkyl_basis *basis, const struct gkyl_range *range)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;  
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  dg_arrayMin_range_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis, *range);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE> 
__global__ void
dg_arraySum_blockRedAtomic_cub(const struct gkyl_array* inp, double* out, int comp, const struct gkyl_basis *basis)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  const double *inp_d = (const double*) inp->data;

  double f = 0;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k=0; k<num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[linc*nComp+comp*basis->num_basis], arr_nodal, k);
      f += arr_nodal[k];
    }
  }
  double bResult = 0;
  bResult = BlockReduceT(temp).Reduce(f, 
#if CUDART_VERSION > 12090
    ::cuda::std::plus()
#else
    cub::Sum()
#endif
    );
  if (threadIdx.x == 0) {
    atomicAdd(&out[0], bResult);
  }
}

template <unsigned int BLOCKSIZE> 
__global__ void
dg_arraySum_range_blockRedAtomic_cub(const struct gkyl_array* inp, double* out,
int comp, const struct gkyl_basis *basis, struct gkyl_range range)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  int idx[GKYL_MAX_DIM];
  gkyl_sub_range_inv_idx(&range, linc, idx);
  long start = gkyl_range_idx(&range, idx);
  const double* inp_d = (const double*) gkyl_array_cfetch(inp, start);

  double f = 0;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k=0; k<num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[comp*basis->num_basis], arr_nodal, k);
      f += arr_nodal[k];
    }
  }
  double bResult = 0;
  bResult = BlockReduceT(temp).Reduce(f, 
#if CUDART_VERSION > 12090
    ::cuda::std::plus()
#else
    cub::Sum()
#endif
    );
  if (threadIdx.x == 0) {
    atomicAdd(&out[0], bResult);
  }
}

void
gkyl_array_dg_reducec_sum_cu(double *out_d, const struct gkyl_array* inp, int comp, const struct gkyl_basis *basis)
{
  gkyl_cu_memset(out_d, 0, sizeof(double));

  const int nthreads = GKYL_DEFAULT_NUM_THREADS;  
  int nblocks = gkyl_int_div_up(inp->size, nthreads);
  dg_arraySum_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_dg_reducec_range_sum_cu(double *out_d, const struct gkyl_array* inp,
  int comp, const struct gkyl_basis *basis, const struct gkyl_range *range)
{
  gkyl_cu_memset(out_d, 0, sizeof(double));

  const int nthreads = GKYL_DEFAULT_NUM_THREADS;  
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  dg_arraySum_range_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis, *range);
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

