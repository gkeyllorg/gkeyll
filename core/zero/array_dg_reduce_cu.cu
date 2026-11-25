/* -*- c++ -*- */

#include <gkyl_array_reduce_cu_priv.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_array_dg_reduce.h>
#include <gkyl_array_dg_reduce_priv.h>
#include <gkyl_basis.h>
}

template <unsigned int BLOCKSIZE>
__global__ void
dg_arrayMax_blockRedAtomic_cub(
  const struct gkyl_array *inp, double *out, int comp, const struct gkyl_basis *basis
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  const double *inp_d = (const double *)inp->data;

  double f = -DBL_MAX;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k = 0; k < num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[linc * nComp + comp * basis->num_basis], arr_nodal, k);
      f = fmax(f, arr_nodal[k]);
    }
  }
  array_reduce_block<BLOCKSIZE, GKYL_MAX, false>(temp, f, &out[0]);
}

template <unsigned int BLOCKSIZE>
__global__ void
dg_arrayMax_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, double *out, int comp, const struct gkyl_basis *basis,
  struct gkyl_range range
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  int idx[GKYL_MAX_DIM];
  const double *inp_d = NULL;
  if (linc < nCells) {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    inp_d = (const double *)gkyl_array_cfetch(inp, start);
  }

  double f = -DBL_MAX;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k = 0; k < num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[comp * basis->num_basis], arr_nodal, k);
      f = fmax(f, arr_nodal[k]);
    }
  }
  array_reduce_block<BLOCKSIZE, GKYL_MAX, false>(temp, f, &out[0]);
}

void
gkyl_array_dg_reducec_max_cu(
  double *out_d, const struct gkyl_array *inp, int comp, const struct gkyl_basis *basis
)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<1, 1>>>(out_d, 1, -DBL_MAX);
  if (inp->size > 0) {
    int nblocks = gkyl_int_div_up(inp->size, nthreads);
    dg_arrayMax_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_dg_reducec_range_max_cu(
  double *out_d, const struct gkyl_array *inp, int comp, const struct gkyl_basis *basis,
  const struct gkyl_range *range
)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<1, 1>>>(out_d, 1, -DBL_MAX);
  if (range->volume > 0) {
    int nblocks = gkyl_int_div_up(range->volume, nthreads);
    dg_arrayMax_range_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis, *range);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE>
__global__ void
dg_arrayMin_blockRedAtomic_cub(
  const struct gkyl_array *inp, double *out, int comp, const struct gkyl_basis *basis
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  const double *inp_d = (const double *)inp->data;

  double f = DBL_MAX;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k = 0; k < num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[linc * nComp + comp * basis->num_basis], arr_nodal, k);
      f = fmin(f, arr_nodal[k]);
    }
  }
  array_reduce_block<BLOCKSIZE, GKYL_MIN, false>(temp, f, &out[0]);
}

template <unsigned int BLOCKSIZE>
__global__ void
dg_arrayMin_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, double *out, int comp, const struct gkyl_basis *basis,
  struct gkyl_range range
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  int idx[GKYL_MAX_DIM];
  const double *inp_d = NULL;
  if (linc < nCells) {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    inp_d = (const double *)gkyl_array_cfetch(inp, start);
  }

  double f = DBL_MAX;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k = 0; k < num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[comp * basis->num_basis], arr_nodal, k);
      f = fmin(f, arr_nodal[k]);
    }
  }
  array_reduce_block<BLOCKSIZE, GKYL_MIN, false>(temp, f, &out[0]);
}

void
gkyl_array_dg_reducec_min_cu(
  double *out_d, const struct gkyl_array *inp, int comp, const struct gkyl_basis *basis
)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<1, 1>>>(out_d, 1, DBL_MAX);
  if (inp->size > 0) {
    int nblocks = gkyl_int_div_up(inp->size, nthreads);
    dg_arrayMin_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_dg_reducec_range_min_cu(
  double *out_d, const struct gkyl_array *inp, int comp, const struct gkyl_basis *basis,
  const struct gkyl_range *range
)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<1, 1>>>(out_d, 1, DBL_MAX);
  if (range->volume > 0) {
    int nblocks = gkyl_int_div_up(range->volume, nthreads);
    dg_arrayMin_range_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis, *range);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

template <unsigned int BLOCKSIZE>
__global__ void
dg_arraySum_blockRedAtomic_cub(
  const struct gkyl_array *inp, double *out, int comp, const struct gkyl_basis *basis
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = inp->size;
  size_t nComp = inp->ncomp;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  const double *inp_d = (const double *)inp->data;

  double f = 0;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k = 0; k < num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[linc * nComp + comp * basis->num_basis], arr_nodal, k);
      f += arr_nodal[k];
    }
  }
  array_reduce_block<BLOCKSIZE, GKYL_SUM, false>(temp, f, &out[0]);
}

template <unsigned int BLOCKSIZE>
__global__ void
dg_arraySum_range_blockRedAtomic_cub(
  const struct gkyl_array *inp, double *out, int comp, const struct gkyl_basis *basis,
  struct gkyl_range range
)
{
  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume;
  int num_nodes = basis->num_quad;
  const int num_nodes_max = 27; // MF 2025/01/15: hard coded to p=2 3x for now.

  int idx[GKYL_MAX_DIM];
  const double *inp_d = NULL;
  if (linc < nCells) {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    inp_d = (const double *)gkyl_array_cfetch(inp, start);
  }

  double f = 0;
  if (linc < nCells) {
    double arr_nodal[num_nodes_max];
    for (int k = 0; k < num_nodes; ++k) {
      basis->modal_to_quad_nodal(&inp_d[comp * basis->num_basis], arr_nodal, k);
      f += arr_nodal[k];
    }
  }
  array_reduce_block<BLOCKSIZE, GKYL_SUM, false>(temp, f, &out[0]);
}

void
gkyl_array_dg_reducec_sum_cu(
  double *out_d, const struct gkyl_array *inp, int comp, const struct gkyl_basis *basis
)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<1, 1>>>(out_d, 1, 0.0);
  if (inp->size > 0) {
    int nblocks = gkyl_int_div_up(inp->size, nthreads);
    dg_arraySum_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}

void
gkyl_array_dg_reducec_range_sum_cu(
  double *out_d, const struct gkyl_array *inp, int comp, const struct gkyl_basis *basis,
  const struct gkyl_range *range
)
{
  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  array_reduce_init_cu_kernel<<<1, 1>>>(out_d, 1, 0.0);
  if (range->volume > 0) {
    int nblocks = gkyl_int_div_up(range->volume, nthreads);
    dg_arraySum_range_blockRedAtomic_cub<nthreads>
      <<<nblocks, nthreads>>>(inp->on_dev, out_d, comp, basis, *range);
  }
  // device synchronize required because out_d may be host pinned memory
  cudaDeviceSynchronize();
}
