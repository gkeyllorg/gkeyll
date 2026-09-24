/* -*- c++ -*- */

// CUB for reductions.
#include <cub/cub.cuh>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_integrate_priv.h>
}

__global__ static void
gkyl_array_integrate_set_ker_cu(struct gkyl_array_integrate *up, enum gkyl_array_integrate_op op, struct gkyl_basis basis)
{
  int ndim = basis.ndim, poly_order = basis.poly_order;

  if (op == GKYL_ARRAY_INTEGRATE_OP_NONE) {
    up->kernel = gkyl_array_integrate_none_ker_list_ser[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_ABS) {
    up->kernel = gkyl_array_integrate_abs_ker_list_ser[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_SQ) {
    up->kernel = gkyl_array_integrate_sq_ker_list_ser[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_SQ_WEIGHTED) {
    if (basis.b_type == GKYL_BASIS_MODAL_SERENDIPITY)
      up->kernel = gkyl_array_integrate_sq_weighted_ker_list_ser[ndim-1].kernels[poly_order-1];
    else if (basis.b_type == GKYL_BASIS_MODAL_GKHYBRID)
      up->kernel = gkyl_array_integrate_sq_weighted_ker_list_gkhyb[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_GRAD_SQ) {
    up->kernel = gkyl_array_integrate_gradsq_ker_list[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_GRADPERP_SQ) {
    up->kernel = gkyl_array_integrate_gradperpsq_ker_list[ndim-1].kernels[poly_order-1];
  }
  else if (op == GKYL_ARRAY_INTEGRATE_OP_EPS_GRADPERP_SQ) {
    up->kernel = gkyl_array_integrate_epsgradperpsq_ker_list[ndim-1].kernels[poly_order-1];
  } else {
    assert(false);
  }
}

struct gkyl_array_integrate*
gkyl_array_integrate_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  int num_comp, enum gkyl_array_integrate_op op)
{
  // Allocate space for new updater.
  struct gkyl_array_integrate *up = (struct gkyl_array_integrate*) gkyl_malloc(sizeof(struct gkyl_array_integrate));

  up->op = op;
  up->num_basis = basis->num_basis;
  up->num_comp = num_comp;
  up->use_gpu = true;
  for (int d=0; d<grid->ndim; ++d) up->dxSq[d] = grid->dx[d]*grid->dx[d];

  assert(basis->poly_order > 0); // Need to check normalization for p=0.

  int ndim = basis->ndim;
  up->vol = 1.0;
  for (unsigned d=0; d<ndim; ++d)
    up->vol *= grid->dx[d]/2.0;

  // Copy struct to device.
  struct gkyl_array_integrate *up_cu = (struct gkyl_array_integrate*) gkyl_cu_malloc(sizeof(struct gkyl_array_integrate));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_array_integrate), GKYL_CU_MEMCPY_H2D);

  // Set the kernel.
  gkyl_array_integrate_set_ker_cu<<<1,1>>>(up_cu, op, *basis);

  up->on_dev = up_cu;

  return up;
}

template <unsigned int BLOCKSIZE>
__global__ void
array_integrate_blockRedAtomic_cub(struct gkyl_array_integrate *up, const struct gkyl_array *inp,
  double factor, const struct gkyl_array *weight, const struct gkyl_range range, struct gkyl_range weight_range, double *out)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  int idx[GKYL_MAX_DIM];
  gkyl_sub_range_inv_idx(&range, linc, idx);
  long start = gkyl_range_idx(&range, idx);

  double outLocal[10]; // Set to max of 10 (e.g. heat flux tensor).
  for (unsigned int k=0; k<up->num_comp; ++k)
    outLocal[k] = 0.0;

  // Integrate in this cell
  if (linc < range.volume) {
    const double *fptr = (const double*) gkyl_array_cfetch(inp, start);

    const double *wptr = 0;
    if (weight) {
      int widx[GKYL_MAX_DIM];
      for (int d=0; d<weight_range.ndim; d++) widx[d] = idx[d]; 
      long linidx_w = gkyl_range_idx(&weight_range, widx);
      wptr = (const double*) gkyl_array_cfetch(weight, linidx_w);
    }

    up->kernel(up->dxSq, up->vol*factor, up->num_comp, up->num_basis, wptr, fptr, outLocal);
  }

  for (size_t k = 0; k < up->num_comp; ++k) {
    double bResult = 0;
    bResult = BlockReduceT(temp).Reduce(outLocal[k],
#if CUDART_VERSION > 12090
    ::cuda::std::plus()
#else
    cub::Sum()
#endif
    );
    if (threadIdx.x == 0)
      atomicAdd(&out[k], bResult);
  }
}

void gkyl_array_integrate_advance_cu(gkyl_array_integrate *up, const struct gkyl_array *fin,
  double factor, const struct gkyl_array *weight, const struct gkyl_range *range, const struct gkyl_range *weight_range, double *out)
{
  gkyl_cu_memset(out, 0, up->num_comp*sizeof(double));

  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  struct gkyl_array *weight_on_dev = NULL;
  struct gkyl_range weight_range_copy;
  if (weight) {
    weight_on_dev = weight->on_dev;
    weight_range_copy = *weight_range;
  }
  array_integrate_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(up->on_dev, fin->on_dev, factor, 
    weight_on_dev, *range, weight_range_copy, out);
  // device synchronize required because out may be host pinned memory
  cudaDeviceSynchronize();
}
