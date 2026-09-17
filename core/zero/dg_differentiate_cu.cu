/* -*- c++ -*- */

extern "C" {
#include <gkyl_dg_differentiate.h>
#include <gkyl_dg_differentiate_priv.h>
#include <gkyl_util.h>
}

#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)
// Compute size of 'arr' 
#define NSIZE(arr) (arr->size)

__global__ void
gkyl_dg_differentiate_op_local_cu_kernel(struct gkyl_basis basis, int dir, int diff_order,
  double dx, int c_oop, struct gkyl_array *out, int c_iop, const struct gkyl_array *inp)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;

  differentiate_op_t diff_op;
  switch (basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      diff_op = choose_ser_differentiate_kern(ndim, dir, poly_order, diff_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      diff_op = choose_ten_differentiate_kern(ndim, dir, poly_order, diff_order);
      break;
    default:
      return;
  }

  for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x) {
    const double *inp_d = (const double *) gkyl_array_cfetch(inp, linc);
    double *out_d = (double *) gkyl_array_fetch(out, linc);
    diff_op(dx, inp_d+c_iop*num_basis, out_d+c_oop*num_basis);
  }
}

void
gkyl_dg_differentiate_op_local_cu(const struct gkyl_basis *basis, int dir, int diff_order,
  double dx, int c_oop, struct gkyl_array *out, int c_iop, const struct gkyl_array *inp)
{
  gkyl_dg_differentiate_op_local_cu_kernel<<<out->nblocks, out->nthreads>>>(
    *basis, dir, diff_order, dx, c_oop, out->on_dev, c_iop, inp->on_dev);
}

__global__ void
gkyl_dg_differentiate_op_local_range_cu_kernel(struct gkyl_basis basis, int dir, int diff_order,
  double dx, int c_oop, struct gkyl_array *out, int c_iop, const struct gkyl_array *inp,
  struct gkyl_range range)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;

  differentiate_op_t diff_op;
  switch (basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      diff_op = choose_ser_differentiate_kern(ndim, dir, poly_order, diff_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      diff_op = choose_ten_differentiate_kern(ndim, dir, poly_order, diff_order);
      break;
    default:
      return;
  }

  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    long start = gkyl_range_idx(&range, idx);

    const double *inp_d = (const double *) gkyl_array_cfetch(inp, start);
    double *out_d = (double *) gkyl_array_fetch(out, start);
    diff_op(dx, inp_d+c_iop*num_basis, out_d+c_oop*num_basis);
  }
}

void
gkyl_dg_differentiate_op_local_range_cu(const struct gkyl_basis *basis, int dir, int diff_order,
  double dx, int c_oop, struct gkyl_array *out, int c_iop, const struct gkyl_array *inp,
  const struct gkyl_range *range)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_dg_differentiate_op_local_range_cu_kernel<<<nblocks, nthreads>>>(
    *basis, dir, diff_order, dx, c_oop, out->on_dev, c_iop, inp->on_dev, *range);
}
