/* -*- c++ -*- */

extern "C" {
#include <gkyl_prim_cross_m0deltas_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_prim_cross_m0deltas_set_op_range_cu_kernel(struct gkyl_nmat *As, struct gkyl_nmat *xs,
  struct gkyl_basis basis, bool normNu, double betap1T2,
  double massself, struct gkyl_array* m0self, struct gkyl_array* nuself,
  double massother, struct gkyl_array* m0other, struct gkyl_array* nuother,
  struct gkyl_range range, struct gkyl_array* out)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  div_set_op_t div_set_op = choose_ser_div_set_kern(ndim, poly_order);
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

  int idx[GKYL_MAX_DIM];
  // MF 2022/11/19: Hardcoded to a max number of basis for 3x p2 ser.
  double denom[20], numer[20];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&range, idx);

    const double *m0self_d = (const double*) gkyl_array_cfetch(m0self, start);
    const double *nuself_d = (const double*) gkyl_array_cfetch(nuself, start);
    const double *m0other_d = (const double*) gkyl_array_cfetch(m0other, start);
    const double *nuother_d = (const double*) gkyl_array_cfetch(nuother, start);

    if (nuself_d[0] > 0.0 && nuother_d[0] > 0.0) {
      // compute the numerator and denominator, if collision frequency is
      // constant, in
      //   n_s*delta_s*(beta+1) = 2*(beta+1) * n_s * m_r * n_r * nu_rs / (m_s * n_s * nu_sr + m_r * n_r * nu_rs)
      // or, if the collision frequency varies in space and time, in:
      //   nu_sr*n_s*delta_s*(beta+1) = 2*(beta+1) * n_s * nu_sr * m_r * n_r * nu_rs / (m_s * n_s * nu_sr + m_r * n_r * nu_rs)
      mul_op(nuself_d, m0self_d, denom);
      mul_op(nuother_d, m0other_d, numer);
  
      for (int k=0; k<num_basis; k++) {
        denom[k] *= massself;
        numer[k] *= betap1T2*massother;
      }
  
      array_acc1(num_basis, denom, 1.0/betap1T2, numer);
  
      mul_op(m0self_d, numer, numer);
      if (normNu)
        mul_op(nuself_d, numer, numer);
    }
    else {
      // Both collision frequencies are zero, so set the numerator and
      // denominator to 1. In this case the collision operator will be turned
      // off anyway, so we just want to avoid a division by 0 here.
      denom[0] = 1.0;
      numer[0] = 1.0;
      for (int k=1; k<num_basis; k++) {
        denom[k] = 0.0;
        numer[k] = 0.0;
      }
    }

    struct gkyl_mat A = gkyl_nmat_get(As, linc1);
    struct gkyl_mat x = gkyl_nmat_get(xs, linc1);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0);

    div_set_op(&A, &x, numer, denom);
  }
}

// Modeled after gkyl_dg_div_copy_sol_op_range_cu_kernel in dg_bin_ops_cu.cu.
__global__ void
gkyl_prim_cross_m0deltas_copy_sol_range_cu_kernel(struct gkyl_nmat *xs,
  struct gkyl_basis basis,
  struct gkyl_array* out, struct gkyl_range range)
{
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
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&range, idx);

    double *out_d = (double*) gkyl_array_fetch(out, start);

    struct gkyl_mat x = gkyl_nmat_get(xs, linc1);

    binop_div_copy_sol(&x, out_d);
  }
}

void gkyl_prim_cross_m0deltas_advance_cu(gkyl_prim_cross_m0deltas *up,
  double massself, const struct gkyl_array* m0self, const struct gkyl_array* nuself,
  double massother, const struct gkyl_array* m0other, const struct gkyl_array* nuother,
  struct gkyl_array* out)
{
  int nblocks = up->range->nblocks;
  int nthreads = up->range->nthreads;
  // Allocate memory for use in kernels.
  struct gkyl_nmat *A_d = up->mem->As;
  struct gkyl_nmat *x_d = up->mem->xs;

  // Construct matrices using CUDA kernel.
  gkyl_prim_cross_m0deltas_set_op_range_cu_kernel<<<nblocks, nthreads>>>(A_d->on_dev,
    x_d->on_dev, *up->basis, up->normNu, up->betap1T2, massself, m0self->on_dev, nuself->on_dev,
    massother, m0other->on_dev, nuother->on_dev, *up->range, out->on_dev);

  // Invert all matrices in batch mode.
  bool status = gkyl_nmat_linsolve_lu_pa(up->mem->lu_mem, A_d, x_d);
  assert(status);

  // Copy solution into array (also lives on the device).
  gkyl_prim_cross_m0deltas_copy_sol_range_cu_kernel<<<nblocks, nthreads>>>(x_d->on_dev,
    *up->basis, out->on_dev, *up->range);
}
