#include <gkyl_prim_cross_m0deltas.h>
#include <gkyl_prim_cross_m0deltas_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_alloc.h>

gkyl_prim_cross_m0deltas*
gkyl_prim_cross_m0deltas_new(bool normNu, const struct gkyl_basis *basis,
  const struct gkyl_range *range, double betap1, bool use_gpu)
{
  gkyl_prim_cross_m0deltas *up = gkyl_malloc(sizeof(gkyl_prim_cross_m0deltas));

  // MF 2022/11/19: hardcoded arrays for a max of 3x p2 Ser basis below.
  assert(basis->num_basis <= 20);

  up->normNu = normNu;
  up->basis = basis;
  up->range = range;
  up->betap1T2 = betap1*2.0;
  up->use_gpu = use_gpu;

  // Preallocate memory for the weak division.
  up->mem = use_gpu ? gkyl_dg_bin_op_mem_cu_dev_new(range->volume, basis->num_basis)
	             : gkyl_dg_bin_op_mem_new(range->volume, basis->num_basis);

  return up;
}

void
gkyl_prim_cross_m0deltas_advance(gkyl_prim_cross_m0deltas *up,
  double massself, const struct gkyl_array* m0self, const struct gkyl_array* nuself,
  double massother, const struct gkyl_array* m0other, const struct gkyl_array* nuother,
  struct gkyl_array* out)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    return gkyl_prim_cross_m0deltas_advance_cu(up, massself, m0self, nuself,
      massother, m0other, nuother, out);
#endif

  int num_basis = up->basis->num_basis;
  int ndim = up->basis->ndim;
  int poly_order = up->basis->poly_order;
  div_set_op_t div_set_op = choose_ser_div_set_kern(ndim, poly_order);
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

  // Allocate memory for use in kernels.
  struct gkyl_nmat *As = up->mem->As;
  struct gkyl_nmat *xs = up->mem->xs;
  // MF 2022/11/19: Hardcoded to a max number of basis for 3x p2 ser.
  double denom[20], numer[20];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, up->range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(up->range, iter.idx);

    const double *m0self_d = gkyl_array_cfetch(m0self, loc);
    const double *nuself_d = gkyl_array_cfetch(nuself, loc);
    const double *m0other_d = gkyl_array_cfetch(m0other, loc);
    const double *nuother_d = gkyl_array_cfetch(nuother, loc);

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
        numer[k] *= up->betap1T2*massother;
      }
  
      array_acc1(num_basis, denom, 1.0/up->betap1T2, numer);
  
      mul_op(m0self_d, numer, numer);
      if (up->normNu)
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

    struct gkyl_mat A = gkyl_nmat_get(As, count);
    struct gkyl_mat x = gkyl_nmat_get(xs, count);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0);

    div_set_op(&A, &x, numer, denom);

    count += 1;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem->lu_mem, As, xs);
  assert(status);

  gkyl_range_iter_init(&iter, up->range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(up->range, iter.idx);

    double *out_d = gkyl_array_fetch(out, loc);
    struct gkyl_mat x = gkyl_nmat_get(xs, count);
    binop_div_copy_sol(&x, out_d);

    count += 1;
  }
}

void
gkyl_prim_cross_m0deltas_release(gkyl_prim_cross_m0deltas* up)
{
  gkyl_dg_bin_op_mem_release(up->mem);
  gkyl_free(up);
}
