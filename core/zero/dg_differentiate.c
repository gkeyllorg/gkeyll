#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_dg_differentiate.h>
#include <gkyl_dg_differentiate_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

void
gkyl_dg_differentiate_op_local(const struct gkyl_basis *basis, int dir, int diff_order,
  double dx, int c_oop, struct gkyl_array *out, int c_iop, const struct gkyl_array *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_differentiate_op_local_cu(basis, dir, diff_order, dx, c_oop, out, c_iop, inp);
  }
#endif
  int num_basis = basis->num_basis;
  int ndim = basis->ndim;
  int poly_order = basis->poly_order;

  differentiate_op_t diff_op;
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      diff_op = choose_ser_differentiate_kern(ndim, dir, poly_order, diff_order);
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      diff_op = choose_ten_differentiate_kern(ndim, dir, poly_order, diff_order);
      break;

    default:
      assert(false);
      break;
  }
  assert(diff_op);

  for (size_t i=0; i<out->size; ++i) {
    const double *inp_d = gkyl_array_cfetch(inp, i);
    double *out_d = gkyl_array_fetch(out, i);
    diff_op(dx, inp_d+c_iop*num_basis, out_d+c_oop*num_basis);
  }
}

void
gkyl_dg_differentiate_op_local_range(const struct gkyl_basis *basis, int dir, int diff_order,
  double dx, int c_oop, struct gkyl_array *out, int c_iop, const struct gkyl_array *inp,
  const struct gkyl_range *range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_differentiate_op_local_range_cu(basis, dir, diff_order, dx, c_oop, out, c_iop, inp, range);
  }
#endif
  int num_basis = basis->num_basis;
  int ndim = basis->ndim;
  int poly_order = basis->poly_order;
  differentiate_op_t diff_op;

  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      diff_op = choose_ser_differentiate_kern(ndim, dir, poly_order, diff_order);
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      diff_op = choose_ten_differentiate_kern(ndim, dir, poly_order, diff_order);
      break;

    default:
      assert(false);
      break;
  }
  assert(diff_op);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);
    const double *inp_d = gkyl_array_cfetch(inp, loc);
    double *out_d = gkyl_array_fetch(out, loc);
    diff_op(dx, inp_d+c_iop*num_basis, out_d+c_oop*num_basis);
  }
}
