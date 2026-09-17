#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>

gkyl_dg_bin_op_mem*
gkyl_dg_bin_op_mem_new(size_t nbatch, size_t neqn)
{
  struct gkyl_dg_bin_op_mem *mem = gkyl_malloc(sizeof(struct gkyl_dg_bin_op_mem));

  mem->on_gpu = false;
  mem->batch_sz = nbatch;
  mem->ncols = mem->nrows = neqn;
  
  mem->As = gkyl_nmat_new(nbatch, neqn, neqn);
  mem->xs = gkyl_nmat_new(nbatch, neqn, 1);
  mem->lu_mem = gkyl_nmat_linsolve_lu_new(mem->As->num, mem->As->nr);

  return mem;
}

gkyl_dg_bin_op_mem*
gkyl_dg_bin_op_mem_cu_dev_new(size_t nbatch, size_t neqn)
{
  struct gkyl_dg_bin_op_mem *mem = gkyl_malloc(sizeof(struct gkyl_dg_bin_op_mem));

  mem->on_gpu = false;
  mem->batch_sz = nbatch;
  mem->ncols = mem->nrows = neqn;
  
  mem->As = gkyl_nmat_cu_dev_new(nbatch, neqn, neqn);
  mem->xs = gkyl_nmat_cu_dev_new(nbatch, neqn, 1);
  mem->lu_mem = gkyl_nmat_linsolve_lu_cu_dev_new(mem->As->num, mem->As->nr);

  return mem;
}

void
gkyl_dg_bin_op_mem_release(gkyl_dg_bin_op_mem *mem)
{
  gkyl_nmat_release(mem->As);
  gkyl_nmat_release(mem->xs);
  gkyl_nmat_linsolve_lu_release(mem->lu_mem);
  
  if (mem->on_gpu)
    gkyl_cu_free(mem);
  else
    gkyl_free(mem);
}

// multiplication
void
gkyl_dg_mul_op(const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_mul_op_cu(basis, c_oop, out, c_lop, lop, c_rop, rop);
  }
#endif

  int num_basis = basis->num_basis;
  int ndim = basis->ndim;
  int poly_order = basis->poly_order;
  mul_op_t mul_op;
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mul_op = choose_ser_mul_kern(ndim, poly_order);
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      mul_op = choose_ten_mul_kern(ndim, poly_order);
      
      break;

    default:
      assert(false);
      break;    
  }

  for (size_t i=0; i<out->size; ++i) {
    
    const double *lop_d = gkyl_array_cfetch(lop, i);
    const double *rop_d = gkyl_array_cfetch(rop, i);
    double *out_d = gkyl_array_fetch(out, i);

    mul_op(lop_d+c_lop*num_basis, rop_d+c_rop*num_basis, out_d+c_oop*num_basis);
  }
}

void gkyl_dg_mul_op_range(const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, const struct gkyl_range *range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_mul_op_range_cu(basis, c_oop, out, c_lop, lop, c_rop, rop, range);
  }
#endif

  int num_basis = basis->num_basis;
  int ndim = basis->ndim;
  int poly_order = basis->poly_order;
  mul_op_t mul_op;
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mul_op = choose_ser_mul_kern(ndim, poly_order);
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      mul_op = choose_ten_mul_kern(ndim, poly_order);
      break;

    default:
      assert(false);
      break;    
  }
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *lop_d = gkyl_array_cfetch(lop, loc);
    const double *rop_d = gkyl_array_cfetch(rop, loc);
    double *out_d = gkyl_array_fetch(out, loc);

    mul_op(lop_d+c_lop*num_basis, rop_d+c_rop*num_basis, out_d+c_oop*num_basis);
  }
}

// Dot product.
void
gkyl_dg_dot_product_op(const struct gkyl_basis *basis,
  struct gkyl_array* out,
  const struct gkyl_array* lop,
  const struct gkyl_array* rop)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_dot_product_op_cu(basis, out, lop, rop);
  }
#endif

  int ndim = basis->ndim;
  int poly_order = basis->poly_order;
  mul_op_t mul_op;
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mul_op = choose_ser_mul_kern(ndim, poly_order);

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      mul_op = choose_ten_mul_kern(ndim, poly_order);
      
      break;

    default:
      assert(false);
      break;    
  }

  int num_basis = basis->num_basis;
  int vcomp = lop->ncomp/out->ncomp;

  for (size_t i=0; i<out->size; ++i) {
    
    const double *lop_d = gkyl_array_cfetch(lop, i);
    const double *rop_d = gkyl_array_cfetch(rop, i);
    double *out_d = gkyl_array_fetch(out, i);
    for (int k=0; k<num_basis; k++) out_d[k] = 0.;

    for (int d=0; d<vcomp; d++) {
      double comp_out[num_basis];
      mul_op(lop_d+d*num_basis, rop_d+d*num_basis, comp_out);
      for (int k=0; k<num_basis; k++) out_d[k] += comp_out[k]; 
    }
  }
}

void gkyl_dg_dot_product_op_range(const struct gkyl_basis *basis,
  struct gkyl_array* out,
  const struct gkyl_array* lop,
  const struct gkyl_array* rop, const struct gkyl_range *range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_dot_product_op_range_cu(basis, out, lop, rop, range);
  }
#endif

  int ndim = basis->ndim;
  int poly_order = basis->poly_order;
  mul_op_t mul_op;
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mul_op = choose_ser_mul_kern(ndim, poly_order);

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      mul_op = choose_ten_mul_kern(ndim, poly_order);
      
      break;

    default:
      assert(false);
      break;    
  }

  int num_basis = basis->num_basis;
  int vcomp = lop->ncomp/out->ncomp;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *lop_d = gkyl_array_cfetch(lop, loc);
    const double *rop_d = gkyl_array_cfetch(rop, loc);
    double *out_d = gkyl_array_fetch(out, loc);
    for (int k=0; k<num_basis; k++) out_d[k] = 0.;

    for (int d=0; d<vcomp; d++) {
      double comp_out[num_basis];
      mul_op(lop_d+d*num_basis, rop_d+d*num_basis, comp_out);
      for (int k=0; k<num_basis; k++) out_d[k] += comp_out[k]; 
    }
  }
}

// conf*phase multiplication.
void gkyl_dg_mul_conf_phase_op_range(const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, struct gkyl_array* pout,
  const struct gkyl_array* cop, const struct gkyl_array* pop,
  const struct gkyl_range *crange, const struct gkyl_range *prange)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pout)) {
    return gkyl_dg_mul_conf_phase_op_range_cu(cbasis, pbasis, pout, cop, pop, crange, prange);
  }
#endif

  int cnum_basis = cbasis->num_basis, pnum_basis = pbasis->num_basis;
  assert(pnum_basis > cnum_basis);

  int cdim = cbasis->ndim;
  int vdim = pbasis->ndim - cdim;
  int poly_order = cbasis->poly_order;
  mul_op_t mul_op = choose_mul_conf_phase_kern(pbasis->b_type, cdim, vdim, poly_order);

  struct gkyl_range_iter piter;
  gkyl_range_iter_init(&piter, prange);

  while (gkyl_range_iter_next(&piter)) {
    long ploc = gkyl_range_idx(prange, piter.idx);

    const double *pop_d = gkyl_array_cfetch(pop, ploc);
    double *pout_d = gkyl_array_fetch(pout, ploc);

    int cidx[3]; 
    for (int d=0; d<cdim; d++) cidx[d] = piter.idx[d];
    long cloc = gkyl_range_idx(crange, cidx);
    const double *cop_d = gkyl_array_cfetch(cop, cloc);

    mul_op(cop_d, pop_d, pout_d);
  }
}

// conf*phase multiplication with accumulation to output.
void gkyl_dg_mul_conf_phase_op_accumulate_range(const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, struct gkyl_array* pout, double a, 
  const struct gkyl_array* cop, const struct gkyl_array* pop,
  const struct gkyl_range *crange, const struct gkyl_range *prange)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pout)) {
    return gkyl_dg_mul_conf_phase_op_accumulate_range_cu(cbasis, pbasis, 
      pout, a, cop, pop, crange, prange);
  }
#endif

  int cnum_basis = cbasis->num_basis, pnum_basis = pbasis->num_basis;
  assert(pnum_basis > cnum_basis);

  int cdim = cbasis->ndim;
  int vdim = pbasis->ndim - cdim;
  int poly_order = cbasis->poly_order;
  mul_accumulate_op_t mul_accumulate_op = choose_mul_conf_phase_accumulate_kern(pbasis->b_type, cdim, vdim, poly_order);

  struct gkyl_range_iter piter;
  gkyl_range_iter_init(&piter, prange);

  while (gkyl_range_iter_next(&piter)) {
    long ploc = gkyl_range_idx(prange, piter.idx);

    const double *pop_d = gkyl_array_cfetch(pop, ploc);
    double *pout_d = gkyl_array_fetch(pout, ploc);

    int cidx[3]; 
    for (int d=0; d<cdim; d++) cidx[d] = piter.idx[d];
    long cloc = gkyl_range_idx(crange, cidx);
    const double *cop_d = gkyl_array_cfetch(cop, cloc);

    mul_accumulate_op(a, cop_d, pop_d, pout_d);
  }
}

// division
void
gkyl_dg_div_op(gkyl_dg_bin_op_mem *mem, const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_div_op_cu(mem, basis, c_oop, out, c_lop, lop, c_rop, rop);
  }
#endif

  int num_basis = basis->num_basis;
  int ndim = basis->ndim;
  int poly_order = basis->poly_order;
  div_set_op_t div_set_op;
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      div_set_op = choose_ser_div_set_kern(ndim, poly_order);

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      div_set_op = choose_ten_div_set_kern(ndim, poly_order);
      
      break;

    default:
      assert(false);
      break;    
  }

  struct gkyl_nmat *As = mem->As;
  struct gkyl_nmat *xs = mem->xs;

  for (size_t i=0; i<out->size; ++i) {
    
    const double *lop_d = gkyl_array_cfetch(lop, i);
    const double *rop_d = gkyl_array_cfetch(rop, i);

    struct gkyl_mat A = gkyl_nmat_get(As, i);
    struct gkyl_mat x = gkyl_nmat_get(xs, i);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0);
    div_set_op(&A, &x, lop_d+c_lop*num_basis, rop_d+c_rop*num_basis);
  }

  bool status = gkyl_nmat_linsolve_lu_pa(mem->lu_mem, As, xs);
  assert(status);

  for (size_t i=0; i<out->size; ++i) {
    double *out_d = gkyl_array_fetch(out, i);
    struct gkyl_mat x = gkyl_nmat_get(xs, i);
    binop_div_copy_sol(&x, out_d+c_oop*num_basis);
  }
}

void gkyl_dg_div_op_range(gkyl_dg_bin_op_mem *mem, const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, const struct gkyl_range *range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_div_op_range_cu(mem, basis, c_oop, out, c_lop, lop, c_rop, rop, range);
  }
#endif

  int num_basis = basis->num_basis;
  int ndim = basis->ndim;
  int poly_order = basis->poly_order;
  div_set_op_t div_set_op;
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      div_set_op = choose_ser_div_set_kern(ndim, poly_order);

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      div_set_op = choose_ten_div_set_kern(ndim, poly_order);
      
      break;

    default:
      assert(false);
      break;    
  }

  // allocate memory for use in kernels
  struct gkyl_nmat *As = mem->As;
  struct gkyl_nmat *xs = mem->xs;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *lop_d = gkyl_array_cfetch(lop, loc);
    const double *rop_d = gkyl_array_cfetch(rop, loc);

    struct gkyl_mat A = gkyl_nmat_get(As, count);
    struct gkyl_mat x = gkyl_nmat_get(xs, count);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0); 

    div_set_op(&A, &x, lop_d+c_lop*num_basis, rop_d+c_rop*num_basis);

    count += 1;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(mem->lu_mem, As, xs);
  assert(status);

  gkyl_range_iter_init(&iter, range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    double *out_d = gkyl_array_fetch(out, loc);
    struct gkyl_mat x = gkyl_nmat_get(xs, count);
    binop_div_copy_sol(&x, out_d+c_oop*num_basis);

    count += 1;
  }
}

void gkyl_dg_inv_op(const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out, int c_iop, const struct gkyl_array* iop)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_inv_op_cu(basis, c_oop, out, c_iop, iop);
  }
#endif

  int num_basis = basis->num_basis;
  int ndim = basis->ndim;
  int poly_order = basis->poly_order;
  inv_op_t inv_op;
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      inv_op = choose_ser_inv_kern(ndim, poly_order);
      break;

    default:
      assert(false);
      break;    
  }
  assert(inv_op);

  for (size_t i=0; i<out->size; ++i) {
    const double *iop_d = gkyl_array_cfetch(iop, i);
    double *out_d = gkyl_array_fetch(out, i);

    inv_op(iop_d+c_iop*num_basis, out_d+c_oop*num_basis);
  }
}

void gkyl_dg_inv_op_range(const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out, int c_iop, const struct gkyl_array* iop,
  const struct gkyl_range *range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_inv_op_range_cu(basis, c_oop, out, c_iop, iop, range);
  }
#endif

  int num_basis = basis->num_basis;
  int ndim = basis->ndim;
  int poly_order = basis->poly_order;
  inv_op_t inv_op;
  switch (basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      inv_op = choose_ser_inv_kern(ndim, poly_order);
      break;

    default:
      assert(false);
      break;    
  }
  assert(inv_op);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *iop_d = gkyl_array_cfetch(iop, loc);
    double *out_d = gkyl_array_fetch(out, loc);

    inv_op(iop_d+c_iop*num_basis, out_d+c_oop*num_basis);
  }
}

void
gkyl_dg_calc_op_range(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out,
  int c_iop, const struct gkyl_array *iop,
  struct gkyl_range range, enum gkyl_dg_op op)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {
    return gkyl_dg_calc_op_range_cu(basis, c_oop, out, c_iop, iop, range, op);
  }
#endif
  
  int num_basis = basis->num_basis;
  int ndim = basis->ndim;
  int poly_order = basis->poly_order;

  dp_op_t op_func = dg_get_op_func(op);
  double fact = // factor for rescaling return value of op_func
    op == GKYL_DG_OP_MEAN ? sqrt(pow(2,ndim)) : pow(2,ndim);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    const double *iop_d = gkyl_array_cfetch(iop, loc);
    double *out_d = gkyl_array_fetch(out, loc);

    out_d[c_oop] =
      op_func(num_basis, iop_d+c_iop*num_basis)/fact;
  }  
}

void
gkyl_dg_calc_average_range(const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_iop, const struct gkyl_array* iop, struct gkyl_range range)
{
  gkyl_dg_calc_op_range(basis, c_oop, out, c_iop, iop, range, GKYL_DG_OP_MEAN);
}

void
gkyl_dg_calc_l2_range(const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_iop, const struct gkyl_array* iop, struct gkyl_range range)
{
  gkyl_dg_calc_op_range(basis, c_oop, out, c_iop, iop, range, GKYL_DG_OP_MEAN_L2);
}