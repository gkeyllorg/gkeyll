/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>
}

// start ID for use in various loops
#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)

__global__ void
gkyl_dg_mul_op_cu_kernel(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

  for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x) {
    
    const double *lop_d = (const double*) gkyl_array_cfetch(lop, linc);
    const double *rop_d = (const double*) gkyl_array_cfetch(rop, linc);
    double *out_d = (double*) gkyl_array_fetch(out, linc);

    mul_op(lop_d+c_lop*num_basis, rop_d+c_rop*num_basis, out_d+c_oop*num_basis);
  }  
}

// Host-side wrapper for dg multiplication operation
void
gkyl_dg_mul_op_cu(const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
  gkyl_dg_mul_op_cu_kernel<<<out->nblocks, out->nthreads>>>(*basis, c_oop, out->on_dev,
    c_lop, lop->on_dev, c_rop, rop->on_dev);
}

__global__ void
gkyl_dg_mul_op_range_cu_kernel(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range range)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

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

    const double *lop_d = (const double*) gkyl_array_cfetch(lop, start);
    const double *rop_d = (const double*) gkyl_array_cfetch(rop, start);
    double *out_d = (double*) gkyl_array_fetch(out, start);

    mul_op(lop_d+c_lop*num_basis, rop_d+c_rop*num_basis, out_d+c_oop*num_basis);
  }
}

// Host-side wrapper for range-based dg multiplication operation
void
gkyl_dg_mul_op_range_cu(const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, const struct gkyl_range *range)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_dg_mul_op_range_cu_kernel<<<nblocks, nthreads>>>(*basis, c_oop, out->on_dev,
    c_lop, lop->on_dev, c_rop, rop->on_dev, *range);
}

__global__ void
gkyl_dg_mul_conf_phase_op_range_cu_kernel(struct gkyl_basis cbasis,
  struct gkyl_basis pbasis, struct gkyl_array* pout,
  const struct gkyl_array* cop, const struct gkyl_array* pop,
  struct gkyl_range crange, struct gkyl_range prange)
{
  int cdim = cbasis.ndim;
  int vdim = pbasis.ndim - cdim;
  int poly_order = cbasis.poly_order;
  mul_op_t mul_op = choose_mul_conf_phase_kern(pbasis.b_type, cdim, vdim, poly_order);

  int pidx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < prange.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&prange, linc1, pidx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&prange, pidx);

    const double *pop_d = (const double*) gkyl_array_cfetch(pop, start);
    double *pout_d = (double*) gkyl_array_fetch(pout, start);

    int cidx[3];
    for (int d=0; d<cdim; d++) cidx[d] = pidx[d];
    long cstart = gkyl_range_idx(&crange, cidx);
    const double *cop_d = (const double*) gkyl_array_cfetch(cop, cstart);

    mul_op(cop_d, pop_d, pout_d);
  }
}

// Host-side wrapper for range-based dg conf*phase multiplication.
void
gkyl_dg_mul_conf_phase_op_range_cu(const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, struct gkyl_array* pout,
  const struct gkyl_array* cop, const struct gkyl_array* pop,
  const struct gkyl_range *crange, const struct gkyl_range *prange)
{
  int nblocks = prange->nblocks;
  int nthreads = prange->nthreads;
  gkyl_dg_mul_conf_phase_op_range_cu_kernel<<<nblocks, nthreads>>>(*cbasis, *pbasis,
    pout->on_dev, cop->on_dev, pop->on_dev, *crange, *prange);
}

static void
gkyl_parallelize_components_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads 
  // so we can parallelize over components too
  dimBlock->y = ncomp; // ncomp *must* be less than 256
  dimGrid->y = 1;
  dimBlock->x = GKYL_DEFAULT_NUM_THREADS/ncomp;
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
}

__global__ void
gkyl_dg_mul_conf_phase_op_accumulate_range_cu_kernel(struct gkyl_basis cbasis,
  struct gkyl_basis pbasis, struct gkyl_array* pout, double a, 
  const struct gkyl_array* cop, const struct gkyl_array* pop,
  struct gkyl_range crange, struct gkyl_range prange)
{
  int cdim = cbasis.ndim;
  int vdim = pbasis.ndim - cdim;
  int poly_order = cbasis.poly_order;
  // On GPU, choose kernels which parallelize over components
  mul_accumulate_comp_par_op_t mul_accumulate_op = choose_mul_conf_phase_accumulate_comp_par_kern(pbasis.b_type, 
    cdim, vdim, poly_order);

  int pidx[GKYL_MAX_DIM];
  long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < prange.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&prange, linc1, pidx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&prange, pidx);

    const double *pop_d = (const double*) gkyl_array_cfetch(pop, start);
    double *pout_d = (double*) gkyl_array_fetch(pout, start);

    int cidx[3];
    for (int d=0; d<cdim; d++) cidx[d] = pidx[d];
    long cstart = gkyl_range_idx(&crange, cidx);
    const double *cop_d = (const double*) gkyl_array_cfetch(cop, cstart);

    mul_accumulate_op(a, cop_d, pop_d, pout_d, linc2);
  }
}

// Host-side wrapper for range-based dg conf*phase multiplication with accumulation to output.
void
gkyl_dg_mul_conf_phase_op_accumulate_range_cu(const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, struct gkyl_array* pout, double a, 
  const struct gkyl_array* cop, const struct gkyl_array* pop,
  const struct gkyl_range *crange, const struct gkyl_range *prange)
{
  dim3 dimGrid, dimBlock;
  int num_phase_basis = pbasis->num_basis;
  gkyl_parallelize_components_kernel_launch_dims(&dimGrid, &dimBlock, *prange, num_phase_basis);
  gkyl_dg_mul_conf_phase_op_accumulate_range_cu_kernel<<<dimGrid, dimBlock>>>(*cbasis, *pbasis,
    pout->on_dev, a, cop->on_dev, pop->on_dev, *crange, *prange);
}

__global__ void
gkyl_dg_dot_product_op_cu_kernel(struct gkyl_basis basis,
  struct gkyl_array* out, const struct gkyl_array* lop,
  const struct gkyl_array* rop)
{
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

  int num_basis = basis.num_basis;
  int vcomp = lop->ncomp/out->ncomp;

  for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x) {
    
    const double *lop_d = (const double*) gkyl_array_cfetch(lop, linc);
    const double *rop_d = (const double*) gkyl_array_cfetch(rop, linc);
    double *out_d = (double*) gkyl_array_fetch(out, linc);
    for (int k=0; k<num_basis; k++) out_d[k] = 0.;

    for (int d=0; d<vcomp; d++) {
      double comp_out[20];  // MF 2022/09/08: Hardcoded to number of basis in 3x p=2.
      mul_op(lop_d+d*num_basis, rop_d+d*num_basis, comp_out);
      for (int k=0; k<num_basis; k++) out_d[k] += comp_out[k];
    }
  }  
}

// Host-side wrapper for dg dot product operation.
void
gkyl_dg_dot_product_op_cu(const struct gkyl_basis *basis,
  struct gkyl_array* out, const struct gkyl_array* lop,
  const struct gkyl_array* rop)
{
  assert(basis->num_basis <= 20); // MF 2022/09/08: see hardcode in kernel above.
  gkyl_dg_dot_product_op_cu_kernel<<<out->nblocks, out->nthreads>>>(*basis, out->on_dev,
    lop->on_dev, rop->on_dev);
}

__global__ void
gkyl_dg_dot_product_op_range_cu_kernel(struct gkyl_basis basis,
  struct gkyl_array* out, const struct gkyl_array* lop,
  const struct gkyl_array* rop, struct gkyl_range range)
{
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  mul_op_t mul_op = choose_ser_mul_kern(ndim, poly_order);

  int num_basis = basis.num_basis;
  int vcomp = lop->ncomp/out->ncomp;

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

    const double *lop_d = (const double*) gkyl_array_cfetch(lop, start);
    const double *rop_d = (const double*) gkyl_array_cfetch(rop, start);
    double *out_d = (double*) gkyl_array_fetch(out, start);
    for (int k=0; k<num_basis; k++) out_d[k] = 0.;

    for (int d=0; d<vcomp; d++) {
      double comp_out[20];  // MF 2022/09/08: Hardcoded to number of basis in 3x p=2.
      mul_op(lop_d+d*num_basis, rop_d+d*num_basis, comp_out);
      for (int k=0; k<num_basis; k++) out_d[k] += comp_out[k];
    }
  }
}

// Host-side wrapper for range-based dg dot product operation.
void
gkyl_dg_dot_product_op_range_cu(const struct gkyl_basis *basis,
  struct gkyl_array* out, const struct gkyl_array* lop,
  const struct gkyl_array* rop, const struct gkyl_range *range)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  assert(basis->num_basis <= 20); // MF 2022/09/08: see hardcode in kernel above.
  gkyl_dg_dot_product_op_range_cu_kernel<<<nblocks, nthreads>>>(*basis, out->on_dev,
    lop->on_dev, rop->on_dev, *range);
}

__global__ void
gkyl_dg_div_set_op_cu_kernel(struct gkyl_nmat *As, struct gkyl_nmat *xs,
  struct gkyl_basis basis, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  div_set_op_t div_set_op = choose_ser_div_set_kern(ndim, poly_order);

  for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x) {
    
    const double *lop_d = (const double*) gkyl_array_cfetch(lop, linc);
    const double *rop_d = (const double*) gkyl_array_cfetch(rop, linc);

    struct gkyl_mat A = gkyl_nmat_get(As, linc);
    struct gkyl_mat x = gkyl_nmat_get(xs, linc);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0);
    div_set_op(&A, &x, lop_d+c_lop*num_basis, rop_d+c_rop*num_basis);
  }
}

__global__ void
gkyl_dg_div_copy_sol_op_cu_kernel(struct gkyl_nmat *xs,
  struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out)
{
  int num_basis = basis.num_basis;
  for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x) {
    double *out_d = (double*) gkyl_array_fetch(out, linc);
    struct gkyl_mat x = gkyl_nmat_get(xs, linc);
    binop_div_copy_sol(&x, out_d+c_oop*num_basis);
  }  
}

// Host-side wrapper for dg division operation
void
gkyl_dg_div_op_cu(gkyl_dg_bin_op_mem *mem, const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop)
{
  // allocate memory for use in kernels
  struct gkyl_nmat *A_d = mem->As;
  struct gkyl_nmat *x_d = mem->xs;

  // construct matrices using CUDA kernel
  gkyl_dg_div_set_op_cu_kernel<<<out->nblocks, out->nthreads>>>(A_d->on_dev, x_d->on_dev,
    *basis, out->on_dev, c_lop, lop->on_dev, c_rop, rop->on_dev);
  // invert all matrices in batch mode
  bool status = gkyl_nmat_linsolve_lu_pa(mem->lu_mem, A_d, x_d);
  assert(status);
  // copy solution into array (also lives on the device)
  gkyl_dg_div_copy_sol_op_cu_kernel<<<out->nblocks, out->nthreads>>>(x_d->on_dev, *basis, c_oop, out->on_dev);

}

__global__ void
gkyl_dg_div_set_op_range_cu_kernel(struct gkyl_nmat *As, struct gkyl_nmat *xs,
  struct gkyl_basis basis, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range range)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  div_set_op_t div_set_op = choose_ser_div_set_kern(ndim, poly_order);

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

    const double *lop_d = (const double*) gkyl_array_cfetch(lop, start);
    const double *rop_d = (const double*) gkyl_array_cfetch(rop, start);

    struct gkyl_mat A = gkyl_nmat_get(As, linc1);
    struct gkyl_mat x = gkyl_nmat_get(xs, linc1);
    gkyl_mat_clear(&A, 0.0); gkyl_mat_clear(&x, 0.0);  

    div_set_op(&A, &x, lop_d+c_lop*num_basis, rop_d+c_rop*num_basis);
  }
}

__global__ void
gkyl_dg_div_copy_sol_op_range_cu_kernel(struct gkyl_nmat *xs,
  struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out, struct gkyl_range range)
{
  int num_basis = basis.num_basis;

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

    binop_div_copy_sol(&x, out_d+c_oop*num_basis);
  }  
}

// Host-side wrapper for range-based dg division operation
void
gkyl_dg_div_op_range_cu(gkyl_dg_bin_op_mem *mem, const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, const struct gkyl_range *range)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  // allocate memory for use in kernels
  struct gkyl_nmat *A_d = mem->As;
  struct gkyl_nmat *x_d = mem->xs;

  // construct matrices using CUDA kernel  
  gkyl_dg_div_set_op_range_cu_kernel<<<nblocks, nthreads>>>(A_d->on_dev,
    x_d->on_dev, *basis, out->on_dev, c_lop, lop->on_dev, c_rop, rop->on_dev, *range);
  // invert all matrices in batch mode
  bool status = gkyl_nmat_linsolve_lu_pa(mem->lu_mem, A_d, x_d);
  assert(status);
  // copy solution into array (also lives on the device)
  gkyl_dg_div_copy_sol_op_range_cu_kernel<<<nblocks, nthreads>>>(x_d->on_dev,
    *basis, c_oop, out->on_dev, *range);
}

__global__ void
gkyl_dg_inv_op_cu_kernel(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_iop, const struct gkyl_array* iop)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  inv_op_t inv_op = choose_ser_inv_kern(ndim, poly_order);
  assert(inv_op);

  for (unsigned long linc = START_ID; linc < NSIZE(out); linc += blockDim.x*gridDim.x) {
    const double *iop_d = (const double*) gkyl_array_cfetch(iop, linc);
    double *out_d = (double*) gkyl_array_fetch(out, linc);

    inv_op(iop_d+c_iop*num_basis, out_d+c_oop*num_basis);
  }  
}

// Host-side wrapper for dg inversion operation.
void
gkyl_dg_inv_op_cu(const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_iop, const struct gkyl_array* iop)
{
  gkyl_dg_inv_op_cu_kernel<<<out->nblocks, out->nthreads>>>(*basis, c_oop, out->on_dev,
    c_iop, iop->on_dev);
}

__global__ void
gkyl_dg_inv_op_range_cu_kernel(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_iop, const struct gkyl_array* iop, struct gkyl_range range)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;
  int poly_order = basis.poly_order;
  inv_op_t inv_op = choose_ser_inv_kern(ndim, poly_order);
  assert(inv_op);

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

    const double *iop_d = (const double*) gkyl_array_cfetch(iop, start);
    double *out_d = (double*) gkyl_array_fetch(out, start);

    inv_op(iop_d+c_iop*num_basis, out_d+c_oop*num_basis);
  }
}

// Host-side wrapper for range-based dg invtiplication operation
void
gkyl_dg_inv_op_range_cu(const struct gkyl_basis *basis,
  int c_oop, struct gkyl_array* out,
  int c_iop, const struct gkyl_array* iop, const struct gkyl_range *range)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_dg_inv_op_range_cu_kernel<<<nblocks, nthreads>>>(*basis, c_oop, out->on_dev,
    c_iop, iop->on_dev, *range);
}

__global__ void
gkyl_dg_calc_op_range_cu_kernel(struct gkyl_basis basis, int c_oop, struct gkyl_array *out,
  int c_iop, const struct gkyl_array *iop,
  struct gkyl_range range, enum gkyl_dg_op op)
{
  int num_basis = basis.num_basis;
  int ndim = basis.ndim;

  dp_op_t op_func = dg_get_op_func(op);
  double fact = // factor for rescaling return value of op_func
    op == GKYL_DG_OP_MEAN ? sqrt(pow(2,ndim)) : pow(2,ndim);

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

    const double *iop_d = (const double*) gkyl_array_cfetch(iop, start);
    double *out_d = (double*) gkyl_array_fetch(out, start);

    out_d[c_oop] =
      op_func(num_basis, iop_d+c_iop*num_basis)/fact;
  }
}

void
gkyl_dg_calc_op_range_cu(const struct gkyl_basis *basis, int c_oop, struct gkyl_array *out,
  int c_iop, const struct gkyl_array *iop,
  struct gkyl_range range, enum gkyl_dg_op op)
{
  gkyl_dg_calc_op_range_cu_kernel<<<out->nblocks, out->nthreads>>>(*basis, c_oop, out->on_dev,
    c_iop, iop->on_dev, range, op);
}
