extern "C" {
#include <assert.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_fpo_vlasov_coeff_correct.h>
#include <gkyl_fpo_vlasov_coeff_correct_priv.h>
}

__global__ static void
gkyl_fpo_coeff_correct_set_mat_cu_ker(gkyl_fpo_coeff_correct *up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array *fpo_moms, const struct gkyl_array *boundary_corrections,
  const struct gkyl_array *moms)
{
  int cidx[GKYL_MAX_DIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume;
      tid += gridDim.x*blockDim.x)
  {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);
    long linc = gkyl_range_idx(&conf_range, cidx);

    const double *fpo_moms_d = (const double *)gkyl_array_cfetch(fpo_moms, linc);
    const double *boundary_corrections_d = (const double *)gkyl_array_cfetch(boundary_corrections, linc);
    const double *moms_d = (const double *)gkyl_array_cfetch(moms, linc); 

    struct gkyl_mat lhs = gkyl_nmat_get(As, tid);
    struct gkyl_mat rhs = gkyl_nmat_get(xs, tid);
    gkyl_mat_clear(&lhs, 0.0);
    gkyl_mat_clear(&rhs, 0.0);

    // Set matrix elements
    up->mat_set_kernel(&lhs, &rhs, fpo_moms_d, boundary_corrections_d, moms_d);
  }
}

__global__ static void
gkyl_fpo_coeff_correct_copy_sol_cu_ker(gkyl_fpo_coeff_correct *up,
  struct gkyl_range conf_range, struct gkyl_nmat *xs, struct gkyl_array *drag_diff_coeff_corrs)
{
  int nc = up->num_conf_basis;
  int vdim = 3;
  int N = nc*(vdim + 1);

  int cidx[GKYL_MAX_DIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_range.volume;
      tid += gridDim.x*blockDim.x)
  {
    gkyl_sub_range_inv_idx(&conf_range, tid, cidx);
    long linc = gkyl_range_idx(&conf_range, cidx);

    double *drag_diff_coeff_corrs_d  = (double *)gkyl_array_fetch(drag_diff_coeff_corrs, linc);
    const struct gkyl_mat out = gkyl_nmat_get(xs, tid);
    for (size_t i=0; i<N; ++i) {
      drag_diff_coeff_corrs_d[i] = gkyl_mat_get(&out, i, 0);
    }
  }
}

__global__ static void
gkyl_fpo_coeff_correct_accum_cu_ker(gkyl_fpo_coeff_correct *up,
  struct gkyl_nmat *xs, struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *drag_diff_coeff_corrs, struct gkyl_array *drag_coeff,
  struct gkyl_array *drag_coeff_surf, struct gkyl_array *diff_coeff, struct gkyl_array *diff_coeff_surf)
{
  int pidx[GKYL_MAX_DIM];
  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume;
      tid += gridDim.x*blockDim.x)
  {
    gkyl_sub_range_inv_idx(&phase_range, tid, pidx);
    long linc = gkyl_range_idx(&conf_range, pidx);
    long linp = gkyl_range_idx(&phase_range, pidx);

    double *drag_coeff_d = (double *)gkyl_array_fetch(drag_coeff, linp);
    double *drag_coeff_surf_d = (double *)gkyl_array_fetch(drag_coeff_surf, linp);
    double *diff_coeff_d = (double *)gkyl_array_fetch(diff_coeff, linp);
    double *diff_coeff_surf_d = (double *)gkyl_array_fetch(diff_coeff_surf, linp);

    const double *drag_diff_coeff_corrs_d  = (const double *)gkyl_array_cfetch(drag_diff_coeff_corrs, linc);

    // Call to kernels to accumulate corrections
    up->accum_kernel(drag_diff_coeff_corrs_d, drag_coeff_d, drag_coeff_surf_d, 
      diff_coeff_d, diff_coeff_surf_d); 
  }
}


__global__ static void
gkyl_fpo_coeff_correct_set_cu_dev_ptrs(gkyl_fpo_coeff_correct *up, gkyl_basis_type b_type, int cdim, int poly_order) 
{
  // Kernels for setting linear system matrices
  const gkyl_fpo_coeff_correct_mat_set_kern_list *mat_set_kern_list;
  const gkyl_fpo_coeff_correct_accum_kern_list *accum_kern_list;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mat_set_kern_list = ser_fpo_coeff_correct_mat_set_kernels;
      accum_kern_list = ser_fpo_coeff_correct_accum_kernels;
      break;

    default:
      assert(false);
      break;
  }

  up->mat_set_kernel = mat_set_kern_list[cdim].kernels[poly_order];
  up->accum_kernel = accum_kern_list[cdim].kernels[poly_order];
}

gkyl_fpo_coeff_correct*
gkyl_fpo_coeff_correct_cu_dev_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_range)
{
  struct gkyl_fpo_coeff_correct *up = (struct gkyl_fpo_coeff_correct *)gkyl_malloc(sizeof(struct gkyl_fpo_coeff_correct)); 

  int cdim = conf_basis->ndim;
  int poly_order = conf_basis->poly_order;
  
  up->grid = grid;
  up->conf_basis = conf_basis;
  up->num_conf_basis = conf_basis->num_basis;

  up->is_first = true;

  // Matrices and memory for linear solve
  up->As = 0;
  up->xs = 0;
  up->mem = 0;

  struct gkyl_fpo_coeff_correct *up_cu = (struct gkyl_fpo_coeff_correct *)
    gkyl_cu_malloc(sizeof(struct gkyl_fpo_coeff_correct)); 
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_fpo_coeff_correct), GKYL_CU_MEMCPY_H2D);

  gkyl_fpo_coeff_correct_set_cu_dev_ptrs<<<1,1>>>(up_cu, conf_basis->b_type, cdim, poly_order);

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  // set pointer to device struct
  up->on_dev = up_cu;

  return up;
}

void gkyl_fpo_coeff_correct_advance_cu(gkyl_fpo_coeff_correct *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_array *fpo_moms, const struct gkyl_array *boundary_corrections,
  const struct gkyl_array *moms, struct gkyl_array *drag_diff_coeff_corrs,
  struct gkyl_array *drag_coeff, struct gkyl_array *drag_coeff_surf,
  struct gkyl_array *diff_coeff, struct gkyl_array *diff_coeff_surf)
{
  // allocate memory for use in kernels
  int nc = up->num_conf_basis;
  int vdim = 3;
  int N = nc*(vdim + 1);

  // Initialize matrices if this is the first call to advance
  if (up->is_first) {
    up->As = gkyl_nmat_cu_dev_new(conf_range->volume, N, N);
    up->xs = gkyl_nmat_cu_dev_new(conf_range->volume, N, 1);
    up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);
    up->is_first = false;
  }

  // Set matrices
  gkyl_fpo_coeff_correct_set_mat_cu_ker<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, up->As->on_dev, up->xs->on_dev, *conf_range,
    fpo_moms->on_dev, boundary_corrections->on_dev, moms->on_dev);

  // Solve linear system
  bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
  assert(status);

  // Copy corrections from matrix solve memory
  gkyl_fpo_coeff_correct_copy_sol_cu_ker<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, up->xs->on_dev, drag_diff_coeff_corrs->on_dev);

  // Accumulate corrections onto drag/diffusion coefficients
  gkyl_fpo_coeff_correct_accum_cu_ker<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, up->xs->on_dev, *conf_range, *phase_range, drag_diff_coeff_corrs->on_dev,
    drag_coeff->on_dev, drag_coeff_surf->on_dev, diff_coeff->on_dev, diff_coeff_surf->on_dev);
}
