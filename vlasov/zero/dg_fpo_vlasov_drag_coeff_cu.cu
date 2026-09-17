extern "C" {
  #include <gkyl_alloc.h>
  #include <gkyl_array.h>
  #include <gkyl_array_ops.h>
  #include <gkyl_array_ops_priv.h>
  #include <gkyl_fpo_vlasov_coeff_recovery.h>
  #include <gkyl_fpo_vlasov_coeff_recovery_priv.h>
  #include <gkyl_dg_fpo_vlasov_drag_coeff.h>
  #include <gkyl_dg_fpo_vlasov_drag_coeff_priv.h>
  #include <gkyl_util.h>
}

__device__ static
int idx_to_inloup_ker(int dim, const int *idx, const int *dirs, const int *num_cells) {
  int iout = 0;

  for (int d=0; d<dim; ++d) {
    if (idx[dirs[d]] == 1) {
      iout = 2*iout+(int)(pow(3,d)+0.5);
    } else if (idx[dirs[d]] == num_cells[dirs[d]]) {
      iout = 2*iout+(int)(pow(3,d)+0.5)+1;
    }
  }
  return iout;
}

__global__ static void
gkyl_calc_fpo_drag_coeff_recovery_cu_kernel(const struct gkyl_fpo_vlasov_coeff_recovery *coeff_recovery,
  const struct gkyl_rect_grid grid,
  struct gkyl_basis pbasis, const struct gkyl_range phase_range, 
  const struct gkyl_range conf_range,
  const struct gkyl_array *gamma, const struct gkyl_array* fpo_h,
  const struct gkyl_array* fpo_dhdv_surf, struct gkyl_array* fpo_drag_coeff,
  struct gkyl_array* fpo_drag_coeff_surf)
{
  int cdim = coeff_recovery->cdim;
  int vdim = coeff_recovery->vdim;

  int idxp[GKYL_MAX_DIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume;
      tid += gridDim.x*blockDim.x)
  {
    gkyl_sub_range_inv_idx(&phase_range, tid, idxp);
    long linp = gkyl_range_idx(&phase_range, idxp);
    long linc = gkyl_range_idx(&conf_range, idxp);

    const double *fpo_dhdv_surf_c = (const double*)gkyl_array_cfetch(fpo_dhdv_surf, linp);
    double *fpo_drag_coeff_c = (double *)gkyl_array_fetch(fpo_drag_coeff, linp);
    double *fpo_drag_coeff_surf_c = (double *)gkyl_array_fetch(fpo_drag_coeff_surf, linp);

    const double *gamma_c = (const double *)gkyl_array_cfetch(gamma, linc);

    for (int d=0; d<vdim; ++d) {
      int dir = d + cdim;

      // Always a 1D, 3-cell stencil.
      int update_dir[] = {dir};

      // Index into kernel list.
      const long *offsets = &(coeff_recovery->offsets[d*3]);
      int keri = idx_to_inloup_ker(1, idxp, update_dir, phase_range.upper);

      const double* fpo_h_stencil[3];

      for (int i=0; i<3; ++i) {
        fpo_h_stencil[i] = (const double *)gkyl_array_cfetch(fpo_h, linp+offsets[i]);
      }

      coeff_recovery->drag_coeff_recovery_stencil[d][keri](
        grid.dx, gamma_c, fpo_h_stencil, 
        fpo_dhdv_surf_c, fpo_drag_coeff_c, fpo_drag_coeff_surf_c);
    }
  }
}

__global__ static void 
gkyl_calc_fpo_sgn_drag_coeff_cu_kernel(const struct gkyl_fpo_vlasov_coeff_recovery *coeff_recovery,
  struct gkyl_basis pbasis, const struct gkyl_range phase_range,
  struct gkyl_array* fpo_drag_coeff_surf, struct gkyl_array* sgn_drag_coeff_surf, struct gkyl_array* const_sgn_drag_coeff_surf)
{
  int cdim = coeff_recovery->cdim;
  int vdim = coeff_recovery->pdim - cdim;

  int idxp[GKYL_MAX_DIM];

  for (unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume;
      tid += gridDim.x*blockDim.x)
  {
    gkyl_sub_range_inv_idx(&phase_range, tid, idxp);
    long linp = gkyl_range_idx(&phase_range, idxp);

    const double *fpo_drag_coeff_surf_d = (const double *)gkyl_array_cfetch(fpo_drag_coeff_surf, linp);
    double *sgn_drag_coeff_surf_d = (double *)gkyl_array_fetch(sgn_drag_coeff_surf, linp);
    int *const_sgn_drag_coeff_surf_d = (int *)gkyl_array_fetch(const_sgn_drag_coeff_surf, linp);
 
    for (int d=0; d<vdim; ++d) {
      int dir = d + cdim;
      int update_dir[] = {dir};

      int keri = idx_to_inloup_ker(1, idxp, update_dir, phase_range.upper);

      coeff_recovery->sgn_drag_coeff_stencil[d][keri](fpo_drag_coeff_surf_d, sgn_drag_coeff_surf_d, &const_sgn_drag_coeff_surf_d[d]);
    } 
  }
}

// Host-side wrapper
void gkyl_calc_fpo_drag_coeff_recovery_cu(const struct gkyl_fpo_vlasov_coeff_recovery *coeff_recovery,
  const struct gkyl_rect_grid *grid, struct gkyl_basis pbasis, 
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array* gamma, const struct gkyl_array* fpo_h, 
  const struct gkyl_array* fpo_dhdv_surf, struct gkyl_array* fpo_drag_coeff,
  struct gkyl_array* fpo_drag_coeff_surf)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;

  gkyl_calc_fpo_drag_coeff_recovery_cu_kernel<<<nblocks, nthreads>>>(coeff_recovery->on_dev, 
    *grid, pbasis, *phase_range, *conf_range, 
    gamma->on_dev, fpo_h->on_dev, fpo_dhdv_surf->on_dev, 
    fpo_drag_coeff->on_dev, fpo_drag_coeff_surf->on_dev);
}

// Host-side wrapper
void gkyl_calc_fpo_sgn_drag_coeff_cu(const struct gkyl_fpo_vlasov_coeff_recovery *coeff_recovery,
  struct gkyl_basis pbasis, const struct gkyl_range *phase_range, 
  struct gkyl_array* fpo_drag_coeff_surf, struct gkyl_array* sgn_drag_coeff_surf, struct gkyl_array* const_sgn_drag_coeff_surf) 
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;

  gkyl_calc_fpo_sgn_drag_coeff_cu_kernel<<<nblocks, nthreads>>>(coeff_recovery->on_dev, 
    pbasis, *phase_range, fpo_drag_coeff_surf->on_dev, sgn_drag_coeff_surf->on_dev, const_sgn_drag_coeff_surf->on_dev);
}
