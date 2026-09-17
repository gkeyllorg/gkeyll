#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_fpo_vlasov_coeff_recovery.h>
#include <gkyl_fpo_vlasov_coeff_recovery_priv.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff.h>
#include <gkyl_dg_fpo_vlasov_drag_coeff_priv.h>
#include <gkyl_util.h>

static int idx_to_inloup_ker(int dim, const int *idx, const int *dirs, const int *num_cells) {
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

void gkyl_calc_fpo_drag_coeff_recovery(const struct gkyl_fpo_vlasov_coeff_recovery *coeff_recovery, 
  const struct gkyl_rect_grid *grid, struct gkyl_basis pbasis, 
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  const struct gkyl_array *gamma, 
  const struct gkyl_array *fpo_h, const struct gkyl_array* fpo_dhdv_surf, 
  struct gkyl_array *fpo_drag_coeff, struct gkyl_array *fpo_drag_coeff_surf, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_calc_fpo_drag_coeff_recovery_cu(coeff_recovery, grid, pbasis, phase_range, conf_range, gamma, fpo_h, fpo_dhdv_surf, fpo_drag_coeff, fpo_drag_coeff_surf);
#endif

  int pdim = coeff_recovery->pdim;
  int vdim = coeff_recovery->vdim;
  int cdim = coeff_recovery->cdim;

  // Indices in each direction
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM], conf_idxc[GKYL_MAX_DIM];
 
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idxc);
    gkyl_copy_int_arr(cdim, iter.idx, conf_idxc);

    long linp = gkyl_range_idx(phase_range, idxc);
    long linc = gkyl_range_idx(conf_range, conf_idxc);
 
    const double *fpo_dhdv_surf_c = gkyl_array_cfetch(fpo_dhdv_surf, linp);
    double *fpo_drag_coeff_c = gkyl_array_fetch(fpo_drag_coeff, linp);
    double *fpo_drag_coeff_surf_c = gkyl_array_fetch(fpo_drag_coeff_surf, linp);

    const double *gamma_c = gkyl_array_cfetch(gamma, linc);

    // Iterate through velocity space directions
    for (int d=0; d<vdim; ++d) {
      int dir = d + cdim;
      int update_dir[] = {dir};

      // Index into kernel list.
      const long *offsets = &(coeff_recovery->offsets[d*3]);
      int keri = idx_to_inloup_ker(1, idxc, update_dir, phase_range->upper);

      const double* fpo_h_stencil[3];

      for (int i=0; i<3; ++i) {
        fpo_h_stencil[i] = gkyl_array_cfetch(fpo_h, linp+offsets[i]);
      }

      coeff_recovery->drag_coeff_recovery_stencil[d][keri](
        grid->dx, gamma_c, fpo_h_stencil, 
        fpo_dhdv_surf_c, fpo_drag_coeff_c, fpo_drag_coeff_surf_c);
    }
  }
}

GKYL_CU_DH
void gkyl_calc_fpo_sgn_drag_coeff(const struct gkyl_fpo_vlasov_coeff_recovery *coeff_recovery, 
  struct gkyl_basis pbasis, const struct gkyl_range *phase_range, 
  struct gkyl_array* fpo_drag_coeff_surf, struct gkyl_array* sgn_drag_coeff_surf, 
  struct gkyl_array* const_sgn_drag_coeff_surf, bool use_gpu) 
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_calc_fpo_sgn_drag_coeff_cu(coeff_recovery, pbasis, phase_range, fpo_drag_coeff_surf, sgn_drag_coeff_surf, const_sgn_drag_coeff_surf);
  }
#endif

  int pdim = coeff_recovery->pdim;
  int vdim = coeff_recovery->vdim;
  int cdim = coeff_recovery->cdim;

  int poly_order = pbasis.poly_order; 

  fpo_sgn_drag_coeff_t sgn_drag_coeff_stencil[3][3];
 
  // Fetch kernels in each direction
  for (int d=0; d<vdim; ++d) {
    for (int idx=0; idx<3; ++idx) {
      sgn_drag_coeff_stencil[d][idx] = 
       choose_ser_fpo_sgn_drag_coeff_recovery_kern(d, cdim, poly_order, idx);
    }
  }
  
  int idxp[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idxp);

    long linp = gkyl_range_idx(phase_range, idxp);

    const double *fpo_drag_coeff_surf_c = gkyl_array_cfetch(fpo_drag_coeff_surf, linp);
    double *sgn_drag_coeff_surf_c = gkyl_array_fetch(sgn_drag_coeff_surf, linp);
    int *const_sgn_drag_coeff_surf_c = gkyl_array_fetch(const_sgn_drag_coeff_surf, linp);

    // Iterate through velocity space directions
    for (int d=0; d<vdim; ++d) {
      int dir = d + cdim;
      int update_dir[] = {dir};

      int keri = idx_to_inloup_ker(1, idxp, update_dir, phase_range->upper);

      sgn_drag_coeff_stencil[d][keri](
        fpo_drag_coeff_surf_c, sgn_drag_coeff_surf_c, &const_sgn_drag_coeff_surf_c[d]);
    }
  }
}
