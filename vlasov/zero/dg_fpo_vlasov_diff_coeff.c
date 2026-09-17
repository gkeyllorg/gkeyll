#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_fpo_vlasov_coeff_recovery.h>
#include <gkyl_fpo_vlasov_coeff_recovery_priv.h>
#include <gkyl_dg_fpo_vlasov_diff_coeff.h>
#include <gkyl_dg_fpo_vlasov_diff_coeff_priv.h>

#include <gkyl_array_rio.h>

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

void gkyl_calc_fpo_diff_coeff_recovery(const struct gkyl_fpo_vlasov_coeff_recovery* coeff_recovery,
  const struct gkyl_rect_grid *grid, 
  struct gkyl_basis pbasis, const struct gkyl_range *phase_range, const struct gkyl_range *conf_range, 
  const struct gkyl_array *gamma, 
  const struct gkyl_array *fpo_g, const struct gkyl_array *fpo_g_surf, 
  const struct gkyl_array *fpo_dgdv_surf, const struct gkyl_array *fpo_d2gdv2_surf, 
  struct gkyl_array *fpo_diff_coeff, struct gkyl_array *fpo_diff_coeff_surf, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_calc_fpo_diff_coeff_recovery_cu(coeff_recovery, grid, pbasis, phase_range, conf_range, gamma, fpo_g, fpo_g_surf, fpo_dgdv_surf, fpo_d2gdv2_surf, fpo_diff_coeff, fpo_diff_coeff_surf);
#endif

  int cdim = coeff_recovery->cdim;
  int pdim = coeff_recovery->pdim; 
  int vdim = pdim - cdim;

  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], conf_idxc[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idxc);
    gkyl_copy_int_arr(cdim, iter.idx, conf_idxc);

    long linp = gkyl_range_idx(phase_range, idxc);
    long linc = gkyl_range_idx(conf_range, conf_idxc);

    const double *fpo_dgdv_surf_c = gkyl_array_cfetch(fpo_dgdv_surf, linp);
    const double *fpo_d2gdv2_surf_c = gkyl_array_cfetch(fpo_d2gdv2_surf, linp);
    double *fpo_diff_coeff_c = gkyl_array_fetch(fpo_diff_coeff, linp);

    const double *gamma_c = gkyl_array_cfetch(gamma, linc);

    // Iterate over velocity space directions
    for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;

      // Diagonal terms of the diffusion tensor.
      // Always a 1D, 3-cell stencil.
      int update_dirs[2] = {0};
      update_dirs[0] = dir1;

      // Index into kernel list.
      const long *offsets = &(coeff_recovery->offsets[3*d1]);
      int keri = idx_to_inloup_ker(1, idxc, update_dirs, phase_range->upper);

      const double* fpo_g_stencil[3];
      for (int i=0; i<3; ++i) {
        fpo_g_stencil[i] = gkyl_array_cfetch(fpo_g, linp+offsets[i]);
      }
     
      // Compute diagonal element of diffusion tensor
      coeff_recovery->diff_coeff_diag_recovery_stencil[d1][keri](grid->dx, gamma_c, 
        fpo_g_stencil, fpo_d2gdv2_surf_c, fpo_diff_coeff_c);

      for (int d2=0; d2<vdim; ++d2) {
        if (d1 == d2) continue;
        int dir2 = d2+cdim;

        // Off-diagonal terms of the diffusion tensor.
        // Offsets that would be outside the grid will point to center cell.
        // Always 2D and we need 9 cell stencil for 2D recovery.
        int update_dirs[2] = {0};
        update_dirs[0] = dir1 < dir2 ? dir1 : dir2;
        update_dirs[1] = dir1 < dir2 ? dir2 : dir1;

        // Index into kernel list
        const long *offsets = &(coeff_recovery->offsets[(d1+d2)*9]);
        int keri = idx_to_inloup_ker(2, idxc, update_dirs, phase_range->upper);

        const double *fpo_g_stencil[9], *fpo_g_surf_stencil[9];
        for (int i=0; i<9; ++i) {
          fpo_g_stencil[i] = gkyl_array_cfetch(fpo_g, linp+offsets[i]);
          fpo_g_surf_stencil[i] = gkyl_array_cfetch(fpo_g_surf, linp+offsets[i]);
        }

        coeff_recovery->diff_coeff_cross_recovery_stencil[d1][d2][keri](grid->dx, 
          gamma_c, fpo_g_stencil,
          fpo_g_surf_stencil, fpo_dgdv_surf_c, fpo_diff_coeff_c);
      }
    }
  }

  // Loop back over phase space to calculate surface expansions on LOWER cell boundary
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idxc);
    long linp_c = gkyl_range_idx(phase_range, idxc);
    double *fpo_diff_coeff_surf_c = gkyl_array_fetch(fpo_diff_coeff_surf, linp_c);
    
    // Iterate over primary direction for recovery, the kernel will handle
    // populating the three directions for the derivative across/along that boundary
    for (int d1=0; d1<vdim; ++d1) {
      int dir1 = d1 + cdim;

      // Do nothing if we're at a lower boundary; this surface expansion won't be used 
      if (idxc[dir1] == 1) {
        continue; 
      }
      
      gkyl_copy_int_arr(pdim, idxc, idxl); 
      idxl[dir1] = idxc[dir1]-1; 
      long linp_l = gkyl_range_idx(phase_range, idxl);
      
      const double *fpo_diff_coeff_l = gkyl_array_cfetch(fpo_diff_coeff, linp_l);
      const double *fpo_diff_coeff_c = gkyl_array_cfetch(fpo_diff_coeff, linp_c);

      coeff_recovery->diff_coeff_surf_recovery[d1](fpo_diff_coeff_l, fpo_diff_coeff_c,
        fpo_diff_coeff_surf_c);
    }
  }
}
