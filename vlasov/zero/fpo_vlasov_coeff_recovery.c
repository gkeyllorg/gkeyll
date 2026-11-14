
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_fpo_vlasov_coeff_recovery.h>
#include <gkyl_fpo_vlasov_coeff_recovery_priv.h>
#include <gkyl_util.h>

gkyl_fpo_vlasov_coeff_recovery *
gkyl_fpo_vlasov_coeff_recovery_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *phase_basis,
  const struct gkyl_range *phase_range, long offsets[GKYL_HYPER_DG_GEN_NUM_OFFSETS], bool use_gpu
)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_fpo_vlasov_coeff_recovery_cu_dev_new(grid, phase_basis, phase_range, offsets);
  }
#endif

  gkyl_fpo_vlasov_coeff_recovery *up = gkyl_malloc(sizeof(struct gkyl_fpo_vlasov_coeff_recovery));

  // FPO always 3V
  int vdim = 3;
  int pdim = phase_basis->ndim;
  int cdim = pdim - vdim;

  int poly_order = phase_basis->poly_order;

  up->cdim = cdim;
  up->vdim = vdim;
  up->pdim = pdim;
  up->poly_order = poly_order;

  // Relative offsets of the 3-cell stencil in each velocity direction (vx, vy, vz) and
  // of the 9-cell stencil in each pair of velocity directions (vxvy, vxvz, vyvz). These
  // are the stencils of the generic stencil hyper DG update, which the FPO diffusion
  // term shares with the drag and diffusion coefficient recovery.
  int vel_dirs[GKYL_MAX_DIM] = {0};
  for (int d = 0; d < vdim; ++d) {
    vel_dirs[d] = cdim + d;
  }
  long offsets_arr[GKYL_HYPER_DG_GEN_NUM_OFFSETS] = {0};
  gkyl_hyper_dg_gen_stencil_offsets(vdim, vel_dirs, phase_range, offsets_arr);
  gkyl_copy_long_arr(GKYL_HYPER_DG_GEN_NUM_OFFSETS, offsets_arr, up->offsets);
  gkyl_copy_long_arr(GKYL_HYPER_DG_GEN_NUM_OFFSETS, offsets_arr, offsets);

  // Set pointers to kernels
  for (int d1 = 0; d1 < vdim; ++d1) {
    up->diff_coeff_surf_recovery[d1] =
      choose_ser_fpo_diff_coeff_surf_recovery_kern(d1, cdim, poly_order);

    // 3-cell stencil pointers
    for (int i = 0; i < 3; ++i) {
      up->drag_coeff_recovery_stencil[d1][i] =
        choose_ser_fpo_drag_coeff_recovery_kern(d1, cdim, poly_order, i);

      up->sgn_drag_coeff_stencil[d1][i] =
        choose_ser_fpo_sgn_drag_coeff_recovery_kern(d1, cdim, poly_order, i);

      up->diff_coeff_diag_recovery_stencil[d1][i] =
        choose_ser_fpo_diff_coeff_diag_recovery_kern(d1, cdim, poly_order, i);
    }

    for (int d2 = 0; d2 < vdim; ++d2) {
      if (d1 != d2) {
        // 9-cell stencil pointers
        for (int idx = 0; idx < 9; ++idx) {
          up->diff_coeff_cross_recovery_stencil[d1][d2][idx] =
            choose_ser_fpo_diff_coeff_cross_recovery_kern(d1, d2, cdim, poly_order, idx);
        }
      }
    }
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);

  up->on_dev = up;

  return up;
}

void
gkyl_fpo_vlasov_coeff_recovery_release(struct gkyl_fpo_vlasov_coeff_recovery *up)
{
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}
