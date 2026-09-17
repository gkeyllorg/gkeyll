extern "C" {
  #include <gkyl_alloc.h>
  #include <gkyl_alloc_flags_priv.h>
  #include <gkyl_array_ops.h>
  #include <gkyl_array_ops_priv.h>
  #include <gkyl_fpo_vlasov_coeff_recovery.h>
  #include <gkyl_fpo_vlasov_coeff_recovery_priv.h>
  #include <gkyl_util.h>
}


__global__ static void fpo_vlasov_coeff_recovery_set_cu_dev_ptrs(gkyl_fpo_vlasov_coeff_recovery *up, 
  int cdim, int poly_order)
{
  int vdim = 3;

  for (int d1=0; d1<vdim; ++d1) {
    up->diff_coeff_surf_recovery[d1] = 
      choose_ser_fpo_diff_coeff_surf_recovery_kern(d1, cdim, poly_order);

    // 3-cell stencil pointers
    for (int i=0; i<3; ++i) {
      up->drag_coeff_recovery_stencil[d1][i] = 
       choose_ser_fpo_drag_coeff_recovery_kern(d1, cdim, poly_order, i);

      up->sgn_drag_coeff_stencil[d1][i] = 
       choose_ser_fpo_sgn_drag_coeff_recovery_kern(d1, cdim, poly_order, i);

      up->diff_coeff_diag_recovery_stencil[d1][i] = 
        choose_ser_fpo_diff_coeff_diag_recovery_kern(d1, cdim, poly_order, i);
    }

    for (int d2=0; d2<vdim; ++d2) {
      if (d1 != d2) {
        for (int idx=0; idx<9; ++idx) {
          up->diff_coeff_cross_recovery_stencil[d1][d2][idx] = 
            choose_ser_fpo_diff_coeff_cross_recovery_kern(d1, d2, cdim, poly_order, idx);
        }
      }
    }
  }

}

gkyl_fpo_vlasov_coeff_recovery*
gkyl_fpo_vlasov_coeff_recovery_cu_dev_new(const struct gkyl_rect_grid *grid, 
    const struct gkyl_basis *phase_basis, const struct gkyl_range *phase_range, long offsets[36])
{
  gkyl_fpo_vlasov_coeff_recovery *up = (gkyl_fpo_vlasov_coeff_recovery *)gkyl_malloc(sizeof(gkyl_fpo_vlasov_coeff_recovery));

  // FPO always 3V
  int vdim = 3;
  int pdim = phase_basis->ndim;
  int cdim = pdim - vdim;

  int poly_order = phase_basis->poly_order;

  up->cdim = cdim;
  up->vdim = vdim;
  up->pdim = pdim;
  up->poly_order = poly_order;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  // Create array of 36 relative offsets.
  // 3-cell stencil in each velocity direction, (vx, vy, vz)
  // 9-cell stencil in some pairs of velocity directions (vxvy, vxvz, vyvz)
  int idxc[GKYL_MAX_DIM] = {0};

  for (int i=cdim; i<pdim; ++i) idxc[i] = 2;

  long offsets_arr[36] = {0};
  for (int d1=0; d1<vdim; ++d1) {
    int dir1 = d1 + cdim;
    int num_update_dir = 1;
    int update_dir[] = {dir1};
    create_offsets(num_update_dir, update_dir, phase_range, idxc, &offsets_arr[d1*3]);

    for (int d2=d1+1; d2<vdim; ++d2) {
      int dir2 = d2 + cdim;
      int num_update_dirs = 2;
      int update_dirs[] = {dir1, dir2};
      create_offsets(num_update_dirs, update_dirs, phase_range, idxc,
        &offsets_arr[(d1+d2)*9]);
    }
  }
  gkyl_copy_long_arr(36, offsets_arr, up->offsets);  
  gkyl_copy_long_arr(36, offsets_arr, offsets);

  gkyl_fpo_vlasov_coeff_recovery *up_cu = (gkyl_fpo_vlasov_coeff_recovery *)
    gkyl_cu_malloc(sizeof(gkyl_fpo_vlasov_coeff_recovery));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_fpo_vlasov_coeff_recovery), GKYL_CU_MEMCPY_H2D);
  fpo_vlasov_coeff_recovery_set_cu_dev_ptrs<<<1,1>>>(up_cu, cdim, poly_order);

  // set pointer to device struct
  up->on_dev = up_cu;

  return up;
}
