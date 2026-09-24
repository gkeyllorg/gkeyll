#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_fpo_vlasov_diff.h>
#include <gkyl_dg_fpo_vlasov_diff_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polynomial order
#define CK(lst, cdim, poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_fpo_vlasov_diff_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_fpo_vlasov_diff* fpo_vlasov_diff = container_of(base->on_dev, struct dg_fpo_vlasov_diff, eqn);
    gkyl_cu_free(fpo_vlasov_diff);
  }
  
  struct dg_fpo_vlasov_diff* fpo_vlasov_diff = container_of(base, struct dg_fpo_vlasov_diff, eqn);
  gkyl_free(fpo_vlasov_diff);
}

void
gkyl_fpo_vlasov_diff_set_auxfields(const struct gkyl_dg_eqn* eqn, struct gkyl_dg_fpo_vlasov_diff_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.g)) {
    gkyl_fpo_vlasov_diff_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif
  
  struct dg_fpo_vlasov_diff* fpo_vlasov_diff = container_of(eqn, struct dg_fpo_vlasov_diff, eqn);
  fpo_vlasov_diff->auxfields.g = auxin.g;
}

struct gkyl_dg_eqn*
gkyl_dg_fpo_vlasov_diff_new(const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu)
    return gkyl_dg_fpo_vlasov_diff_cu_dev_new(pbasis, phase_range);
#endif
  
  struct dg_fpo_vlasov_diff* fpo_vlasov_diff = gkyl_malloc(sizeof(struct dg_fpo_vlasov_diff));

  // Vlasov Fokker-Planck operator only defined in 3 velocity dimensions
  int pdim = pbasis->ndim, vdim = 3, cdim = pdim - vdim;
  int poly_order = pbasis->poly_order;

  fpo_vlasov_diff->cdim = cdim;
  fpo_vlasov_diff->pdim = pdim;

  fpo_vlasov_diff->eqn.num_equations = 1;
  fpo_vlasov_diff->eqn.gen_surf_term = surf;
  fpo_vlasov_diff->eqn.gen_boundary_surf_term = boundary_surf;

  const gkyl_dg_fpo_vlasov_diff_vol_kern_list* vol_kernels;
  const gkyl_dg_fpo_vlasov_diff_surf_kern_list* surf_xx_kernels;
  const gkyl_dg_fpo_vlasov_diff_surf_kern_list* surf_xy_kernels;
  const gkyl_dg_fpo_vlasov_diff_surf_kern_list* surf_xz_kernels;
  const gkyl_dg_fpo_vlasov_diff_surf_kern_list* surf_yx_kernels;
  const gkyl_dg_fpo_vlasov_diff_surf_kern_list* surf_yy_kernels;
  const gkyl_dg_fpo_vlasov_diff_surf_kern_list* surf_yz_kernels;
  const gkyl_dg_fpo_vlasov_diff_surf_kern_list* surf_zx_kernels;
  const gkyl_dg_fpo_vlasov_diff_surf_kern_list* surf_zy_kernels;
  const gkyl_dg_fpo_vlasov_diff_surf_kern_list* surf_zz_kernels; 

  const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list* boundary_surf_xx_kernels;
  const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list* boundary_surf_xy_kernels;
  const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list* boundary_surf_xz_kernels;
  const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list* boundary_surf_yx_kernels;
  const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list* boundary_surf_yy_kernels;
  const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list* boundary_surf_yz_kernels;
  const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list* boundary_surf_zx_kernels;
  const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list* boundary_surf_zy_kernels;
  const gkyl_dg_fpo_vlasov_diff_boundary_surf_kern_list* boundary_surf_zz_kernels; 

  switch (pbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_xx_kernels = ser_surf_xx_kernels;
      surf_xy_kernels = ser_surf_xy_kernels;
      surf_xz_kernels = ser_surf_xz_kernels;
      surf_yx_kernels = ser_surf_yx_kernels;
      surf_yy_kernels = ser_surf_yy_kernels;
      surf_yz_kernels = ser_surf_yz_kernels;
      surf_zx_kernels = ser_surf_zx_kernels;
      surf_zy_kernels = ser_surf_zy_kernels;
      surf_zz_kernels = ser_surf_zz_kernels;

      boundary_surf_xx_kernels = ser_boundary_surf_xx_kernels;
      boundary_surf_xy_kernels = ser_boundary_surf_xy_kernels;
      boundary_surf_xz_kernels = ser_boundary_surf_xz_kernels;
      boundary_surf_yx_kernels = ser_boundary_surf_yx_kernels;
      boundary_surf_yy_kernels = ser_boundary_surf_yy_kernels;
      boundary_surf_yz_kernels = ser_boundary_surf_yz_kernels;
      boundary_surf_zx_kernels = ser_boundary_surf_zx_kernels;
      boundary_surf_zy_kernels = ser_boundary_surf_zy_kernels;
      boundary_surf_zz_kernels = ser_boundary_surf_zz_kernels;
      break;

    default:
      assert(false);
      break;    
  } 

  fpo_vlasov_diff->eqn.vol_term = CK(vol_kernels, cdim, poly_order);

  fpo_vlasov_diff->surf[0][0] = CK(surf_xx_kernels, cdim, poly_order);
  fpo_vlasov_diff->surf[0][1] = CK(surf_xy_kernels, cdim, poly_order);
  fpo_vlasov_diff->surf[0][2] = CK(surf_xz_kernels, cdim, poly_order);
  fpo_vlasov_diff->surf[1][0] = CK(surf_yx_kernels, cdim, poly_order);
  fpo_vlasov_diff->surf[1][1] = CK(surf_yy_kernels, cdim, poly_order);
  fpo_vlasov_diff->surf[1][2] = CK(surf_yz_kernels, cdim, poly_order);
  fpo_vlasov_diff->surf[2][0] = CK(surf_zx_kernels, cdim, poly_order);
  fpo_vlasov_diff->surf[2][1] = CK(surf_zy_kernels, cdim, poly_order);
  fpo_vlasov_diff->surf[2][2] = CK(surf_zz_kernels, cdim, poly_order);

  fpo_vlasov_diff->boundary_surf[0][0] = CK(boundary_surf_xx_kernels, cdim, poly_order);
  fpo_vlasov_diff->boundary_surf[0][1] = CK(boundary_surf_xy_kernels, cdim, poly_order);
  fpo_vlasov_diff->boundary_surf[0][2] = CK(boundary_surf_xz_kernels, cdim, poly_order);
  fpo_vlasov_diff->boundary_surf[1][0] = CK(boundary_surf_yx_kernels, cdim, poly_order);
  fpo_vlasov_diff->boundary_surf[1][1] = CK(boundary_surf_yy_kernels, cdim, poly_order);
  fpo_vlasov_diff->boundary_surf[1][2] = CK(boundary_surf_yz_kernels, cdim, poly_order);
  fpo_vlasov_diff->boundary_surf[2][0] = CK(boundary_surf_zx_kernels, cdim, poly_order);
  fpo_vlasov_diff->boundary_surf[2][1] = CK(boundary_surf_zy_kernels, cdim, poly_order);
  fpo_vlasov_diff->boundary_surf[2][2] = CK(boundary_surf_zz_kernels, cdim, poly_order);

  // ensure non-NULL pointers
  for (int i=0; i<vdim; ++i) 
    for (int j=0; j<vdim; ++j) 
      assert(fpo_vlasov_diff->surf[i][j]);

  for (int i=0; i<vdim; ++i) 
    for (int j=0; j<vdim; ++j) 
      assert(fpo_vlasov_diff->boundary_surf[i][j]);

  fpo_vlasov_diff->auxfields.g = 0;
  fpo_vlasov_diff->phase_range = *phase_range;

  fpo_vlasov_diff->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(fpo_vlasov_diff->eqn.flags);
  fpo_vlasov_diff->eqn.ref_count = gkyl_ref_count_init(gkyl_fpo_vlasov_diff_free);
  fpo_vlasov_diff->eqn.on_dev = &fpo_vlasov_diff->eqn;
  
  return &fpo_vlasov_diff->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_fpo_vlasov_diff_cu_dev_new(const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range)
{
  assert(false);
  return 0;
}

#endif
