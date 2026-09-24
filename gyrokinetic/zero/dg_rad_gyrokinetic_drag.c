#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_rad_gyrokinetic_drag_priv.h>
#include <gkyl_util.h>

void
gkyl_rad_gyrokinetic_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_rad_gyrokinetic_drag *grad = container_of(base, struct dg_rad_gyrokinetic_drag, eqn);

  gkyl_velocity_map_release(grad->vel_map);

  if (GKYL_IS_CU_ALLOC(grad->eqn.flags))
    gkyl_cu_free(grad->eqn.on_dev);
  
  gkyl_free(grad);
}

void
gkyl_rad_gyrokinetic_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_rad_gyrokinetic_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    gkyl_rad_gyrokinetic_drag_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_rad_gyrokinetic_drag *grad = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  grad->auxfields.nvnu_surf = auxin.nvnu_surf;
  grad->auxfields.nvnu = auxin.nvnu;
  grad->auxfields.nvsqnu_surf = auxin.nvsqnu_surf;
  grad->auxfields.nvsqnu = auxin.nvsqnu;
}

struct gkyl_dg_eqn*
gkyl_dg_rad_gyrokinetic_drag_new(const struct gkyl_basis *conf_basis, 
  const struct gkyl_basis *phase_basis, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_velocity_map *vel_map, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)    
    return gkyl_dg_rad_gyrokinetic_drag_cu_dev_new(conf_basis, phase_basis, phase_range, conf_range, vel_map);
#endif
    
  struct dg_rad_gyrokinetic_drag *grad = gkyl_malloc(sizeof(*grad));

  int cdim = conf_basis->ndim, pdim = phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = conf_basis->poly_order;
  
  grad->cdim = cdim;
  grad->pdim = pdim;

  grad->cellav_norm_conf = 1.0/pow(sqrt(2.0),cdim);

  grad->eqn.num_equations = 1;
  grad->eqn.surf_term = surf;
  grad->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_rad_gyrokinetic_vol_kern_list *vol_kernels;
  const gkyl_dg_rad_gyrokinetic_surf_kern_list *surf_vpar_kernels, *surf_mu_kernels;
  const gkyl_dg_rad_gyrokinetic_boundary_surf_kern_list *boundary_surf_vpar_kernels, *boundary_surf_mu_kernels;
  
  switch (conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_vpar_kernels = ser_surf_vpar_kernels;
      surf_mu_kernels = ser_surf_mu_kernels;
      boundary_surf_vpar_kernels = ser_boundary_surf_vpar_kernels;
      boundary_surf_mu_kernels = ser_boundary_surf_mu_kernels;
      break;

    default:
      assert(false);
      break;    
  }  

  grad->eqn.vol_term = CK(vol_kernels, cdim, vdim, poly_order);

  grad->surf[0] = CK(surf_vpar_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    grad->surf[1] = CK(surf_mu_kernels, cdim, vdim, poly_order);

  grad->boundary_surf[0] = CK(boundary_surf_vpar_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    grad->boundary_surf[1] = CK(boundary_surf_mu_kernels, cdim, vdim, poly_order);
  
  // ensure non-NULL pointers
  for (int i=0; i<vdim; ++i) assert(grad->surf[i]);
  for (int i=0; i<vdim; ++i) assert(grad->boundary_surf[i]);

  grad->auxfields.nvnu_surf = 0;
  grad->auxfields.nvnu = 0;
  grad->auxfields.nvsqnu_surf = 0;
  grad->auxfields.nvsqnu = 0;
  grad->phase_range = *phase_range;
  grad->conf_range = *conf_range;
  grad->vel_map = gkyl_velocity_map_acquire(vel_map);

  grad->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(grad->eqn.flags);
  grad->eqn.ref_count = gkyl_ref_count_init(gkyl_rad_gyrokinetic_free);
  grad->eqn.on_dev = &grad->eqn;

  return &grad->eqn;
}

