#include <assert.h>
#include <stdio.h>
#include <float.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_lbo_gyrokinetic_drag.h>
#include <gkyl_dg_lbo_gyrokinetic_drag_priv.h>
#include <gkyl_util.h>

void
gkyl_lbo_gyrokinetic_drag_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_lbo_gyrokinetic_drag *lbo = container_of(base, struct dg_lbo_gyrokinetic_drag, eqn);
  gkyl_gk_geometry_release(lbo->gk_geom);
  gkyl_velocity_map_release(lbo->vel_map);

  if (GKYL_IS_CU_ALLOC(lbo->eqn.flags))
    gkyl_cu_free(lbo->eqn.on_dev);
  
  gkyl_free(lbo);
}

void
gkyl_lbo_gyrokinetic_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_gyrokinetic_drag_auxfields auxin)
{

#ifdef GKYL_HAVE_CUDA
 if (gkyl_array_is_cu_dev(auxin.nuSum) &&
     gkyl_array_is_cu_dev(auxin.nuPrimMomsSum) && gkyl_array_is_cu_dev(auxin.m2self)) {
   gkyl_lbo_gyrokinetic_drag_set_auxfields_cu(eqn->on_dev, auxin);
   return;
 }
#endif

  struct dg_lbo_gyrokinetic_drag *lbo = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  lbo->auxfields.nuSum = auxin.nuSum;
  lbo->auxfields.nuPrimMomsSum = auxin.nuPrimMomsSum;
  lbo->auxfields.m2self = auxin.m2self;
}


struct gkyl_dg_eqn*
gkyl_dg_lbo_gyrokinetic_drag_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range,  const struct gkyl_rect_grid *pgrid,
  double mass, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map, 
  bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_dg_lbo_gyrokinetic_drag_cu_dev_new(cbasis, pbasis,
      conf_range, pgrid, mass, gk_geom, vel_map);
#endif
  struct dg_lbo_gyrokinetic_drag* lbo = gkyl_malloc(sizeof(struct dg_lbo_gyrokinetic_drag));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  lbo->cdim = cdim;
  lbo->pdim = pdim;

  lbo->eqn.num_equations = 1;
  lbo->eqn.surf_term = surf;
  lbo->eqn.boundary_surf_term = boundary_surf;

  lbo->vparMax = GKYL_MAX2(fabs(vel_map->vbounds[0]),vel_map->vbounds[vdim]);
  lbo->vparMaxSq = pow(lbo->vparMax,2);
  lbo->num_cbasis = cbasis->num_basis;

  const gkyl_dg_lbo_gyrokinetic_drag_vol_kern_list *vol_kernels;
  const gkyl_dg_lbo_gyrokinetic_drag_surf_kern_list *surf_vpar_kernels, *surf_mu_kernels;
  const gkyl_dg_lbo_gyrokinetic_drag_boundary_surf_kern_list *boundary_surf_vpar_kernels, *boundary_surf_mu_kernels;
  
  switch (cbasis->b_type) {
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

  lbo->eqn.vol_term = CK(vol_kernels, cdim, vdim, poly_order);

  lbo->surf[0] = CK(surf_vpar_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    lbo->surf[1] = CK(surf_mu_kernels, cdim, vdim, poly_order);

  lbo->boundary_surf[0] = CK(boundary_surf_vpar_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    lbo->boundary_surf[1] = CK(boundary_surf_mu_kernels, cdim, vdim, poly_order);

  // ensure non-NULL pointers
  for (int i=0; i<vdim; ++i) assert(lbo->surf[i]);
  for (int i=0; i<vdim; ++i) assert(lbo->boundary_surf[i]);

  lbo->mass = mass;
  lbo->conf_range = *conf_range;
  lbo->gk_geom = gkyl_gk_geometry_acquire(gk_geom);
  lbo->vel_map = gkyl_velocity_map_acquire(vel_map);
  
  lbo->auxfields.nuSum = 0;
  lbo->auxfields.nuPrimMomsSum = 0;
  lbo->auxfields.m2self = 0;

  lbo->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(lbo->eqn.flags);
  lbo->eqn.ref_count = gkyl_ref_count_init(gkyl_lbo_gyrokinetic_drag_free);
  lbo->eqn.on_dev = &lbo->eqn;
  
  return &lbo->eqn;
}
