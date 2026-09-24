#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_lbo_pkpm_drag.h>
#include <gkyl_dg_lbo_pkpm_drag_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_lbo_pkpm_drag_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_lbo_pkpm_drag *lbo_pkpm_drag  = container_of(base, struct dg_lbo_pkpm_drag, eqn);

  if (GKYL_IS_CU_ALLOC(lbo_pkpm_drag->eqn.flags))
    gkyl_cu_free(lbo_pkpm_drag->eqn.on_dev);
  
  gkyl_free(lbo_pkpm_drag);
}

void
gkyl_lbo_pkpm_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, const struct gkyl_dg_lbo_pkpm_drag_auxfields auxin)
{

#ifdef GKYL_HAVE_CUDA
 if (gkyl_array_is_cu_dev(auxin.nuSum) && gkyl_array_is_cu_dev(auxin.nuPrimMomsSum)) {
   gkyl_lbo_pkpm_drag_set_auxfields_cu(eqn->on_dev, auxin);
   return;
 }
#endif

  struct dg_lbo_pkpm_drag *lbo_pkpm_drag = container_of(eqn, struct dg_lbo_pkpm_drag, eqn);
  lbo_pkpm_drag->auxfields.nuSum = auxin.nuSum;
  lbo_pkpm_drag->auxfields.nuPrimMomsSum = auxin.nuPrimMomsSum;
}

struct gkyl_dg_eqn*
gkyl_dg_lbo_pkpm_drag_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_lbo_pkpm_drag_cu_dev_new(cbasis, pbasis, conf_range, pgrid);
  } 
#endif
  struct dg_lbo_pkpm_drag* lbo_pkpm_drag = gkyl_malloc(sizeof(struct dg_lbo_pkpm_drag));

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;

  lbo_pkpm_drag->cdim = cdim;
  lbo_pkpm_drag->pdim = pdim;

  lbo_pkpm_drag->eqn.num_equations = 2;
  lbo_pkpm_drag->eqn.surf_term = surf;
  lbo_pkpm_drag->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_lbo_pkpm_drag_vol_kern_list *vol_kernels;
  const gkyl_dg_lbo_pkpm_drag_surf_kern_list *surf_vpar_kernels;
  const gkyl_dg_lbo_pkpm_drag_boundary_surf_kern_list *boundary_surf_vpar_kernels;
  
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_vpar_kernels = ser_surf_vpar_kernels;
      boundary_surf_vpar_kernels = ser_boundary_surf_vpar_kernels;
      
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = ten_vol_kernels;
      surf_vpar_kernels = ten_surf_vpar_kernels;
      boundary_surf_vpar_kernels = ten_boundary_surf_vpar_kernels;
      
      break;

    default:
      assert(false);
      break;    
  }  

  lbo_pkpm_drag->eqn.vol_term = CK(vol_kernels, cdim, poly_order);

  lbo_pkpm_drag->surf = CK(surf_vpar_kernels, cdim, poly_order);

  lbo_pkpm_drag->boundary_surf = CK(boundary_surf_vpar_kernels, cdim, poly_order);

  // ensure non-NULL pointers
  assert(lbo_pkpm_drag->surf);
  assert(lbo_pkpm_drag->boundary_surf);

  lbo_pkpm_drag->auxfields.nuSum = 0;
  lbo_pkpm_drag->auxfields.nuPrimMomsSum = 0;
  lbo_pkpm_drag->conf_range = *conf_range;
  lbo_pkpm_drag->vMaxSq = pow(pgrid->upper[cdim],2);
  lbo_pkpm_drag->num_cbasis = cbasis->num_basis;
  
  lbo_pkpm_drag->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(lbo_pkpm_drag->eqn.flags);
  lbo_pkpm_drag->eqn.ref_count = gkyl_ref_count_init(gkyl_lbo_pkpm_drag_free);
  lbo_pkpm_drag->eqn.on_dev = &lbo_pkpm_drag->eqn;
  
  return &lbo_pkpm_drag->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_lbo_pkpm_drag_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid)
{
  assert(false);
  return 0;
}

#endif
