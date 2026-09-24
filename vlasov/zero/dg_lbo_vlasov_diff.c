#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_lbo_vlasov_diff.h>
#include <gkyl_dg_lbo_vlasov_diff_priv.h>
#include <gkyl_util.h>

void
gkyl_lbo_vlasov_diff_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff  = container_of(base, struct dg_lbo_vlasov_diff, eqn);

  if (GKYL_IS_CU_ALLOC(lbo_vlasov_diff->eqn.flags))
    gkyl_cu_free(lbo_vlasov_diff->eqn.on_dev);
  
  gkyl_free(lbo_vlasov_diff);
}

void
gkyl_lbo_vlasov_diff_set_auxfields(const struct gkyl_dg_eqn *eqn, const struct gkyl_dg_lbo_vlasov_diff_auxfields auxin)
{

#ifdef GKYL_HAVE_CUDA
 if (gkyl_array_is_cu_dev(auxin.nuSum) && gkyl_array_is_cu_dev(auxin.nuPrimMomsSum)) {
   gkyl_lbo_vlasov_diff_set_auxfields_cu(eqn->on_dev, auxin);
   return;
 }
#endif

  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  lbo_vlasov_diff->auxfields.nuSum = auxin.nuSum;
  lbo_vlasov_diff->auxfields.nuPrimMomsSum = auxin.nuPrimMomsSum;
}


struct gkyl_dg_eqn*
gkyl_dg_lbo_vlasov_diff_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_lbo_vlasov_diff_cu_dev_new(cbasis, pbasis, conf_range, pgrid);
  } 
#endif
  struct dg_lbo_vlasov_diff* lbo_vlasov_diff = gkyl_malloc(sizeof(struct dg_lbo_vlasov_diff));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  lbo_vlasov_diff->cdim = cdim;
  lbo_vlasov_diff->vdim = vdim;
  lbo_vlasov_diff->pdim = pdim;

  lbo_vlasov_diff->eqn.num_equations = 1;
  lbo_vlasov_diff->eqn.surf_term = surf;
  lbo_vlasov_diff->eqn.boundary_surf_term = boundary_surf;

  lbo_vlasov_diff->vMaxSq = -1.;
  for (int d=0; d<vdim; d++) {
    lbo_vlasov_diff->viMax[d] = pgrid->upper[cdim+d];
    lbo_vlasov_diff->vMaxSq = fmax(lbo_vlasov_diff->vMaxSq, pow(pgrid->upper[cdim+d],2));
  }
  lbo_vlasov_diff->num_cbasis = cbasis->num_basis;

  const gkyl_dg_lbo_vlasov_diff_vol_kern_list *vol_kernels;
  const gkyl_dg_lbo_vlasov_diff_surf_kern_list *surf_vx_kernels, *surf_vy_kernels, *surf_vz_kernels;
  const gkyl_dg_lbo_vlasov_diff_boundary_surf_kern_list *boundary_surf_vx_kernels, *boundary_surf_vy_kernels,
    *boundary_surf_vz_kernels;
  
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_vx_kernels = ser_surf_vx_kernels;
      surf_vy_kernels = ser_surf_vy_kernels;
      surf_vz_kernels = ser_surf_vz_kernels;
      boundary_surf_vx_kernels = ser_boundary_surf_vx_kernels;
      boundary_surf_vy_kernels = ser_boundary_surf_vy_kernels;
      boundary_surf_vz_kernels = ser_boundary_surf_vz_kernels;
      break;

    default:
      assert(false);
      break;    
  }  

  lbo_vlasov_diff->eqn.vol_term = CK(vol_kernels, cdim, vdim, poly_order);

  lbo_vlasov_diff->surf[0] = CK(surf_vx_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    lbo_vlasov_diff->surf[1] = CK(surf_vy_kernels, cdim, vdim, poly_order);
  if (vdim>2)
    lbo_vlasov_diff->surf[2] = CK(surf_vz_kernels, cdim, vdim, poly_order);

  lbo_vlasov_diff->boundary_surf[0] = CK(boundary_surf_vx_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    lbo_vlasov_diff->boundary_surf[1] = CK(boundary_surf_vy_kernels, cdim, vdim, poly_order);
  if (vdim>2)
    lbo_vlasov_diff->boundary_surf[2] = CK(boundary_surf_vz_kernels, cdim, vdim, poly_order);

  // ensure non-NULL pointers
  for (int i=0; i<vdim; ++i) assert(lbo_vlasov_diff->surf[i]);
  for (int i=0; i<vdim; ++i) assert(lbo_vlasov_diff->boundary_surf[i]);

  lbo_vlasov_diff->auxfields.nuSum = 0;
  lbo_vlasov_diff->auxfields.nuPrimMomsSum = 0;
  lbo_vlasov_diff->conf_range = *conf_range;

  lbo_vlasov_diff->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(lbo_vlasov_diff->eqn.flags);
  lbo_vlasov_diff->eqn.ref_count = gkyl_ref_count_init(gkyl_lbo_vlasov_diff_free);
  lbo_vlasov_diff->eqn.on_dev = &lbo_vlasov_diff->eqn;
  
  return &lbo_vlasov_diff->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_lbo_vlasov_diff_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid)
{
  assert(false);
  return 0;
}

#endif
