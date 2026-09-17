#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_fpo_vlasov_drag.h>
#include <gkyl_dg_fpo_vlasov_drag_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polynomial order
#define CK(lst, cdim, poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_fpo_vlasov_drag_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_fpo_vlasov_drag *fpo_vlasov_drag  = container_of(base, struct dg_fpo_vlasov_drag, eqn);

  if (GKYL_IS_CU_ALLOC(fpo_vlasov_drag->eqn.flags))
    gkyl_cu_free(fpo_vlasov_drag->eqn.on_dev);
 
  gkyl_free(fpo_vlasov_drag);
}

void
gkyl_fpo_vlasov_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, const struct gkyl_dg_fpo_vlasov_drag_auxfields auxin)
{

#ifdef GKYL_HAVE_CUDA
 if (gkyl_array_is_cu_dev(auxin.drag_coeff)) {
   gkyl_fpo_vlasov_drag_set_auxfields_cu(eqn->on_dev, auxin);
   return;
 }
#endif

  struct dg_fpo_vlasov_drag *fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  fpo_vlasov_drag->auxfields.drag_coeff = auxin.drag_coeff;
  fpo_vlasov_drag->auxfields.drag_coeff_surf = auxin.drag_coeff_surf;
  fpo_vlasov_drag->auxfields.sgn_drag_coeff_surf = auxin.sgn_drag_coeff_surf;
  fpo_vlasov_drag->auxfields.const_sgn_drag_coeff_surf = auxin.const_sgn_drag_coeff_surf;
}

struct gkyl_dg_eqn*
gkyl_dg_fpo_vlasov_drag_new(const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu)
    return gkyl_dg_fpo_vlasov_drag_cu_dev_new(pbasis, phase_range);
#endif

  struct dg_fpo_vlasov_drag* fpo_vlasov_drag = gkyl_malloc(sizeof(struct dg_fpo_vlasov_drag));

  // Vlasov Fokker-Planck operator only defined in 3 velocity dimensions
  int pdim = pbasis->ndim, vdim = 3, cdim = pdim - vdim;
  int poly_order = pbasis->poly_order;

  fpo_vlasov_drag->cdim = cdim;
  fpo_vlasov_drag->pdim = pdim;

  fpo_vlasov_drag->eqn.num_equations = 1;
  fpo_vlasov_drag->eqn.surf_term = surf;
  fpo_vlasov_drag->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_fpo_vlasov_drag_vol_kern_list *vol_kernels;
  const gkyl_dg_fpo_vlasov_drag_surf_kern_list *surf_vx_kernel_list;
  const gkyl_dg_fpo_vlasov_drag_surf_kern_list *surf_vy_kernel_list;
  const gkyl_dg_fpo_vlasov_drag_surf_kern_list *surf_vz_kernel_list;
  const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list *boundary_surf_vx_kernel_list;
  const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list *boundary_surf_vy_kernel_list;
  const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list *boundary_surf_vz_kernel_list;
 
  switch (pbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      assert(poly_order == 2);
      vol_kernels = ser_vol_kernels;
      surf_vx_kernel_list = ser_surf_vx_kernels;
      surf_vy_kernel_list = ser_surf_vy_kernels;
      surf_vz_kernel_list = ser_surf_vz_kernels;
      boundary_surf_vx_kernel_list = ser_boundary_surf_vx_kernels;
      boundary_surf_vy_kernel_list = ser_boundary_surf_vy_kernels;
      boundary_surf_vz_kernel_list = ser_boundary_surf_vz_kernels;
      break;

    case GKYL_BASIS_MODAL_HYBRID:
      assert(poly_order == 1);
      vol_kernels = ser_vol_kernels;
      surf_vx_kernel_list = ser_surf_vx_kernels;
      surf_vy_kernel_list = ser_surf_vy_kernels;
      surf_vz_kernel_list = ser_surf_vz_kernels;
      boundary_surf_vx_kernel_list = ser_boundary_surf_vx_kernels;
      boundary_surf_vy_kernel_list = ser_boundary_surf_vy_kernels;
      boundary_surf_vz_kernel_list = ser_boundary_surf_vz_kernels;
      break;

    default:
      assert(false);
      break;
  }

  fpo_vlasov_drag->eqn.vol_term = CK(vol_kernels, cdim, poly_order);
  fpo_vlasov_drag->surf[0] = CK(surf_vx_kernel_list, cdim, poly_order);
  fpo_vlasov_drag->surf[1] = CK(surf_vy_kernel_list, cdim, poly_order);
  fpo_vlasov_drag->surf[2] = CK(surf_vz_kernel_list, cdim, poly_order);
  fpo_vlasov_drag->boundary_surf[0] = CK(boundary_surf_vx_kernel_list, cdim, poly_order);
  fpo_vlasov_drag->boundary_surf[1] = CK(boundary_surf_vy_kernel_list, cdim, poly_order);
  fpo_vlasov_drag->boundary_surf[2] = CK(boundary_surf_vz_kernel_list, cdim, poly_order);

  // ensure non-NULL pointers
  for (int i=0; i<vdim; ++i) assert(fpo_vlasov_drag->surf[i]);
  for (int i=0; i<vdim; ++i) assert(fpo_vlasov_drag->boundary_surf[i]);

  fpo_vlasov_drag->auxfields.drag_coeff = 0;
  fpo_vlasov_drag->auxfields.drag_coeff_surf = 0;
  fpo_vlasov_drag->auxfields.sgn_drag_coeff_surf = 0;
  fpo_vlasov_drag->auxfields.const_sgn_drag_coeff_surf = 0;
  fpo_vlasov_drag->phase_range = *phase_range;

  fpo_vlasov_drag->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(fpo_vlasov_drag->eqn.flags);
  fpo_vlasov_drag->eqn.ref_count = gkyl_ref_count_init(gkyl_fpo_vlasov_drag_free);
  fpo_vlasov_drag->eqn.on_dev = &fpo_vlasov_drag->eqn;
 
  return &fpo_vlasov_drag->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_fpo_vlasov_drag_cu_dev_new(const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range)
{
  assert(false);
  return 0;
}

#endif

