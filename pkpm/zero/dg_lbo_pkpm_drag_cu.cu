/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_lbo_pkpm_drag.h>    
#include <gkyl_dg_lbo_pkpm_drag_priv.h>
}

#include <cassert>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to nuSum and nuPrimMomsSum (collision frequency * primitive moments)
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_lbo_pkpm_drag_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuSum, const struct gkyl_array *nuPrimMomsSum)
{
  struct dg_lbo_pkpm_drag *lbo_pkpm_drag = container_of(eqn, struct dg_lbo_pkpm_drag, eqn);
  lbo_pkpm_drag->auxfields.nuSum = nuSum;
  lbo_pkpm_drag->auxfields.nuPrimMomsSum = nuPrimMomsSum;
}

// Host-side wrapper for device kernels setting nuSum and nuPrimMomsSum.
void
gkyl_lbo_pkpm_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_pkpm_drag_auxfields auxin)
{
  gkyl_lbo_pkpm_drag_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.nuSum->on_dev, auxin.nuPrimMomsSum->on_dev);
}

// CUDA kernel to set device pointers to range object and Vlasov PKPM LBO drag kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_lbo_pkpm_drag_set_cu_dev_ptrs(struct dg_lbo_pkpm_drag *lbo_pkpm_drag, enum gkyl_basis_type b_type,
  int cdim, int poly_order)
{
  lbo_pkpm_drag->auxfields.nuSum = 0; 
  lbo_pkpm_drag->auxfields.nuPrimMomsSum = 0; 

  lbo_pkpm_drag->eqn.surf_term = surf;
  lbo_pkpm_drag->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_lbo_pkpm_drag_vol_kern_list *vol_kernels;
  const gkyl_dg_lbo_pkpm_drag_surf_kern_list *surf_vpar_kernels;
  const gkyl_dg_lbo_pkpm_drag_boundary_surf_kern_list *boundary_surf_vpar_kernels;
  
  switch (b_type) {
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
}

struct gkyl_dg_eqn*
gkyl_dg_lbo_pkpm_drag_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid)
{
  struct dg_lbo_pkpm_drag *lbo_pkpm_drag =
    (struct dg_lbo_pkpm_drag*) gkyl_malloc(sizeof(struct dg_lbo_pkpm_drag));

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;

  lbo_pkpm_drag->cdim = cdim;
  lbo_pkpm_drag->pdim = pdim;

  lbo_pkpm_drag->eqn.num_equations = 2;
  lbo_pkpm_drag->conf_range = *conf_range;
  lbo_pkpm_drag->vMaxSq = pow(pgrid->upper[cdim],2);
  lbo_pkpm_drag->num_cbasis = cbasis->num_basis;
  
  lbo_pkpm_drag->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(lbo_pkpm_drag->eqn.flags);
  lbo_pkpm_drag->eqn.ref_count = gkyl_ref_count_init(gkyl_lbo_pkpm_drag_free);

  // copy the host struct to device struct
  struct dg_lbo_pkpm_drag *lbo_pkpm_drag_cu =
    (struct dg_lbo_pkpm_drag*) gkyl_cu_malloc(sizeof(struct dg_lbo_pkpm_drag));

  gkyl_cu_memcpy(lbo_pkpm_drag_cu, lbo_pkpm_drag,
    sizeof(struct dg_lbo_pkpm_drag), GKYL_CU_MEMCPY_H2D);

  dg_lbo_pkpm_drag_set_cu_dev_ptrs<<<1,1>>>(lbo_pkpm_drag_cu,
    cbasis->b_type, cdim, poly_order);

  lbo_pkpm_drag->eqn.on_dev = &lbo_pkpm_drag_cu->eqn;  
  
  return &lbo_pkpm_drag->eqn;
}
