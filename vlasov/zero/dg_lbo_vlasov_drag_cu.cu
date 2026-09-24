/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_lbo_vlasov_drag.h>    
#include <gkyl_dg_lbo_vlasov_drag_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to nuSum, sum of collisionalities
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_lbo_vlasov_drag_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuSum,
  const struct gkyl_array *nuPrimMomsSum)
{
  struct dg_lbo_vlasov_drag *lbo_vlasov_drag = container_of(eqn, struct dg_lbo_vlasov_drag, eqn);
  lbo_vlasov_drag->auxfields.nuSum = nuSum;
  lbo_vlasov_drag->auxfields.nuPrimMomsSum = nuPrimMomsSum;
}

//// Host-side wrapper for device kernels setting nuSum, nuUSum and nuVtSqSum.
void
gkyl_lbo_vlasov_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_vlasov_drag_auxfields auxin)
{
  gkyl_lbo_vlasov_drag_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.nuSum->on_dev,
    auxin.nuPrimMomsSum->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov LBO kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_lbo_vlasov_drag_set_cu_dev_ptrs(struct dg_lbo_vlasov_drag *lbo_vlasov_drag, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order)
{
  lbo_vlasov_drag->auxfields.nuSum = 0; 
  lbo_vlasov_drag->auxfields.nuPrimMomsSum = 0; 

  lbo_vlasov_drag->eqn.surf_term = surf;
  lbo_vlasov_drag->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_lbo_vlasov_drag_vol_kern_list *vol_kernels;
  const gkyl_dg_lbo_vlasov_drag_surf_kern_list *surf_vx_kernels, *surf_vy_kernels, *surf_vz_kernels;
  const gkyl_dg_lbo_vlasov_drag_boundary_surf_kern_list *boundary_surf_vx_kernels, *boundary_surf_vy_kernels,
    *boundary_surf_vz_kernels;
  
  switch (b_type) {
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
 
  lbo_vlasov_drag->eqn.vol_term = vol_kernels[cv_index].kernels[poly_order];

  lbo_vlasov_drag->surf[0] = surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    lbo_vlasov_drag->surf[1] = surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    lbo_vlasov_drag->surf[2] = surf_vz_kernels[cv_index].kernels[poly_order];

  lbo_vlasov_drag->boundary_surf[0] = boundary_surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    lbo_vlasov_drag->boundary_surf[1] = boundary_surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    lbo_vlasov_drag->boundary_surf[2] = boundary_surf_vz_kernels[cv_index].kernels[poly_order];
}

struct gkyl_dg_eqn*
gkyl_dg_lbo_vlasov_drag_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid)
{
  struct dg_lbo_vlasov_drag *lbo_vlasov_drag =
    (struct dg_lbo_vlasov_drag*) gkyl_malloc(sizeof(struct dg_lbo_vlasov_drag));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  lbo_vlasov_drag->cdim = cdim;
  lbo_vlasov_drag->vdim = vdim;
  lbo_vlasov_drag->pdim = pdim;

  lbo_vlasov_drag->eqn.num_equations = 1;
  lbo_vlasov_drag->conf_range = *conf_range;

  lbo_vlasov_drag->vMaxSq = -1.;
  for (int d=0; d<vdim; d++) {
    lbo_vlasov_drag->viMax[d] = pgrid->upper[cdim+d];
    lbo_vlasov_drag->vMaxSq = fmax(lbo_vlasov_drag->vMaxSq, pow(pgrid->upper[cdim+d],2));
  }
  lbo_vlasov_drag->num_cbasis = cbasis->num_basis;

  lbo_vlasov_drag->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(lbo_vlasov_drag->eqn.flags);
  lbo_vlasov_drag->eqn.ref_count = gkyl_ref_count_init(gkyl_lbo_vlasov_drag_free);

  // copy the host struct to device struct
  struct dg_lbo_vlasov_drag *lbo_vlasov_drag_cu =
    (struct dg_lbo_vlasov_drag*) gkyl_cu_malloc(sizeof(struct dg_lbo_vlasov_drag));

  gkyl_cu_memcpy(lbo_vlasov_drag_cu, lbo_vlasov_drag,
    sizeof(struct dg_lbo_vlasov_drag), GKYL_CU_MEMCPY_H2D);

  dg_lbo_vlasov_drag_set_cu_dev_ptrs<<<1,1>>>(lbo_vlasov_drag_cu,
    cbasis->b_type, cv_index[cdim].vdim[vdim], cdim, vdim, poly_order);

  lbo_vlasov_drag->eqn.on_dev = &lbo_vlasov_drag_cu->eqn;  
  
  return &lbo_vlasov_drag->eqn;
}
