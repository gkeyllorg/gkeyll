/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_fpo_vlasov_drag.h>    
#include <gkyl_dg_fpo_vlasov_drag_priv.h>
}

#include <cassert>

// "Choose Kernel" based on cdim and polynomial order
#define CK(lst, cdim, poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to drag coefficient.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_fpo_vlasov_drag_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn,
  const struct gkyl_array *drag_coeff, const struct gkyl_array *drag_coeff_surf,
  const struct gkyl_array *sgn_drag_coeff_surf, const struct gkyl_array *const_sgn_drag_coeff_surf)
{
  struct dg_fpo_vlasov_drag *fpo_vlasov_drag = container_of(eqn, struct dg_fpo_vlasov_drag, eqn);
  fpo_vlasov_drag->auxfields.drag_coeff = drag_coeff;
  fpo_vlasov_drag->auxfields.drag_coeff_surf = drag_coeff_surf;
  fpo_vlasov_drag->auxfields.sgn_drag_coeff_surf = sgn_drag_coeff_surf;
  fpo_vlasov_drag->auxfields.const_sgn_drag_coeff_surf = const_sgn_drag_coeff_surf;
}

//// Host-side wrapper for device kernels setting g (second Rosenbluth potential).
void
gkyl_fpo_vlasov_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_fpo_vlasov_drag_auxfields auxin)
{
  gkyl_fpo_vlasov_drag_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.drag_coeff->on_dev, auxin.drag_coeff_surf->on_dev,
    auxin.sgn_drag_coeff_surf->on_dev, auxin.const_sgn_drag_coeff_surf->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov fpo kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_fpo_vlasov_drag_set_cu_dev_ptrs(struct dg_fpo_vlasov_drag *fpo_vlasov_drag, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  fpo_vlasov_drag->auxfields.drag_coeff = 0; 
  fpo_vlasov_drag->auxfields.drag_coeff_surf = 0;
  fpo_vlasov_drag->auxfields.sgn_drag_coeff_surf = 0;
  fpo_vlasov_drag->auxfields.const_sgn_drag_coeff_surf = 0;

  fpo_vlasov_drag->eqn.surf_term = surf;
  fpo_vlasov_drag->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_fpo_vlasov_drag_vol_kern_list *vol_kernels;
  const gkyl_dg_fpo_vlasov_drag_surf_kern_list *surf_vx_kernel_list;
  const gkyl_dg_fpo_vlasov_drag_surf_kern_list *surf_vy_kernel_list;
  const gkyl_dg_fpo_vlasov_drag_surf_kern_list *surf_vz_kernel_list;
  const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list *boundary_surf_vx_kernel_list;
  const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list *boundary_surf_vy_kernel_list;
  const gkyl_dg_fpo_vlasov_drag_boundary_surf_kern_list *boundary_surf_vz_kernel_list;
  
  switch (b_type) {
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
}

struct gkyl_dg_eqn*
gkyl_dg_fpo_vlasov_drag_cu_dev_new(const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range)
{
  struct dg_fpo_vlasov_drag *fpo_vlasov_drag =
    (struct dg_fpo_vlasov_drag*) gkyl_malloc(sizeof(struct dg_fpo_vlasov_drag));

  // Vlasov Fokker-Planck operator only defined in 3 velocity dimensions
  int pdim = pbasis->ndim, vdim = 3, cdim = pdim - vdim;
  int poly_order = pbasis->poly_order;

  fpo_vlasov_drag->cdim = cdim;
  fpo_vlasov_drag->pdim = pdim;

  fpo_vlasov_drag->eqn.num_equations = 1;
  fpo_vlasov_drag->phase_range = *phase_range;

  fpo_vlasov_drag->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(fpo_vlasov_drag->eqn.flags);
  fpo_vlasov_drag->eqn.ref_count = gkyl_ref_count_init(gkyl_fpo_vlasov_drag_free);

  // copy the host struct to device struct
  struct dg_fpo_vlasov_drag *fpo_vlasov_drag_cu =
    (struct dg_fpo_vlasov_drag*) gkyl_cu_malloc(sizeof(struct dg_fpo_vlasov_drag));

  gkyl_cu_memcpy(fpo_vlasov_drag_cu, fpo_vlasov_drag,
    sizeof(struct dg_fpo_vlasov_drag), GKYL_CU_MEMCPY_H2D);

  dg_fpo_vlasov_drag_set_cu_dev_ptrs<<<1,1>>>(fpo_vlasov_drag_cu,
    pbasis->b_type, cdim, poly_order);

  fpo_vlasov_drag->eqn.on_dev = &fpo_vlasov_drag_cu->eqn;
  
  return &fpo_vlasov_drag->eqn;
}
