/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_fpo_vlasov_diff.h>    
#include <gkyl_dg_fpo_vlasov_diff_priv.h>
}

#include <cassert>

// "Choose Kernel" based on cdim and polynomial order
#define CK(lst, cdim, poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to g (second Rosenbluth potential).
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_fpo_vlasov_diff_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *g)
{
  struct dg_fpo_vlasov_diff *fpo_vlasov_diff = container_of(eqn, struct dg_fpo_vlasov_diff, eqn);
  fpo_vlasov_diff->auxfields.g = g;
}

//// Host-side wrapper for device kernels setting g (second Rosenbluth potential).
void
gkyl_fpo_vlasov_diff_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_fpo_vlasov_diff_auxfields auxin)
{
  gkyl_fpo_vlasov_diff_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.g->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov fpo kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_fpo_vlasov_diff_set_cu_dev_ptrs(struct dg_fpo_vlasov_diff *fpo_vlasov_diff, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  fpo_vlasov_diff->auxfields.g = 0; 

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
  
  switch (b_type) {
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
}

struct gkyl_dg_eqn*
gkyl_dg_fpo_vlasov_diff_cu_dev_new(const struct gkyl_basis* pbasis, const struct gkyl_range* phase_range)
{
  struct dg_fpo_vlasov_diff *fpo_vlasov_diff =
    (struct dg_fpo_vlasov_diff*) gkyl_malloc(sizeof(struct dg_fpo_vlasov_diff));

  // Vlasov Fokker-Planck operator only defined in 3 velocity dimensions
  int pdim = pbasis->ndim, vdim = 3, cdim = pdim - vdim;
  int poly_order = pbasis->poly_order;

  fpo_vlasov_diff->cdim = cdim;
  fpo_vlasov_diff->pdim = pdim;

  fpo_vlasov_diff->eqn.num_equations = 1;
  fpo_vlasov_diff->phase_range = *phase_range;

  fpo_vlasov_diff->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(fpo_vlasov_diff->eqn.flags);
  fpo_vlasov_diff->eqn.ref_count = gkyl_ref_count_init(gkyl_fpo_vlasov_diff_free);

  // copy the host struct to device struct
  struct dg_fpo_vlasov_diff *fpo_vlasov_diff_cu =
    (struct dg_fpo_vlasov_diff*) gkyl_cu_malloc(sizeof(struct dg_fpo_vlasov_diff));

  gkyl_cu_memcpy(fpo_vlasov_diff_cu, fpo_vlasov_diff,
    sizeof(struct dg_fpo_vlasov_diff), GKYL_CU_MEMCPY_H2D);

  dg_fpo_vlasov_diff_set_cu_dev_ptrs<<<1,1>>>(fpo_vlasov_diff_cu,
    pbasis->b_type, cdim, poly_order);

  fpo_vlasov_diff->eqn.on_dev = &fpo_vlasov_diff_cu->eqn;
  
  return &fpo_vlasov_diff->eqn;
}
