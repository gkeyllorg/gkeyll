/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_fpo_vlasov_diff.h>    
#include <gkyl_dg_fpo_vlasov_diff_priv.h>
}

#include <cassert>

// "Choose Kernel" based on cdim and polynomial order
#define CK(lst, cdim, poly_order) lst[cdim-1].kernels[poly_order-1]

// CUDA kernel to set pointer to diffusion tensor.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_fpo_vlasov_diff_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *diff_coeff, const struct gkyl_array *diff_coeff_surf)
{
  struct dg_fpo_vlasov_diff *fpo_vlasov_diff = container_of(eqn, struct dg_fpo_vlasov_diff, eqn);
  fpo_vlasov_diff->auxfields.diff_coeff = diff_coeff;
  fpo_vlasov_diff->auxfields.diff_coeff_surf = diff_coeff_surf;
}

//// Host-side wrapper for device kernels setting g (second Rosenbluth potential).
void
gkyl_fpo_vlasov_diff_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_fpo_vlasov_diff_auxfields auxin)
{
  gkyl_fpo_vlasov_diff_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.diff_coeff->on_dev, auxin.diff_coeff_surf->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov fpo kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_fpo_vlasov_diff_set_cu_dev_ptrs(struct dg_fpo_vlasov_diff *fpo_vlasov_diff, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  fpo_vlasov_diff->auxfields.diff_coeff = 0; 

  fpo_vlasov_diff->eqn.gen_surf_term = fpo_diff_gen_surf_term;

  const gkyl_dg_fpo_vlasov_diff_vol_kern_list* vol_kernels;
  const fpo_vlasov_diff_surf_stencil_list* surf_vxvx_kernel_list;
  const fpo_vlasov_diff_surf_stencil_list* surf_vxvy_kernel_list;
  const fpo_vlasov_diff_surf_stencil_list* surf_vxvz_kernel_list;
  const fpo_vlasov_diff_surf_stencil_list* surf_vyvx_kernel_list;
  const fpo_vlasov_diff_surf_stencil_list* surf_vyvy_kernel_list;
  const fpo_vlasov_diff_surf_stencil_list* surf_vyvz_kernel_list;
  const fpo_vlasov_diff_surf_stencil_list* surf_vzvx_kernel_list;
  const fpo_vlasov_diff_surf_stencil_list* surf_vzvy_kernel_list;
  const fpo_vlasov_diff_surf_stencil_list* surf_vzvz_kernel_list;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      assert(poly_order == 2);
      vol_kernels = ser_vol_kernels;
      surf_vxvx_kernel_list = ser_surf_vxvx_kernels;
      surf_vxvy_kernel_list = ser_surf_vxvy_kernels;
      surf_vxvz_kernel_list = ser_surf_vxvz_kernels;
      surf_vyvx_kernel_list = ser_surf_vyvx_kernels;
      surf_vyvy_kernel_list = ser_surf_vyvy_kernels;
      surf_vyvz_kernel_list = ser_surf_vyvz_kernels;
      surf_vzvx_kernel_list = ser_surf_vzvx_kernels;
      surf_vzvy_kernel_list = ser_surf_vzvy_kernels;
      surf_vzvz_kernel_list = ser_surf_vzvz_kernels;
      break;

    case GKYL_BASIS_MODAL_HYBRID:
      assert(poly_order == 1);
      vol_kernels = ser_vol_kernels;
      surf_vxvx_kernel_list = ser_surf_vxvx_kernels;
      surf_vxvy_kernel_list = ser_surf_vxvy_kernels;
      surf_vxvz_kernel_list = ser_surf_vxvz_kernels;
      surf_vyvx_kernel_list = ser_surf_vyvx_kernels;
      surf_vyvy_kernel_list = ser_surf_vyvy_kernels;
      surf_vyvz_kernel_list = ser_surf_vyvz_kernels;
      surf_vzvx_kernel_list = ser_surf_vzvx_kernels;
      surf_vzvy_kernel_list = ser_surf_vzvy_kernels;
      surf_vzvz_kernel_list = ser_surf_vzvz_kernels;
      break;

    default:
      assert(false);
      break;
  }  
 
  fpo_vlasov_diff->eqn.vol_term = CK(vol_kernels, cdim, poly_order);

  fpo_vlasov_diff->surf[0][0] = surf_vxvx_kernel_list[cdim-1].list[poly_order-1];
  fpo_vlasov_diff->surf[0][1] = surf_vxvy_kernel_list[cdim-1].list[poly_order-1];
  fpo_vlasov_diff->surf[0][2] = surf_vxvz_kernel_list[cdim-1].list[poly_order-1];
  fpo_vlasov_diff->surf[1][0] = surf_vyvx_kernel_list[cdim-1].list[poly_order-1];
  fpo_vlasov_diff->surf[1][1] = surf_vyvy_kernel_list[cdim-1].list[poly_order-1];
  fpo_vlasov_diff->surf[1][2] = surf_vyvz_kernel_list[cdim-1].list[poly_order-1];
  fpo_vlasov_diff->surf[2][0] = surf_vzvx_kernel_list[cdim-1].list[poly_order-1];
  fpo_vlasov_diff->surf[2][1] = surf_vzvy_kernel_list[cdim-1].list[poly_order-1];
  fpo_vlasov_diff->surf[2][2] = surf_vzvz_kernel_list[cdim-1].list[poly_order-1];
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
