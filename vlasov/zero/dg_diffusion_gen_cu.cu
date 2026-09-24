/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_diffusion_gen.h>    
#include <gkyl_dg_diffusion_gen_priv.h>
}

#include <cassert>

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_diffusion_gen_set_auxfields_cu_kernel(const struct gkyl_dg_eqn* eqn, const struct gkyl_array* Dij)
{
  struct dg_diffusion_gen* diffusion_gen = container_of(eqn, struct dg_diffusion_gen, eqn);
  diffusion_gen->auxfields.Dij = Dij;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_diffusion_gen_set_auxfields_cu(const struct gkyl_dg_eqn* eqn, struct gkyl_dg_diffusion_gen_auxfields auxin)
{
  gkyl_diffusion_gen_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.Dij->on_dev);
}

__global__ void static
dg_diffusion_gen_set_cu_dev_ptrs(struct dg_diffusion_gen* diffusion_gen, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  diffusion_gen->auxfields.Dij = 0; 

  const gkyl_dg_diffusion_gen_vol_kern_list* vol_kernels;
  const gkyl_dg_diffusion_gen_surf_kern_list* surf_xx_kernels;
  const gkyl_dg_diffusion_gen_surf_kern_list* surf_xy_kernels;
  const gkyl_dg_diffusion_gen_surf_kern_list* surf_xz_kernels;
  const gkyl_dg_diffusion_gen_surf_kern_list* surf_yx_kernels;
  const gkyl_dg_diffusion_gen_surf_kern_list* surf_yy_kernels;
  const gkyl_dg_diffusion_gen_surf_kern_list* surf_yz_kernels;
  const gkyl_dg_diffusion_gen_surf_kern_list* surf_zx_kernels;
  const gkyl_dg_diffusion_gen_surf_kern_list* surf_zy_kernels;
  const gkyl_dg_diffusion_gen_surf_kern_list* surf_zz_kernels; 

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
      break;

    default:
      assert(false);
      break;    
  } 
  
  diffusion_gen->eqn.gen_surf_term = surf;
  
  diffusion_gen->eqn.vol_term = CK(vol_kernels, cdim, poly_order);

  diffusion_gen->surf[0][0] = CK(surf_xx_kernels, cdim, poly_order);
  if (cdim>1) {
    diffusion_gen->surf[0][1] = CK(surf_xy_kernels, cdim, poly_order);
    diffusion_gen->surf[1][0] = CK(surf_yx_kernels, cdim, poly_order);
    diffusion_gen->surf[1][1] = CK(surf_yy_kernels, cdim, poly_order);
  }
  if (cdim>2) {
    diffusion_gen->surf[0][2] = CK(surf_xz_kernels, cdim, poly_order);
    diffusion_gen->surf[1][2] = CK(surf_yz_kernels, cdim, poly_order);
    diffusion_gen->surf[2][0] = CK(surf_zx_kernels, cdim, poly_order);
    diffusion_gen->surf[2][1] = CK(surf_zy_kernels, cdim, poly_order);
    diffusion_gen->surf[2][2] = CK(surf_zz_kernels, cdim, poly_order);
  }
}

struct gkyl_dg_eqn*
gkyl_dg_diffusion_gen_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range)
{
  struct dg_diffusion_gen* diffusion_gen = (struct dg_diffusion_gen*) gkyl_malloc(sizeof(struct dg_diffusion_gen));

  // set basic parameters
  diffusion_gen->eqn.num_equations = 1;
  diffusion_gen->conf_range = *conf_range;

  diffusion_gen->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(diffusion_gen->eqn.flags);
  diffusion_gen->eqn.ref_count = gkyl_ref_count_init(gkyl_diffusion_gen_free);

  // copy the host struct to device struct
  struct dg_diffusion_gen* diffusion_gen_cu = (struct dg_diffusion_gen*) gkyl_cu_malloc(sizeof(struct dg_diffusion_gen));
  gkyl_cu_memcpy(diffusion_gen_cu, diffusion_gen, sizeof(struct dg_diffusion_gen), GKYL_CU_MEMCPY_H2D);
  dg_diffusion_gen_set_cu_dev_ptrs<<<1,1>>>(diffusion_gen_cu, cbasis->b_type, cbasis->ndim, cbasis->poly_order);

  // set parent on_dev pointer
  diffusion_gen->eqn.on_dev = &diffusion_gen_cu->eqn;

  return &diffusion_gen->eqn;
}
