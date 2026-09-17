#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_diffusion_gen.h>
#include <gkyl_dg_diffusion_gen_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polynomial order
#define CK(lst, cdim, poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_diffusion_gen_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_diffusion_gen* diffusion_gen = container_of(base->on_dev, struct dg_diffusion_gen, eqn);
    gkyl_cu_free(diffusion_gen);
  }
  
  struct dg_diffusion_gen* diffusion_gen = container_of(base, struct dg_diffusion_gen, eqn);
  gkyl_free(diffusion_gen);
}

void
gkyl_diffusion_gen_set_auxfields(const struct gkyl_dg_eqn* eqn, struct gkyl_dg_diffusion_gen_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.Dij)) {
    gkyl_diffusion_gen_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif
  
  struct dg_diffusion_gen* diffusion_gen = container_of(eqn, struct dg_diffusion_gen, eqn);
  diffusion_gen->auxfields.Dij = auxin.Dij;
}

struct gkyl_dg_eqn*
gkyl_dg_diffusion_gen_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu)
    return gkyl_dg_diffusion_gen_cu_dev_new(cbasis, conf_range);
#endif
  
  struct dg_diffusion_gen* diffusion_gen = gkyl_malloc(sizeof(struct dg_diffusion_gen));

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

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

  switch (cbasis->b_type) {
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

  diffusion_gen->eqn.num_equations = 1;
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

  // ensure non-NULL pointers
  for (int i=0; i<cdim; ++i) assert(diffusion_gen->surf[i]);

  diffusion_gen->auxfields.Dij = 0;
  diffusion_gen->conf_range = *conf_range;

  diffusion_gen->eqn.flags = 0;
  diffusion_gen->eqn.ref_count = gkyl_ref_count_init(gkyl_diffusion_gen_free);
  diffusion_gen->eqn.on_dev = &diffusion_gen->eqn;
  
  return &diffusion_gen->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_diffusion_gen_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range)
{
  assert(false);
  return 0;
}

#endif
