#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_advection.h>
#include <gkyl_dg_advection_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void 
gkyl_advection_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_advection *advection = container_of(base->on_dev, struct dg_advection, eqn);
    gkyl_cu_free(advection);
  }  
  
  struct dg_advection *advection = container_of(base, struct dg_advection, eqn);
  gkyl_free(advection);
}

void
gkyl_advection_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_advection_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.u_i)) {
    gkyl_advection_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  advection->auxfields.u_i = auxin.u_i;
}

struct gkyl_dg_eqn*
gkyl_dg_advection_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_advection_cu_dev_new(cbasis, conf_range);
  } 
#endif
  struct dg_advection *advection = gkyl_malloc(sizeof(struct dg_advection));

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

  const gkyl_dg_advection_vol_kern_list *vol_kernels;
  const gkyl_dg_advection_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_x_kernels = ser_surf_x_kernels;
      surf_y_kernels = ser_surf_y_kernels;
      surf_z_kernels = ser_surf_z_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
    
  advection->eqn.num_equations = 1;
  advection->eqn.surf_term = surf;
  advection->eqn.boundary_surf_term = boundary_surf;

  advection->eqn.vol_term =  CK(vol_kernels, cdim, poly_order);

  advection->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    advection->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    advection->surf[2] = CK(surf_z_kernels, cdim, poly_order);

  // ensure non-NULL pointers 
  for (int i=0; i<cdim; ++i) assert(advection->surf[i]);

  advection->auxfields.u_i = 0;  
  advection->conf_range = *conf_range;

  advection->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(advection->eqn.flags);
  advection->eqn.ref_count = gkyl_ref_count_init(gkyl_advection_free);
  advection->eqn.on_dev = &advection->eqn; // CPU eqn obj points to itself
  
  return &advection->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_advection_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range)
{
  assert(false);
  return 0;
}

#endif
