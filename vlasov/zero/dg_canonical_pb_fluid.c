#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_canonical_pb_fluid.h>
#include <gkyl_dg_canonical_pb_fluid_priv.h>
#include <gkyl_util.h>
#include "gkyl_dg_eqn.h"

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_canonical_pb_fluid_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  
  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_canonical_pb_fluid *can_pb_fluid = container_of(base->on_dev, struct dg_canonical_pb_fluid, eqn);
    gkyl_cu_free(can_pb_fluid);
  }
  
  struct dg_canonical_pb_fluid *can_pb_fluid = container_of(base, struct dg_canonical_pb_fluid, eqn);
  gkyl_free(can_pb_fluid);
}

void
gkyl_canonical_pb_fluid_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_canonical_pb_fluid_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    gkyl_canonical_pb_fluid_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_canonical_pb_fluid *can_pb_fluid = container_of(eqn, struct dg_canonical_pb_fluid, eqn);
  can_pb_fluid->auxfields.phi = auxin.phi;
  can_pb_fluid->auxfields.alpha_surf = auxin.alpha_surf;
  can_pb_fluid->auxfields.sgn_alpha_surf = auxin.sgn_alpha_surf;
  can_pb_fluid->auxfields.const_sgn_alpha = auxin.const_sgn_alpha;
}

struct gkyl_dg_eqn*
gkyl_dg_canonical_pb_fluid_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_wv_eqn *wv_eqn, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_canonical_pb_fluid_cu_dev_new(cbasis, conf_range, wv_eqn);
  } 
#endif
  struct dg_canonical_pb_fluid *can_pb_fluid = gkyl_malloc(sizeof(struct dg_canonical_pb_fluid));

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  assert(cdim == 2); // Only defined for cdim = 2 right now

  can_pb_fluid->cdim = cdim;

  can_pb_fluid->eqn.num_equations = wv_eqn->num_equations;
  can_pb_fluid->eqn.surf_term = surf;

  const gkyl_dg_canonical_pb_fluid_vol_kern_list *vol_kernels;
  const gkyl_dg_canonical_pb_fluid_surf_kern_list *surf_x_kernels, *surf_y_kernels;
  
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (can_pb_fluid->eqn.num_equations == 2) {
        vol_kernels = ser_two_fluid_vol_kernels;
        surf_x_kernels = ser_two_fluid_surf_x_kernels;
        surf_y_kernels = ser_two_fluid_surf_y_kernels;        
      } 
      else {
        vol_kernels = ser_vol_kernels;
        surf_x_kernels = ser_surf_x_kernels;
        surf_y_kernels = ser_surf_y_kernels;
      } 
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      if (can_pb_fluid->eqn.num_equations == 2) {
        vol_kernels = tensor_two_fluid_vol_kernels;
        surf_x_kernels = tensor_two_fluid_surf_x_kernels;
        surf_y_kernels = tensor_two_fluid_surf_y_kernels;        
      } 
      else {
        vol_kernels = tensor_vol_kernels;
        surf_x_kernels = tensor_surf_x_kernels;
        surf_y_kernels = tensor_surf_y_kernels;
      } 
      break;

    default:
      assert(false);
      break;    
  }  
  can_pb_fluid->eqn.vol_term = CK(vol_kernels,cdim,poly_order);

  can_pb_fluid->surf[0] = CK(surf_x_kernels,cdim,poly_order);
  can_pb_fluid->surf[1] = CK(surf_y_kernels,cdim,poly_order);

  // ensure non-NULL pointers
  for (int i=0; i<cdim; ++i) assert(can_pb_fluid->surf[i]);

  can_pb_fluid->auxfields.phi = 0;  
  can_pb_fluid->auxfields.alpha_surf = 0;
  can_pb_fluid->auxfields.sgn_alpha_surf = 0;
  can_pb_fluid->auxfields.const_sgn_alpha = 0;
  can_pb_fluid->conf_range = *conf_range;

  can_pb_fluid->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(can_pb_fluid->eqn.flags);

  can_pb_fluid->eqn.ref_count = gkyl_ref_count_init(gkyl_canonical_pb_fluid_free);
  can_pb_fluid->eqn.on_dev = &can_pb_fluid->eqn; // CPU eqn obj points to itself
  
  return &can_pb_fluid->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_canonical_pb_fluid_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_wv_eqn *wv_eqn)
{
  assert(false);
  return 0;
}

#endif
