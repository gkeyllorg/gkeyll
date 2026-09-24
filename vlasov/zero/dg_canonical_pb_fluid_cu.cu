/* -*- c++ -*- */

extern "C" {
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_canonical_pb_fluid.h>
#include <gkyl_dg_canonical_pb_fluid_priv.h>
#include <gkyl_util.h>
#include "gkyl_dg_eqn.h"
}

#include <cassert>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_canonical_pb_fluid_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, 
  const struct gkyl_array *phi, const struct gkyl_array *alpha_surf, 
  const struct gkyl_array *sgn_alpha_surf, const struct gkyl_array *const_sgn_alpha)
{
  struct dg_canonical_pb_fluid *can_pb_fluid = container_of(eqn, struct dg_canonical_pb_fluid, eqn);
  can_pb_fluid->auxfields.phi = phi;
  can_pb_fluid->auxfields.alpha_surf = alpha_surf;
  can_pb_fluid->auxfields.sgn_alpha_surf = sgn_alpha_surf;
  can_pb_fluid->auxfields.const_sgn_alpha = const_sgn_alpha;
}
// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_canonical_pb_fluid_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_canonical_pb_fluid_auxfields auxin)
{
  gkyl_canonical_pb_fluid_set_auxfields_cu_kernel<<<1,1>>>(eqn, 
    auxin.phi->on_dev, auxin.alpha_surf->on_dev, 
    auxin.sgn_alpha_surf->on_dev, auxin.const_sgn_alpha->on_dev);
}

// CUDA kernel to set device pointers to range object and canonical_pb_fluid kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_canonical_pb_fluid_set_cu_dev_ptrs(struct dg_canonical_pb_fluid *can_pb_fluid, enum gkyl_basis_type b_type, 
  int cdim, int poly_order, int num_equations)
{

  can_pb_fluid->auxfields.phi = 0;  
  can_pb_fluid->auxfields.alpha_surf = 0;
  can_pb_fluid->auxfields.sgn_alpha_surf = 0;
  can_pb_fluid->auxfields.const_sgn_alpha = 0;

  can_pb_fluid->eqn.surf_term = surf;

  const gkyl_dg_canonical_pb_fluid_vol_kern_list *vol_kernels;
  const gkyl_dg_canonical_pb_fluid_surf_kern_list *surf_x_kernels, *surf_y_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (num_equations == 2) {
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
      if (num_equations == 2) {
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
}



struct gkyl_dg_eqn*
gkyl_dg_canonical_pb_fluid_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_wv_eqn *wv_eqn)
{
  struct dg_canonical_pb_fluid *can_pb_fluid = (struct dg_canonical_pb_fluid*)  gkyl_malloc(sizeof(struct dg_canonical_pb_fluid));


  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  assert(cdim == 2); // Only defined for cdim = 2 right now

  can_pb_fluid->cdim = cdim;

  int num_equations = wv_eqn->num_equations;
  can_pb_fluid->eqn.num_equations = num_equations;
  can_pb_fluid->conf_range = *conf_range;

  can_pb_fluid->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(can_pb_fluid->eqn.flags);
  can_pb_fluid->eqn.ref_count = gkyl_ref_count_init(gkyl_canonical_pb_fluid_free);

  // copy the host struct to device struct
  struct dg_canonical_pb_fluid *can_pb_fluid_cu = (struct dg_canonical_pb_fluid*) gkyl_cu_malloc(sizeof(struct dg_canonical_pb_fluid));
  gkyl_cu_memcpy(can_pb_fluid_cu, can_pb_fluid, sizeof(struct dg_canonical_pb_fluid), GKYL_CU_MEMCPY_H2D);

  dg_canonical_pb_fluid_set_cu_dev_ptrs<<<1,1>>>(can_pb_fluid_cu, cbasis->b_type, cdim, poly_order, num_equations);

  // set parent on_dev pointer
  can_pb_fluid->eqn.on_dev = &can_pb_fluid_cu->eqn;

  return &can_pb_fluid->eqn;
}
