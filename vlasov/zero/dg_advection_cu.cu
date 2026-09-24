/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_advection.h>    
#include <gkyl_dg_advection_priv.h>
}

#include <cassert>

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_advection_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *u_i)
{
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  advection->auxfields.u_i = u_i;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_advection_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_advection_auxfields auxin)
{
  gkyl_advection_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.u_i->on_dev);
}

__global__ void static
dg_advection_set_cu_dev_ptrs(struct dg_advection* advection, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  advection->auxfields.u_i = 0; 

  const gkyl_dg_advection_vol_kern_list *vol_kernels;
  const gkyl_dg_advection_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;  
  
  switch (b_type) {
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
  
  advection->eqn.surf_term = surf;
  advection->eqn.boundary_surf_term = boundary_surf;

  advection->eqn.vol_term =  CK(vol_kernels, cdim, poly_order);

  advection->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    advection->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    advection->surf[2] = CK(surf_z_kernels, cdim, poly_order);
}

struct gkyl_dg_eqn*
gkyl_dg_advection_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range)
{
  struct dg_advection *advection = (struct dg_advection*) gkyl_malloc(sizeof(struct dg_advection));

  // set basic parameters
  advection->eqn.num_equations = 1;
  advection->conf_range = *conf_range;

  advection->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(advection->eqn.flags);
  advection->eqn.ref_count = gkyl_ref_count_init(gkyl_advection_free);

  // copy the host struct to device struct
  struct dg_advection *advection_cu = (struct dg_advection*) gkyl_cu_malloc(sizeof(struct dg_advection));
  gkyl_cu_memcpy(advection_cu, advection, sizeof(struct dg_advection), GKYL_CU_MEMCPY_H2D);
  dg_advection_set_cu_dev_ptrs<<<1,1>>>(advection_cu, cbasis->b_type, cbasis->ndim, cbasis->poly_order);

  // set parent on_dev pointer
  advection->eqn.on_dev = &advection_cu->eqn;

  return &advection->eqn;
}
