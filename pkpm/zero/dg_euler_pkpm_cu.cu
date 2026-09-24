/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_euler_pkpm.h>    
#include <gkyl_dg_euler_pkpm_priv.h>
#include <gkyl_util.h>
}

#include <cassert>

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_euler_pkpm_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, 
  const struct gkyl_array *vlasov_pkpm_moms, 
  const struct gkyl_array *pkpm_prim, const struct gkyl_array *pkpm_prim_surf, 
  const struct gkyl_array *pkpm_p_ij, const struct gkyl_array *pkpm_lax, 
  const struct gkyl_array *pkpm_penalization)
{
  struct dg_euler_pkpm *euler_pkpm = container_of(eqn, struct dg_euler_pkpm, eqn);
  euler_pkpm->auxfields.vlasov_pkpm_moms = vlasov_pkpm_moms;
  euler_pkpm->auxfields.pkpm_prim = pkpm_prim;
  euler_pkpm->auxfields.pkpm_prim_surf = pkpm_prim_surf;
  euler_pkpm->auxfields.pkpm_p_ij = pkpm_p_ij;
  euler_pkpm->auxfields.pkpm_lax = pkpm_lax;
  euler_pkpm->auxfields.pkpm_penalization = pkpm_penalization;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_euler_pkpm_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_pkpm_auxfields auxin)
{
  gkyl_euler_pkpm_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.vlasov_pkpm_moms->on_dev, 
    auxin.pkpm_prim->on_dev, auxin.pkpm_prim_surf->on_dev, 
    auxin.pkpm_p_ij->on_dev, auxin.pkpm_lax->on_dev, auxin.pkpm_penalization->on_dev);
}

__global__ void static
dg_euler_pkpm_set_cu_dev_ptrs(struct dg_euler_pkpm* euler_pkpm, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  euler_pkpm->auxfields.vlasov_pkpm_moms = 0; 
  euler_pkpm->auxfields.pkpm_prim = 0;
  euler_pkpm->auxfields.pkpm_prim_surf = 0;  
  euler_pkpm->auxfields.pkpm_p_ij = 0;
  euler_pkpm->auxfields.pkpm_lax = 0; 
  euler_pkpm->auxfields.pkpm_penalization = 0; 
  
  const gkyl_dg_euler_pkpm_vol_kern_list *vol_kernels;
  const gkyl_dg_euler_pkpm_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;  
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_x_kernels = ser_surf_x_kernels;
      surf_y_kernels = ser_surf_y_kernels;
      surf_z_kernels = ser_surf_z_kernels;

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = ten_vol_kernels;
      surf_x_kernels = ten_surf_x_kernels;
      surf_y_kernels = ten_surf_y_kernels;
      surf_z_kernels = ten_surf_z_kernels;
      
      break;

    default:
      assert(false);
      break;    
  }  
  
  euler_pkpm->eqn.surf_term = surf;
  euler_pkpm->eqn.boundary_surf_term = boundary_surf;

  euler_pkpm->eqn.vol_term =  CK(vol_kernels, cdim, poly_order);

  euler_pkpm->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    euler_pkpm->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    euler_pkpm->surf[2] = CK(surf_z_kernels, cdim, poly_order);
}

struct gkyl_dg_eqn*
gkyl_dg_euler_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* conf_range)
{
  struct dg_euler_pkpm *euler_pkpm = (struct dg_euler_pkpm*) gkyl_malloc(sizeof(struct dg_euler_pkpm));

  // set basic parameters
  euler_pkpm->eqn.num_equations = 3;

  euler_pkpm->conf_range = *conf_range;

  euler_pkpm->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(euler_pkpm->eqn.flags);
  euler_pkpm->eqn.ref_count = gkyl_ref_count_init(gkyl_euler_pkpm_free);

  // copy the host struct to device struct
  struct dg_euler_pkpm *euler_pkpm_cu = (struct dg_euler_pkpm*) gkyl_cu_malloc(sizeof(struct dg_euler_pkpm));
  gkyl_cu_memcpy(euler_pkpm_cu, euler_pkpm, sizeof(struct dg_euler_pkpm), GKYL_CU_MEMCPY_H2D);
  dg_euler_pkpm_set_cu_dev_ptrs<<<1,1>>>(euler_pkpm_cu, cbasis->b_type, cbasis->ndim, cbasis->poly_order);

  // set parent on_dev pointer
  euler_pkpm->eqn.on_dev = &euler_pkpm_cu->eqn;

  return &euler_pkpm->eqn;
}
