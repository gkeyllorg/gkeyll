/* -*- c++ -*- */

extern "C" {
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_canonical_pb.h>
#include <gkyl_dg_canonical_pb_priv.h>
#include <gkyl_util.h>
#include "gkyl_dg_eqn.h"
}

#include <cassert>

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_canonical_pb_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, 
  const struct gkyl_array *hamil, const struct gkyl_array *alpha_surf, 
  const struct gkyl_array *sgn_alpha_surf, const struct gkyl_array *const_sgn_alpha)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  canonical_pb->auxfields.hamil = hamil;
  canonical_pb->auxfields.alpha_surf = alpha_surf;
  canonical_pb->auxfields.sgn_alpha_surf = sgn_alpha_surf;
  canonical_pb->auxfields.const_sgn_alpha = const_sgn_alpha;
}
// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_canonical_pb_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_canonical_pb_auxfields auxin)
{
  gkyl_canonical_pb_set_auxfields_cu_kernel<<<1,1>>>(eqn, 
    auxin.hamil->on_dev, auxin.alpha_surf->on_dev, 
    auxin.sgn_alpha_surf->on_dev, auxin.const_sgn_alpha->on_dev);
}

// CUDA kernel to set device pointers to range object and canonical_pb kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_canonical_pb_set_cu_dev_ptrs(struct dg_canonical_pb *canonical_pb, enum gkyl_basis_type b_type, int cv_index,
  int cdim, int vdim, int poly_order)
{

  canonical_pb->auxfields.hamil = 0;  
  canonical_pb->auxfields.alpha_surf = 0;
  canonical_pb->auxfields.sgn_alpha_surf = 0;
  canonical_pb->auxfields.const_sgn_alpha = 0;

  canonical_pb->eqn.surf_term = surf;
  canonical_pb->eqn.boundary_surf_term = boundary_surf;
  canonical_pb->eqn.boundary_diag_term = boundary_diag;

  const gkyl_dg_canonical_pb_vol_kern_list *vol_kernels;
  const gkyl_dg_canonical_pb_stream_surf_kern_list *stream_surf_x_kernels, *stream_surf_y_kernels, *stream_surf_z_kernels;
  const gkyl_dg_canonical_pb_accel_surf_kern_list *accel_surf_vx_kernels, *accel_surf_vy_kernels, *accel_surf_vz_kernels;
  const gkyl_dg_canonical_pb_stream_boundary_surf_kern_list *stream_boundary_surf_x_kernels, *stream_boundary_surf_y_kernels,
    *stream_boundary_surf_z_kernels;
  const gkyl_dg_canonical_pb_accel_boundary_surf_kern_list *accel_boundary_surf_vx_kernels, *accel_boundary_surf_vy_kernels,
    *accel_boundary_surf_vz_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Verify that the poly-order is 2 for ser case
      assert(poly_order == 2);
      vol_kernels = ser_vol_kernels;
      stream_surf_x_kernels = ser_stream_surf_x_kernels;
      stream_surf_y_kernels = ser_stream_surf_y_kernels;
      stream_surf_z_kernels = ser_stream_surf_z_kernels;
      accel_surf_vx_kernels = ser_accel_surf_vx_kernels;
      accel_surf_vy_kernels = ser_accel_surf_vy_kernels;
      accel_surf_vz_kernels = ser_accel_surf_vz_kernels;
      stream_boundary_surf_x_kernels = ser_stream_boundary_surf_x_kernels;
      stream_boundary_surf_y_kernels = ser_stream_boundary_surf_y_kernels;
      stream_boundary_surf_z_kernels = ser_stream_boundary_surf_z_kernels;
      accel_boundary_surf_vx_kernels = ser_accel_boundary_surf_vx_kernels;
      accel_boundary_surf_vy_kernels = ser_accel_boundary_surf_vy_kernels;
      accel_boundary_surf_vz_kernels = ser_accel_boundary_surf_vz_kernels;
      break;

    case GKYL_BASIS_MODAL_HYBRID:
      // Verify that the poly-order is 1 for hybrid case
      assert(poly_order == 1);
      vol_kernels = ser_vol_kernels;
      stream_surf_x_kernels = ser_stream_surf_x_kernels;
      stream_surf_y_kernels = ser_stream_surf_y_kernels;
      stream_surf_z_kernels = ser_stream_surf_z_kernels;
      accel_surf_vx_kernels = ser_accel_surf_vx_kernels;
      accel_surf_vy_kernels = ser_accel_surf_vy_kernels;
      accel_surf_vz_kernels = ser_accel_surf_vz_kernels;
      stream_boundary_surf_x_kernels = ser_stream_boundary_surf_x_kernels;
      stream_boundary_surf_y_kernels = ser_stream_boundary_surf_y_kernels;
      stream_boundary_surf_z_kernels = ser_stream_boundary_surf_z_kernels;
      accel_boundary_surf_vx_kernels = ser_accel_boundary_surf_vx_kernels;
      accel_boundary_surf_vy_kernels = ser_accel_boundary_surf_vy_kernels;
      accel_boundary_surf_vz_kernels = ser_accel_boundary_surf_vz_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = tensor_vol_kernels;
      stream_surf_x_kernels = tensor_stream_surf_x_kernels;
      stream_surf_y_kernels = tensor_stream_surf_y_kernels;
      stream_surf_z_kernels = tensor_stream_surf_z_kernels;
      accel_surf_vx_kernels = tensor_accel_surf_vx_kernels;
      accel_surf_vy_kernels = tensor_accel_surf_vy_kernels;
      accel_surf_vz_kernels = tensor_accel_surf_vz_kernels;
      stream_boundary_surf_x_kernels = tensor_stream_boundary_surf_x_kernels;
      stream_boundary_surf_y_kernels = tensor_stream_boundary_surf_y_kernels;
      stream_boundary_surf_z_kernels = tensor_stream_boundary_surf_z_kernels;
      accel_boundary_surf_vx_kernels = tensor_accel_boundary_surf_vx_kernels;
      accel_boundary_surf_vy_kernels = tensor_accel_boundary_surf_vy_kernels;
      accel_boundary_surf_vz_kernels = tensor_accel_boundary_surf_vz_kernels;
      break;

    default:
      assert(false);
      break;    
  }  

  canonical_pb->eqn.vol_term = CK(vol_kernels,cv_index,poly_order);

  canonical_pb->stream_surf[0] = CK(stream_surf_x_kernels,cv_index,poly_order);
  if (cdim>1)
    canonical_pb->stream_surf[1] = CK(stream_surf_y_kernels,cv_index,poly_order);
  if (cdim>2)
    canonical_pb->stream_surf[2] = CK(stream_surf_z_kernels,cv_index,poly_order);

  canonical_pb->accel_surf[0] = CK(accel_surf_vx_kernels,cv_index,poly_order);
  if (vdim>1)
    canonical_pb->accel_surf[1] = CK(accel_surf_vy_kernels,cv_index,poly_order);
  if (vdim>2)
    canonical_pb->accel_surf[2] = CK(accel_surf_vz_kernels,cv_index,poly_order);

  canonical_pb->stream_boundary_surf[0] = CK(stream_boundary_surf_x_kernels,cv_index,poly_order);
  if (cdim>1)
    canonical_pb->stream_boundary_surf[1] = CK(stream_boundary_surf_y_kernels,cv_index,poly_order);
  if (cdim>2)
    canonical_pb->stream_boundary_surf[2] = CK(stream_boundary_surf_z_kernels,cv_index,poly_order);

  canonical_pb->accel_boundary_surf[0] = CK(accel_boundary_surf_vx_kernels,cv_index,poly_order);
  if (vdim>1)
    canonical_pb->accel_boundary_surf[1] = CK(accel_boundary_surf_vy_kernels,cv_index,poly_order);
  if (vdim>2)
    canonical_pb->accel_boundary_surf[2] = CK(accel_boundary_surf_vz_kernels,cv_index,poly_order);
}



struct gkyl_dg_eqn*
gkyl_dg_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* phase_range)
{
  struct dg_canonical_pb *canonical_pb = (struct dg_canonical_pb*)  gkyl_malloc(sizeof(struct dg_canonical_pb));


  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  canonical_pb->cdim = cdim;
  canonical_pb->pdim = pdim;

  canonical_pb->eqn.num_equations = 1;
  canonical_pb->phase_range = *phase_range;

  canonical_pb->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(canonical_pb->eqn.flags);
  canonical_pb->eqn.ref_count = gkyl_ref_count_init(gkyl_canonical_pb_free);

  // copy the host struct to device struct
  struct dg_canonical_pb *canonical_pb_cu = (struct dg_canonical_pb*) gkyl_cu_malloc(sizeof(struct dg_canonical_pb));
  gkyl_cu_memcpy(canonical_pb_cu, canonical_pb, sizeof(struct dg_canonical_pb), GKYL_CU_MEMCPY_H2D);

  dg_canonical_pb_set_cu_dev_ptrs<<<1,1>>>(canonical_pb_cu, pbasis->b_type, cv_index[cdim].vdim[vdim], cdim, vdim, poly_order);

  // set parent on_dev pointer
  canonical_pb->eqn.on_dev = &canonical_pb_cu->eqn;

  return &canonical_pb->eqn;
}
