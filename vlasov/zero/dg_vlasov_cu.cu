/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_vlasov.h>    
#include <gkyl_dg_vlasov_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_vlasov_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, 
  const struct gkyl_array *field, 
  const struct gkyl_array *alpha_surf, const struct gkyl_array *sgn_alpha_surf, const struct gkyl_array *const_sgn_alpha)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);
  vlasov->auxfields.field = field; 
  vlasov->auxfields.alpha_surf = alpha_surf;
  vlasov->auxfields.sgn_alpha_surf = sgn_alpha_surf;
  vlasov->auxfields.const_sgn_alpha = const_sgn_alpha;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_vlasov_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_auxfields auxin)
{
  gkyl_vlasov_set_auxfields_cu_kernel<<<1,1>>>(eqn,
    auxin.field ? auxin.field->on_dev : 0,
    auxin.alpha_surf ? auxin.alpha_surf->on_dev : 0,
    auxin.sgn_alpha_surf ? auxin.sgn_alpha_surf->on_dev : 0,
    auxin.const_sgn_alpha ? auxin.const_sgn_alpha->on_dev : 0);
}

// CUDA kernel to set device pointers to range object and vlasov kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_vlasov_set_cu_dev_ptrs(struct dg_vlasov *vlasov, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order, 
  enum gkyl_model_id model_id, enum gkyl_field_id field_id)
{
  vlasov->auxfields.field = 0;
  vlasov->auxfields.alpha_surf = 0;
  vlasov->auxfields.sgn_alpha_surf = 0;
  vlasov->auxfields.const_sgn_alpha = 0;

  vlasov->eqn.surf_term = surf;
  vlasov->eqn.boundary_surf_term = boundary_surf;
  vlasov->eqn.boundary_diag_term = boundary_diag;

  const gkyl_dg_vlasov_stream_vol_kern_list *stream_vol_kernels;
  const gkyl_dg_vlasov_vol_kern_list *vol_kernels;

  const gkyl_dg_vlasov_stream_surf_kern_list *stream_surf_x_kernels, 
    *stream_surf_y_kernels, 
    *stream_surf_z_kernels;
  const gkyl_dg_vlasov_accel_surf_kern_list *accel_surf_vx_kernels, 
    *accel_surf_vy_kernels, 
    *accel_surf_vz_kernels;
  const gkyl_dg_vlasov_stream_boundary_surf_kern_list *stream_boundary_surf_x_kernels, 
    *stream_boundary_surf_y_kernels,
    *stream_boundary_surf_z_kernels;
  const gkyl_dg_vlasov_accel_boundary_surf_kern_list *accel_boundary_surf_vx_kernels, 
    *accel_boundary_surf_vy_kernels,
    *accel_boundary_surf_vz_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      stream_vol_kernels = ser_stream_vol_kernels;
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
      stream_vol_kernels = tensor_stream_vol_kernels;
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
  if (field_id == GKYL_FIELD_NULL) {
    vlasov->eqn.vol_term = stream_vol_kernels[cv_index].kernels[poly_order];
  }
  else {
    vlasov->eqn.vol_term = vol_kernels[cv_index].kernels[poly_order];

    vlasov->accel_surf[0] = accel_surf_vx_kernels[cv_index].kernels[poly_order];
    if (vdim>1)
      vlasov->accel_surf[1] = accel_surf_vy_kernels[cv_index].kernels[poly_order];
    if (vdim>2)
      vlasov->accel_surf[2] = accel_surf_vz_kernels[cv_index].kernels[poly_order];

    vlasov->accel_boundary_surf[0] = accel_boundary_surf_vx_kernels[cv_index].kernels[poly_order];
    if (vdim>1)
      vlasov->accel_boundary_surf[1] = accel_boundary_surf_vy_kernels[cv_index].kernels[poly_order];
    if (vdim>2)
      vlasov->accel_boundary_surf[2] = accel_boundary_surf_vz_kernels[cv_index].kernels[poly_order];
  }
  // Streaming kernels are the same for each field_id
  vlasov->stream_surf[0] = stream_surf_x_kernels[cv_index].kernels[poly_order];
  if (cdim>1)
    vlasov->stream_surf[1] = stream_surf_y_kernels[cv_index].kernels[poly_order];
  if (cdim>2)
    vlasov->stream_surf[2] = stream_surf_z_kernels[cv_index].kernels[poly_order];

  vlasov->stream_boundary_surf[0] = stream_boundary_surf_x_kernels[cv_index].kernels[poly_order];
  if (cdim>1)
    vlasov->stream_boundary_surf[1] = stream_boundary_surf_y_kernels[cv_index].kernels[poly_order];
  if (cdim>2)
    vlasov->stream_boundary_surf[2] = stream_boundary_surf_z_kernels[cv_index].kernels[poly_order];
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const struct gkyl_range* phase_range,
  enum gkyl_model_id model_id, enum gkyl_field_id field_id)
{
  struct dg_vlasov *vlasov = (struct dg_vlasov*) gkyl_malloc(sizeof(struct dg_vlasov));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov->cdim = cdim;
  vlasov->pdim = pdim;

  vlasov->eqn.num_equations = 1;
  vlasov->conf_range = *conf_range;
  vlasov->phase_range = *phase_range;

  vlasov->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(vlasov->eqn.flags);
  vlasov->eqn.ref_count = gkyl_ref_count_init(gkyl_vlasov_free);

  // copy the host struct to device struct
  struct dg_vlasov *vlasov_cu = (struct dg_vlasov*) gkyl_cu_malloc(sizeof(struct dg_vlasov));
  gkyl_cu_memcpy(vlasov_cu, vlasov, sizeof(struct dg_vlasov), GKYL_CU_MEMCPY_H2D);

  dg_vlasov_set_cu_dev_ptrs<<<1,1>>>(vlasov_cu, cbasis->b_type, cv_index[cdim].vdim[vdim],
    cdim, vdim, poly_order, model_id, field_id);

  // set parent on_dev pointer
  vlasov->eqn.on_dev = &vlasov_cu->eqn;
  
  return &vlasov->eqn;
}
