#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_canonical_pb.h>
#include <gkyl_dg_canonical_pb_priv.h>
#include <gkyl_util.h>
#include "gkyl_dg_eqn.h"

void
gkyl_canonical_pb_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  
  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_canonical_pb *canonical_pb = container_of(base->on_dev, struct dg_canonical_pb, eqn);
    gkyl_cu_free(canonical_pb);
  }
  
  struct dg_canonical_pb *canonical_pb = container_of(base, struct dg_canonical_pb, eqn);
  gkyl_free(canonical_pb);
}

void
gkyl_canonical_pb_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_canonical_pb_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    gkyl_canonical_pb_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  canonical_pb->auxfields.hamil = auxin.hamil;
  canonical_pb->auxfields.alpha_surf = auxin.alpha_surf;
  canonical_pb->auxfields.sgn_alpha_surf = auxin.sgn_alpha_surf;
  canonical_pb->auxfields.const_sgn_alpha = auxin.const_sgn_alpha;
}

struct gkyl_dg_eqn*
gkyl_dg_canonical_pb_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* phase_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_canonical_pb_cu_dev_new(cbasis, pbasis, phase_range);
  } 
#endif
  struct dg_canonical_pb *canonical_pb = gkyl_malloc(sizeof(struct dg_canonical_pb));


  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  canonical_pb->cdim = cdim;
  canonical_pb->pdim = pdim;

  canonical_pb->eqn.num_equations = 1;
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
  
  switch (pbasis->b_type) {
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
  int cv_index_val = cv_index[cdim].vdim[vdim];
  canonical_pb->eqn.vol_term = CK(vol_kernels,cv_index_val,poly_order);

  canonical_pb->stream_surf[0] = CK(stream_surf_x_kernels,cv_index_val,poly_order);
  if (cdim>1)
    canonical_pb->stream_surf[1] = CK(stream_surf_y_kernels,cv_index_val,poly_order);
  if (cdim>2)
    canonical_pb->stream_surf[2] = CK(stream_surf_z_kernels,cv_index_val,poly_order);

  canonical_pb->accel_surf[0] = CK(accel_surf_vx_kernels,cv_index_val,poly_order);
  if (vdim>1)
    canonical_pb->accel_surf[1] = CK(accel_surf_vy_kernels,cv_index_val,poly_order);
  if (vdim>2)
    canonical_pb->accel_surf[2] = CK(accel_surf_vz_kernels,cv_index_val,poly_order);

  canonical_pb->stream_boundary_surf[0] = CK(stream_boundary_surf_x_kernels,cv_index_val,poly_order);
  if (cdim>1)
    canonical_pb->stream_boundary_surf[1] = CK(stream_boundary_surf_y_kernels,cv_index_val,poly_order);
  if (cdim>2)
    canonical_pb->stream_boundary_surf[2] = CK(stream_boundary_surf_z_kernels,cv_index_val,poly_order);

  canonical_pb->accel_boundary_surf[0] = CK(accel_boundary_surf_vx_kernels,cv_index_val,poly_order);
  if (vdim>1)
    canonical_pb->accel_boundary_surf[1] = CK(accel_boundary_surf_vy_kernels,cv_index_val,poly_order);
  if (vdim>2)
    canonical_pb->accel_boundary_surf[2] = CK(accel_boundary_surf_vz_kernels,cv_index_val,poly_order);

  // ensure non-NULL pointers
  for (int i=0; i<cdim; ++i) assert(canonical_pb->stream_surf[i]);
  for (int i=0; i<vdim; ++i) assert(canonical_pb->accel_surf[i]);
  for (int i=0; i<vdim; ++i) assert(canonical_pb->accel_boundary_surf[i]);
  for (int i=0; i<cdim; ++i) assert(canonical_pb->stream_boundary_surf[i]);

  canonical_pb->auxfields.hamil = 0;  
  canonical_pb->auxfields.alpha_surf = 0;
  canonical_pb->auxfields.sgn_alpha_surf = 0;
  canonical_pb->auxfields.const_sgn_alpha = 0;
  canonical_pb->phase_range = *phase_range;

  canonical_pb->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(canonical_pb->eqn.flags);

  canonical_pb->eqn.ref_count = gkyl_ref_count_init(gkyl_canonical_pb_free);
  canonical_pb->eqn.on_dev = &canonical_pb->eqn; // CPU eqn obj points to itself
  
  return &canonical_pb->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
 const struct gkyl_range* phase_range)
{
  assert(false);
  return 0;
}

#endif
