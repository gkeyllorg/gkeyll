#include "gkyl_dg_eqn.h"
#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_vlasov_pkpm.h>
#include <gkyl_dg_vlasov_pkpm_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_vlasov_pkpm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  
  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_vlasov_pkpm *vlasov_pkpm = container_of(base->on_dev, struct dg_vlasov_pkpm, eqn);
    gkyl_cu_free(vlasov_pkpm);
  }
  
  struct dg_vlasov_pkpm *vlasov_pkpm = container_of(base, struct dg_vlasov_pkpm, eqn);
  gkyl_free(vlasov_pkpm);
}

void
gkyl_vlasov_pkpm_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_pkpm_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    gkyl_vlasov_pkpm_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_vlasov_pkpm *vlasov_pkpm = container_of(eqn, struct dg_vlasov_pkpm, eqn);
  vlasov_pkpm->auxfields.bvar = auxin.bvar;
  vlasov_pkpm->auxfields.bvar_surf = auxin.bvar_surf;
  vlasov_pkpm->auxfields.pkpm_prim = auxin.pkpm_prim;
  vlasov_pkpm->auxfields.pkpm_prim_surf = auxin.pkpm_prim_surf;
  vlasov_pkpm->auxfields.max_b = auxin.max_b;
  vlasov_pkpm->auxfields.pkpm_lax = auxin.pkpm_lax;
  vlasov_pkpm->auxfields.div_b = auxin.div_b;  
  vlasov_pkpm->auxfields.pkpm_accel_vars = auxin.pkpm_accel_vars;
  vlasov_pkpm->auxfields.g_dist_source = auxin.g_dist_source;
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_pkpm_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const struct gkyl_range* phase_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_vlasov_pkpm_cu_dev_new(cbasis, pbasis, conf_range, phase_range);
  } 
#endif
  struct dg_vlasov_pkpm *vlasov_pkpm = gkyl_malloc(sizeof(struct dg_vlasov_pkpm));

  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;

  vlasov_pkpm->cdim = cdim;
  vlasov_pkpm->pdim = pdim;

  vlasov_pkpm->eqn.num_equations = 2;
  vlasov_pkpm->eqn.surf_term = surf;
  vlasov_pkpm->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_vlasov_pkpm_vol_kern_list *vol_kernels;
  const gkyl_dg_vlasov_pkpm_stream_surf_kern_list *stream_surf_x_kernels, *stream_surf_y_kernels, *stream_surf_z_kernels;
  const gkyl_dg_vlasov_pkpm_accel_surf_kern_list *accel_surf_vpar_kernels;
  const gkyl_dg_vlasov_pkpm_accel_boundary_surf_kern_list *accel_boundary_surf_vpar_kernels;
  
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      stream_surf_x_kernels = ser_stream_surf_x_kernels;
      stream_surf_y_kernels = ser_stream_surf_y_kernels;
      stream_surf_z_kernels = ser_stream_surf_z_kernels;
      accel_surf_vpar_kernels = ser_accel_surf_vpar_kernels;
      accel_boundary_surf_vpar_kernels = ser_accel_boundary_surf_vpar_kernels;
      
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = ten_vol_kernels;
      stream_surf_x_kernels = ten_stream_surf_x_kernels;
      stream_surf_y_kernels = ten_stream_surf_y_kernels;
      stream_surf_z_kernels = ten_stream_surf_z_kernels;
      accel_surf_vpar_kernels = ten_accel_surf_vpar_kernels;
      accel_boundary_surf_vpar_kernels = ten_accel_boundary_surf_vpar_kernels;
      
      break;      
      
    default:
      assert(false);
      break;    
  }  

  vlasov_pkpm->eqn.vol_term = CK(vol_kernels,cdim,poly_order);

  vlasov_pkpm->stream_surf[0] = CK(stream_surf_x_kernels,cdim,poly_order);
  if (cdim>1)
    vlasov_pkpm->stream_surf[1] = CK(stream_surf_y_kernels,cdim,poly_order);
  if (cdim>2)
    vlasov_pkpm->stream_surf[2] = CK(stream_surf_z_kernels,cdim,poly_order);

  vlasov_pkpm->accel_surf = CK(accel_surf_vpar_kernels,cdim,poly_order);

  vlasov_pkpm->accel_boundary_surf = CK(accel_boundary_surf_vpar_kernels,cdim,poly_order);

  // ensure non-NULL pointers
  for (int i=0; i<cdim; ++i) assert(vlasov_pkpm->stream_surf[i]);
  assert(vlasov_pkpm->accel_surf);
  assert(vlasov_pkpm->accel_boundary_surf);

  vlasov_pkpm->auxfields.bvar = 0;  
  vlasov_pkpm->auxfields.bvar_surf = 0;  
  vlasov_pkpm->auxfields.pkpm_prim = 0;
  vlasov_pkpm->auxfields.pkpm_prim_surf = 0;
  vlasov_pkpm->auxfields.max_b = 0;
  vlasov_pkpm->auxfields.pkpm_lax = 0;
  vlasov_pkpm->auxfields.div_b = 0;
  vlasov_pkpm->auxfields.pkpm_accel_vars = 0;
  vlasov_pkpm->auxfields.g_dist_source = 0;
  vlasov_pkpm->conf_range = *conf_range;
  vlasov_pkpm->phase_range = *phase_range;

  vlasov_pkpm->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(vlasov_pkpm->eqn.flags);

  vlasov_pkpm->eqn.ref_count = gkyl_ref_count_init(gkyl_vlasov_pkpm_free);
  vlasov_pkpm->eqn.on_dev = &vlasov_pkpm->eqn; // CPU eqn obj points to itself
  
  return &vlasov_pkpm->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_vlasov_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* phase_range)
{
  assert(false);
  return 0;
}

#endif
