#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_maxwell_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order] 

void 
gkyl_maxwell_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_maxwell *maxwell = container_of(base->on_dev, struct dg_maxwell, eqn);
    gkyl_cu_free(maxwell);
  }  
  
  struct dg_maxwell *maxwell = container_of(base, struct dg_maxwell, eqn);

  // Free the conf_flux_surf acquisition
  if (maxwell->use_conf_flux_surf) {
    gkyl_array_release(maxwell->conf_flux_surf);
    gkyl_surf_and_vol_node_arrays_release(maxwell->lapse);
    gkyl_surf_and_vol_node_arrays_release(maxwell->shift);
    gkyl_surf_and_vol_node_arrays_release(maxwell->h_ij);
    gkyl_surf_and_vol_node_arrays_release(maxwell->h_ij_inv);
    gkyl_surf_and_vol_node_arrays_release(maxwell->det_h);
  }

  gkyl_free(maxwell);
}

struct gkyl_dg_eqn*
gkyl_dg_maxwell_inew(const struct gkyl_dg_maxwell_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if(inp->use_gpu) {
    return gkyl_dg_maxwell_cu_dev_inew(inp);
  } 
#endif
  struct dg_maxwell *maxwell = gkyl_malloc(sizeof(struct dg_maxwell));

  int cdim = inp->cbasis->ndim;
  int poly_order = inp->cbasis->poly_order;

  // Configuration-space range + position-map Jacobian for the curl kernels.
  // The Jacobian (borrowed) is the per-conf-cell J of the C^0 linear position
  // map. It is required: callers without a nonuniform position map pass an
  // identity position map's Jacobian (J = 1 at every quadrature point), which
  // keeps the kernels bit-identical to the uniform-grid case.
  if (inp->crange) maxwell->crange = *inp->crange;
  assert(inp->jacob_pos);
  maxwell->jacob_pos = inp->jacob_pos;

  // For configuration fluxes in dg gr maxwell
  maxwell->use_conf_flux_surf = false;
  if (inp->field_id == GKYL_FIELD_GR_D_B){
    assert(inp->crange);
    maxwell->use_conf_flux_surf = true;
    maxwell->crange = *inp->crange;
    maxwell->conf_flux_surf = gkyl_array_acquire(inp->conf_flux_surf);
    maxwell->lapse = gkyl_surf_and_vol_node_arrays_acquire(inp->lapse);
    maxwell->shift = gkyl_surf_and_vol_node_arrays_acquire(inp->shift);
    maxwell->h_ij = gkyl_surf_and_vol_node_arrays_acquire(inp->h_ij);
    maxwell->h_ij_inv = gkyl_surf_and_vol_node_arrays_acquire(inp->h_ij_inv);
    maxwell->det_h = gkyl_surf_and_vol_node_arrays_acquire(inp->det_h);
  }

  const gkyl_dg_maxwell_vol_kern_list *vol_kernels;
  const gkyl_dg_maxwell_gr_maxwell_vol_kern_list *vol_gr_kernels;
  const gkyl_dg_maxwell_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;
  const gkyl_dg_maxwell_surf_from_flux_kern_list *surf_from_flux_x_kernels, *surf_from_flux_y_kernels, *surf_from_flux_z_kernels;

  switch (inp->cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if ( maxwell->use_conf_flux_surf ) {
        vol_gr_kernels = ser_gr_maxwell_vol_kernels;
        surf_from_flux_x_kernels = ser_gr_maxwell_surf_from_flux_x_kernels;
        surf_from_flux_y_kernels = ser_gr_maxwell_surf_from_flux_y_kernels;
        surf_from_flux_z_kernels = ser_gr_maxwell_surf_from_flux_z_kernels;
      }
      else {
        vol_kernels = ser_vol_kernels;
        surf_x_kernels = ser_surf_x_kernels;
        surf_y_kernels = ser_surf_y_kernels;
        surf_z_kernels = ser_surf_z_kernels;
      }

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      if ( maxwell->use_conf_flux_surf ) {
        vol_gr_kernels = ten_gr_maxwell_vol_kernels;
        surf_from_flux_x_kernels = tensor_gr_maxwell_surf_from_flux_x_kernels;
        surf_from_flux_y_kernels = tensor_gr_maxwell_surf_from_flux_y_kernels;
        surf_from_flux_z_kernels = tensor_gr_maxwell_surf_from_flux_z_kernels;
      }
      else {
        vol_kernels = ten_vol_kernels;
        surf_x_kernels = ten_surf_x_kernels;
        surf_y_kernels = ten_surf_y_kernels;
        surf_z_kernels = ten_surf_z_kernels;
      }
      
      break;

    default:
      assert(false);
      break;    
  }  
    
  maxwell->eqn.num_equations = 8;
  maxwell->eqn.surf_term = surf;
  maxwell->eqn.boundary_surf_term = boundary_surf;

  if (inp->field_id == GKYL_FIELD_GR_D_B) assert(fabs(inp->lightSpeed - 1.0) < 1e-12);
  maxwell->maxwell_data.c = inp->lightSpeed;
  maxwell->maxwell_data.chi = inp->lightSpeed*inp->elcErrorSpeedFactor;
  maxwell->maxwell_data.gamma = inp->lightSpeed*inp->mgnErrorSpeedFactor;
  maxwell->gr_maxwell_data.chi = inp->elcErrorSpeedFactor;
  maxwell->gr_maxwell_data.gamma = inp->mgnErrorSpeedFactor;
  maxwell->gr_maxwell_data.K_phi = 0.0;
  maxwell->gr_maxwell_data.K_psi = 0.0;

  // For volume kernel selection
  if ( maxwell->use_conf_flux_surf ) {
    maxwell->eqn.vol_term = vol;
    maxwell->vol = CK(vol_gr_kernels, cdim, poly_order);
  }
  else {
    maxwell->eqn.vol_term = CK(vol_kernels, cdim, poly_order);
  }

  // For surf-from-flux kernels
  if ( maxwell->use_conf_flux_surf ) {
    maxwell->surf_from_flux[0] = CK(surf_from_flux_x_kernels, cdim, poly_order);
    if (cdim>1)
      maxwell->surf_from_flux[1] = CK(surf_from_flux_y_kernels, cdim, poly_order);
    if (cdim>2)
      maxwell->surf_from_flux[2] = CK(surf_from_flux_z_kernels, cdim, poly_order);
  } 
  else {
    maxwell->surf[0] = CK(surf_x_kernels, cdim, poly_order);
    if (cdim>1)
      maxwell->surf[1] = CK(surf_y_kernels, cdim, poly_order);
    if (cdim>2)
      maxwell->surf[2] = CK(surf_z_kernels, cdim, poly_order);
  }

  // ensure non-NULL pointers 
  if ( maxwell->use_conf_flux_surf ) {
    for (int i=0; i<cdim; ++i) assert(maxwell->surf_from_flux[i]);
  }
  else {
    for (int i=0; i<cdim; ++i) assert(maxwell->surf[i]);
  }

  maxwell->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(maxwell->eqn.flags);
  maxwell->eqn.ref_count = gkyl_ref_count_init(gkyl_maxwell_free);
  maxwell->eqn.on_dev = &maxwell->eqn; // CPU eqn obj points to itself
  
  return &maxwell->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_maxwell_cu_dev_inew(const struct gkyl_dg_maxwell_inp *inp)
{
  assert(false);
  return 0;
}

#endif
