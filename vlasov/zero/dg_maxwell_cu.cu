/* -*- c++ -*- */

#include <math.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_maxwell.h>    
#include <gkyl_dg_maxwell_priv.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>
}

#include <cassert>

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

__global__ void static
dg_maxwell_set_cu_dev_ptrs(struct dg_maxwell* maxwell, enum gkyl_basis_type b_type, int cdim, 
  bool use_conf_flux_surf, int poly_order)
{
  const gkyl_dg_maxwell_vol_kern_list *vol_kernels;
  const gkyl_dg_maxwell_gr_maxwell_vol_kern_list *vol_gr_kernels;
  const gkyl_dg_maxwell_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;
  const gkyl_dg_maxwell_surf_from_flux_kern_list *surf_from_flux_x_kernels, *surf_from_flux_y_kernels, *surf_from_flux_z_kernels;
 
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if ( use_conf_flux_surf ) {
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
      if ( use_conf_flux_surf ) {
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
  
  maxwell->eqn.surf_term = surf;
  maxwell->eqn.boundary_surf_term = boundary_surf;

  if ( use_conf_flux_surf ) {
    maxwell->eqn.vol_term = vol;  
    maxwell->vol = CK(vol_gr_kernels, cdim, poly_order);  
  }
  else {
    maxwell->eqn.vol_term = CK(vol_kernels, cdim, poly_order);  
  }
  

  if ( use_conf_flux_surf ) {
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
}

struct gkyl_dg_eqn*
gkyl_dg_maxwell_cu_dev_inew(const struct gkyl_dg_maxwell_inp *inp)
{
  struct dg_maxwell *maxwell = (struct dg_maxwell*) gkyl_malloc(sizeof(struct dg_maxwell));

  // set basic parameters
  maxwell->eqn.num_equations = 8;
  if (inp->field_id == GKYL_FIELD_GR_D_B) assert(fabs(inp->lightSpeed - 1.0) < 1e-12);
  maxwell->maxwell_data.c = inp->lightSpeed;
  maxwell->maxwell_data.chi = inp->lightSpeed*inp->elcErrorSpeedFactor;
  maxwell->maxwell_data.gamma = inp->lightSpeed*inp->mgnErrorSpeedFactor;
  maxwell->gr_maxwell_data.chi = inp->elcErrorSpeedFactor;
  maxwell->gr_maxwell_data.gamma = inp->mgnErrorSpeedFactor;
  maxwell->gr_maxwell_data.K_phi = 0.0;
  maxwell->gr_maxwell_data.K_psi = 0.0;

  maxwell->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(maxwell->eqn.flags);
  maxwell->eqn.ref_count = gkyl_ref_count_init(gkyl_maxwell_free);

  maxwell->conf_flux_surf = 0;
  maxwell->lapse = 0;
  maxwell->shift = 0;
  maxwell->h_ij = 0;
  maxwell->h_ij_inv = 0;
  maxwell->det_h = 0;
  struct gkyl_array *conf_flux_surf = 0;
  struct gkyl_surf_and_vol_node_arrays *lapse = 0;
  struct gkyl_surf_and_vol_node_arrays *shift = 0;
  struct gkyl_surf_and_vol_node_arrays *h_ij = 0;
  struct gkyl_surf_and_vol_node_arrays *h_ij_inv = 0;
  struct gkyl_surf_and_vol_node_arrays *det_h = 0;
  maxwell->use_conf_flux_surf = false;
  if (inp->field_id == GKYL_FIELD_GR_D_B) {
    assert(inp->crange);
    maxwell->use_conf_flux_surf = true;
    maxwell->crange = *inp->crange;
    conf_flux_surf = gkyl_array_acquire(inp->conf_flux_surf);
    maxwell->conf_flux_surf = conf_flux_surf->on_dev;

    lapse = gkyl_surf_and_vol_node_arrays_acquire(inp->lapse);
    shift = gkyl_surf_and_vol_node_arrays_acquire(inp->shift);
    h_ij = gkyl_surf_and_vol_node_arrays_acquire(inp->h_ij);
    h_ij_inv = gkyl_surf_and_vol_node_arrays_acquire(inp->h_ij_inv);
    det_h = gkyl_surf_and_vol_node_arrays_acquire(inp->det_h);
    maxwell->lapse = lapse->on_dev;
    maxwell->shift = shift->on_dev;
    maxwell->h_ij = h_ij->on_dev;
    maxwell->h_ij_inv = h_ij_inv->on_dev;
    maxwell->det_h = det_h->on_dev;
  }

  // Configuration-space range + position-map Jacobian (device array) for the
  // mapped curl kernels. The Jacobian is required: callers without a
  // nonuniform position map pass an identity position map's Jacobian (J = 1
  // at every quadrature point), which keeps the kernels bit-identical to the
  // uniform-grid case.
  if (inp->crange) maxwell->crange = *inp->crange;
  assert(inp->jacob_pos);
  maxwell->jacob_pos = inp->jacob_pos->on_dev;

  // copy the host struct to device struct
  struct dg_maxwell *maxwell_cu = (struct dg_maxwell*) gkyl_cu_malloc(sizeof(struct dg_maxwell));
  gkyl_cu_memcpy(maxwell_cu, maxwell, sizeof(struct dg_maxwell), GKYL_CU_MEMCPY_H2D);
  dg_maxwell_set_cu_dev_ptrs<<<1,1>>>(maxwell_cu, inp->cbasis->b_type, inp->cbasis->ndim, 
    maxwell->use_conf_flux_surf, inp->cbasis->poly_order);

  // set parent on_dev pointer
  maxwell->eqn.on_dev = &maxwell_cu->eqn;

  // Host-side equation object should store host pointers.
  maxwell->conf_flux_surf = conf_flux_surf; 
  maxwell->lapse = lapse;
  maxwell->shift = shift;
  maxwell->h_ij = h_ij;
  maxwell->h_ij_inv = h_ij_inv;
  maxwell->det_h = det_h;

  return &maxwell->eqn;
}
