#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_gr_maxwell_surf_and_vol_nodes.h>

// Object type
typedef struct gkyl_dg_gr_maxwell_conf_flux_surf gkyl_dg_gr_maxwell_conf_flux_surf;

// Input packaged as a struct
struct gkyl_dg_gr_maxwell_conf_flux_surf_inp {
  const struct gkyl_rect_grid *conf_grid; // Configuration-space grid. 
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions. 
  const struct gkyl_vlasov_position_map *pos_map; // Configuration-space mapping object (REQUIRED; provides the cell-constant jacob_pos).
  enum gkyl_field_id field_id; // enum to determine what type of Maxwell model (e.g., curved space non-relativistic vs. relativistic). 
  const int *theta_pole_lo; // (lower edge) Default zeros, but 1 if any directions use theta-pole BC's.
  const int *theta_pole_up; // (upper edge) Default zeros, but 1 if any directions use theta-pole BC's.
  double chi; // Electric-field divergence correction speed.
  double gamma; // Magnetic-field divergence correction speed.
  bool use_gpu; // bool to determine if on GPU. 
};

/**
 * Create a new updater to compute the configuration-space fluxes in a modal DG scheme.
 * Updater computes the fluxes at the necessary Gauss-Legendre quadrature points.
 * 
 * Here, fields are evaluated just outside/inside 
 * the cell interface and we utilize the full Lax-Friedrichs form to avoid warp divergence on GPUs. 
 * This nodal array is then converted back to a modal expansion on the *lower* surface of each cell
 * for use in the hyper_dg updater for integrating the surface terms in the DG discretization.    
 *
 * @param inp Input parameters defined in gkyl_dg_gr_maxwell_conf_flux_surf_inp struct.
 * @return Pointer to configuration-space surface flux updater. 
 */
struct gkyl_dg_gr_maxwell_conf_flux_surf* 
gkyl_dg_gr_maxwell_conf_flux_surf_inew(const struct gkyl_dg_gr_maxwell_conf_flux_surf_inp *inp);

/**
 * Create new updater to compute the configuration-space fluxes in a modal DG scheme on 
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_gr_maxwell_conf_flux_surf* 
gkyl_dg_gr_maxwell_conf_flux_surf_cu_dev_inew(const struct gkyl_dg_gr_maxwell_conf_flux_surf_inp *inp);

/**
 * Compute the configuration-space fluxes in a modal DG scheme. 
 *
 * @param up configuration-space surface flux updater. 
 * @param conf_range Configuration-space range for indexing electromagnetic fields. 
 * @param conf_range_ext Configuration-space range for indexing electromagnetic fields (extended). 
 * @param lapse Lapse (ADM alpha) at the nodal surface points. 
 * @param shift shift (ADM \beta^i - contravaraint) at the nodal surface points. 
 * @param h_ij Spatial metric (ADM h_ij - covariant) at the nodal surface points. 
 * @param h_ij_inv Inverse spatial metric (ADM h^ij - contravariant) at the nodal surface points.
 * @param det_h Spatial metric Jacobian (Jc = det_h = sqrt(det(h_ij))) at the nodal surface points.  
 * @param field_con Input contravariant (JD,JB) field components. 
 * @param field_no_J_con Input contravariant (D,B) field components (without determinant factor Jc) Used for edges.
 * @param cflrate Input cflrate. Accumulated to from maximum alpha_v at quadrature points.  
 * @param conf_flux_surf Output modal configuration-space fluxes. 
 */
void 
gkyl_dg_gr_maxwell_conf_flux_surf_advance(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_range_ext, 
  const struct gkyl_surf_and_vol_node_arrays *lapse, const struct gkyl_surf_and_vol_node_arrays *shift, 
  const struct gkyl_surf_and_vol_node_arrays *h_ij, const struct gkyl_surf_and_vol_node_arrays *h_ij_inv,
  const struct gkyl_surf_and_vol_node_arrays *det_h, 
  const struct gkyl_array *field_con, const struct gkyl_array *field_no_J_con, struct gkyl_array *cflrate, 
  struct gkyl_array *conf_flux_surf);

 /**
 * Host-side wrapper for computing configuration-space fluxes on device. 
 */  
void 
gkyl_dg_gr_maxwell_conf_flux_surf_advance_cu(struct gkyl_dg_gr_maxwell_conf_flux_surf *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_range_ext, 
  const struct gkyl_surf_and_vol_node_arrays *lapse, const struct gkyl_surf_and_vol_node_arrays *shift, 
  const struct gkyl_surf_and_vol_node_arrays *h_ij, const struct gkyl_surf_and_vol_node_arrays *h_ij_inv,
  const struct gkyl_surf_and_vol_node_arrays *det_h,
  const struct gkyl_array *field_con, const struct gkyl_array *field_no_J_con, struct gkyl_array *cflrate, 
  struct gkyl_array *conf_flux_surf);

void
gkyl_dg_gr_maxwell_conf_flux_surf_release(struct gkyl_dg_gr_maxwell_conf_flux_surf* up);
