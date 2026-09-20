#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_vlasov_velocity_map.h>
#include <gkyl_vlasov_position_map.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_vlasov_conf_flux_surf gkyl_dg_vlasov_conf_flux_surf;

// Input packaged as a struct
struct gkyl_dg_vlasov_conf_flux_surf_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid. 
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions. 
  const struct gkyl_basis *phase_basis; // Phase-space basis functions. 
  const struct gkyl_range *hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  const struct gkyl_range *vel_range; // Velocity-space range for use in indexing velocity-space Jacobian. 
  const struct gkyl_vlasov_velocity_map *vel_map; // Velocity-space mapping object (REQUIRED; provides vmap/jacob_vel_surf).
  const struct gkyl_vlasov_position_map *pos_map; // Configuration-space mapping object (REQUIRED; provides jacob_pos).
  double skip_cell_thresh; // Phase-space density threshold for skipping cells in the Vlasov equation; by default no cells are skipped.
  enum gkyl_model_id model_id; // enum to determine what type of Vlasov model (e.g., non-relativistic vs. relativistic). 
  enum gkyl_hamil_id hamil_id; // enum for the Hamiltonian representation (sparse/dense velocity-space or phase-space expansion).
  bool use_lo; // bool to determine if using high-order kernels for non-canonical Hamiltonian models.
  bool use_gpu; // bool to determine if on GPU. 
};

/**
 * Create a new updater to compute the configuration-space fluxes in a modal DG scheme.
 * Updater computes the fluxes at the necessary Gauss-Legendre quadrature points to
 * eliminate aliasing errors by evaluating a total phase space characteristic, e.g., 
 * alpha_v = -grad_x(H) + q/m (E + grad_v(H) x B) - grad_x(q/m*phi + m*phi_g) + q/m*grad_v(H) x curl(A) - F_rad
 * at each quadrature point, and then computing F = alpha_v/2*(f^+ + f^-) - |alpha_v|/2*(f^+ - f^-)
 * at quadrature points. Here, f^+/- are the distribution function evaluated just outside/inside 
 * the cell interface and we utilize the full Lax-Friedrichs form to avoid warp divergence on GPUs. 
 * This nodal array is then converted back to a modal expansion on the *lower* surface of each cell
 * for use in the hyper_dg updater for integrating the surface terms in the DG discretization.    
 *
 * @param inp Input parameters defined in gkyl_dg_vlasov_conf_flux_surf_inp struct.
 * @return Pointer to configuration-space surface flux updater. 
 */
struct gkyl_dg_vlasov_conf_flux_surf* 
gkyl_dg_vlasov_conf_flux_surf_inew(const struct gkyl_dg_vlasov_conf_flux_surf_inp *inp);

/**
 * Create new updater to compute the configuration-space fluxes in a modal DG scheme on 
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_dg_vlasov_conf_flux_surf* 
gkyl_dg_vlasov_conf_flux_surf_cu_dev_inew(const struct gkyl_dg_vlasov_conf_flux_surf_inp *inp);

/**
 * Compute the configuration-space fluxes in a modal DG scheme. 
 *
 * @param up configuration-space surface flux updater. 
 * @param conf_range Configuration-space range for indexing electromagnetic fields. 
 * @param phase_range Phase-space range for indexing velocity-space flux array. 
 * @param poisson_tensor_conf Input configuration space expansion of the Poisson Tensor. 
 * @param hamil Input Hamiltonian for computing Hamiltonian forces.  
 * @param fin Input distribution function. 
 * @param cflrate Input cflrate. Accumulated to from maximum alpha_v at quadrature points.  
 * @param conf_flux_surf Output modal configuration-space fluxes. 
 */
void 
gkyl_dg_vlasov_conf_flux_surf_advance(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext, 
  const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil, 
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf);
 
 /**
 * Host-side wrapper for computing configuration-space fluxes on device. 
 */  
void 
gkyl_dg_vlasov_conf_flux_surf_advance_cu(struct gkyl_dg_vlasov_conf_flux_surf *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext, 
  const struct gkyl_array *poisson_tensor_conf, const struct gkyl_array *hamil, 
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *conf_flux_surf);

void
gkyl_dg_vlasov_conf_flux_surf_release(struct gkyl_dg_vlasov_conf_flux_surf* up);