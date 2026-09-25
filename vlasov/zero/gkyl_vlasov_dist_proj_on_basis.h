#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>
#include <gkyl_dist_type.h>

// Object type
typedef struct gkyl_vlasov_dist_proj_on_basis gkyl_vlasov_dist_proj_on_basis;

// input packaged as a struct
struct gkyl_vlasov_dist_proj_on_basis_inp {
  const struct gkyl_dist_proj_type *dist_proj; // Generic distribution projection object
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid on which to compute moments
  const struct gkyl_rect_grid *vel_grid; // Velocity-space grid
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions
  const struct gkyl_basis *vel_basis; // Velocity-space basis functions
  const struct gkyl_basis *phase_basis; // Phase-space basis functions
  const struct gkyl_range *conf_range; // Configuration-space range
  const struct gkyl_range *conf_range_ext; // Extended configuration-space range (for internal memory allocations)
  const struct gkyl_range *vel_range; // velocity space range
  const struct gkyl_range *phase_range; // phase space range
  const struct gkyl_array *gamma; // SR quantitiy: gamma = sqrt(1 + p^2)
  const struct gkyl_array *gamma_inv; // SR quantitiy: 1/gamma = 1/sqrt(1 + p^2)
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  const struct gkyl_array *h_ij; // (Can-bp quantity) metric tensor (covariant components)
  const struct gkyl_array *h_ij_inv; // (Can-bp quantity) inverse of the metric tensor (contravariant components)
  const struct gkyl_array *det_h; // (Can-bp quantity) determinant of the metric tensor
  const struct gkyl_array *hamil; // (Can-bp quantity) Hamiltonian
  enum gkyl_model_id model_id; // Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h)
  double mass; // Mass factor
  bool use_gpu; // bool for gpu useage

  enum gkyl_quad_type quad_type;
};

/**
 * Create new updater to project a selected distribution function onto basis functions. 
 * The distribution is selected through inp->dist_proj.
 * Free using gkyl_vlasov_dist_proj_on_basis_release method.
 *
 * @param inp Input parameters defined in gkyl_vlasov_dist_proj_on_basis_inp struct.
 * @return New updater pointer.
 */
struct gkyl_vlasov_dist_proj_on_basis* 
gkyl_vlasov_dist_proj_on_basis_inew(const struct gkyl_vlasov_dist_proj_on_basis_inp *inp);

/**
 * Evaluate the distribution parameters at configuration-space quadrature nodes.
 *
 * @param up Project on basis updater to run
 * @param conf_range Configuration-space range
 * @param moms_dist Modal distribution parameters
 */
void gkyl_vlasov_dist_proj_on_basis_moments_advance(gkyl_vlasov_dist_proj_on_basis *up,
  const struct gkyl_range *conf_range,
  const struct gkyl_array *moms_dist);

/**
 * Compute projection of a selected distribution on basis. 
 * The distribution parameters must first be evaluated using
 * gkyl_vlasov_dist_proj_on_basis_moments_advance.
 * Further note: We utilize the gkyl_correct_density_moment_vlasov_dist (see gkyl_correct_dist.h)
 * to correct the density moment of the projected distribution before we return f_dist.
 *
 * @param up Project on basis updater to run
 * @param phase_range Phase-space range
 * @param conf_range Configuration-space range
 * @param f_dist Output distribution function
 */
void gkyl_vlasov_dist_proj_on_basis_advance(gkyl_vlasov_dist_proj_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
  struct gkyl_array *f_dist);

/**
 * Host-side wrapper for initial canonical-pb vars
 */
void gkyl_vlasov_dist_proj_on_basis_geom_quad_vars_cu(gkyl_vlasov_dist_proj_on_basis *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array *h_ij,
  const struct gkyl_array *h_ij_inv, const struct gkyl_array *det_h);

/**
 * Host-side wrapper for projection of the selected distribution function on device
 */
void gkyl_vlasov_dist_proj_on_basis_advance_cu(gkyl_vlasov_dist_proj_on_basis *up,
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range,
 struct gkyl_array *f_dist);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_vlasov_dist_proj_on_basis_release(gkyl_vlasov_dist_proj_on_basis* up);
