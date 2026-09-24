#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h> 

// Object type
typedef struct gkyl_vlasov_lte_moments gkyl_vlasov_lte_moments;

// input packaged as a struct
struct gkyl_vlasov_lte_moments_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid on which to compute moments
  const struct gkyl_rect_grid *vel_grid; // Velocity-space grid
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions
  const struct gkyl_basis *vel_basis; // Velocity-space basis functions
  const struct gkyl_basis *phase_basis; // Phase-space basis functions
  const struct gkyl_range *conf_range; // Configuration-space range
  const struct gkyl_range *conf_range_ext; // Extended configuration-space range (for internal memory allocations)
  const struct gkyl_range *vel_range; // Velocity-space range
  const struct gkyl_range *phase_range; // Phase-space range
  const struct gkyl_array *gamma; // SR quantitiy: gamma = sqrt(1 + p^2)
  const struct gkyl_array *gamma_inv; // SR quantitiy: 1/gamma = 1/sqrt(1 + p^2)
  const struct gkyl_array *h_ij; // Can-pb quantity: metric tensor (covariant components)
  const struct gkyl_array *h_ij_inv; // Can-pb quantity: Inverse metric tensor (contravaraint components)
  const struct gkyl_array *det_h; // Can-pb quantity: determinant of the metric tensor
  const struct gkyl_array *hamil; // Can-pb quantity: hamiltonian
  enum gkyl_model_id model_id; // Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h)
  double mass; // Mass factor 
  bool use_gpu; // bool for gpu useage
};


/**
 * Create new updater to compute the moments for the equivalent LTE (local thermodynamic equlibrium) 
 * distribution function (Maxwellian for non-relativistic/Maxwell-Juttner for relativistic)
 * Updater always returns (n, V_drift, T/m) where n and T/m are the stationary frame 
 * density and temperature/mass (the frame moving at velocity V_drift).
 * 
 * @param inp Input parameters defined in gkyl_vlasov_lte_moments_inp struct.
 * @return New updater pointer.
 */
struct gkyl_vlasov_lte_moments*
gkyl_vlasov_lte_moments_inew(const struct gkyl_vlasov_lte_moments_inp *inp);

/**
 * Compute the density moments of an arbitrary distribution function for the equivalent 
 * LTE (local thermodynamic equlibrium) distribution function.
 * (Maxwellian for non-relativistic/Maxwell-Juttner for relativistic)
 * Computes n, the stationary frame density (the frame moving at velocity V_drift).
 *
 * @param lte_moms LTE moments updater
 * @param phase_local Phase-space range on which to compute moments.
 * @param conf_local Configuration-space range on which to compute moments.
 * @param fin Input distribution function
 * @param density Output stationary-frame density
 */
void gkyl_vlasov_lte_density_moment_advance(struct gkyl_vlasov_lte_moments *lte_moms, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *density);

/**
 * Compute the moments of an arbitrary distribution function for the equivalent 
 * LTE (local thermodynamic equlibrium) distribution function.
 * (Maxwellian for non-relativistic/Maxwell-Juttner for relativistic)
 * Computes (n, V_drift, T/m) where n and T/m are the stationary frame 
 * density and temperature/mass (the frame moving at velocity V_drift).
 *
 * @param lte_moms LTE moments updater
 * @param phase_local Phase-space range on which to compute moments.
 * @param conf_local Configuration-space range on which to compute moments.
 * @param fin Input distribution function
 * @param moms Output LTE moments (n, V_drift, T/m)
 */
void gkyl_vlasov_lte_moments_advance(struct gkyl_vlasov_lte_moments *lte_moms, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local, 
  const struct gkyl_array *fin, struct gkyl_array *moms);

/**
 * Delete updater.
 *
 * @param lte_moms Updater to delete.
 */
void gkyl_vlasov_lte_moments_release(gkyl_vlasov_lte_moments* lte_moms);
