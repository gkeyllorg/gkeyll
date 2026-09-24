#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

// Object type
typedef struct gkyl_vlasov_lte_correct gkyl_vlasov_lte_correct;

// input packaged as a struct
struct gkyl_vlasov_lte_correct_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid on which to compute moments
  const struct gkyl_rect_grid *vel_grid; // Velocity-space grid
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions
  const struct gkyl_basis *vel_basis; // Velocity-space basis functions
  const struct gkyl_basis *phase_basis; // Phase-space basis functions
  const struct gkyl_range *conf_range; // Configuration-space range
  const struct gkyl_range *conf_range_ext; // Extended configuration-space range (for internal memory allocations)
  const struct gkyl_range *vel_range; // velocity space range
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  const struct gkyl_range *phase_range; // phase space range
  const struct gkyl_array *gamma; // SR quantitiy: gamma = sqrt(1 + p^2)
  const struct gkyl_array *gamma_inv; // SR quantitiy: 1/gamma = 1/sqrt(1 + p^2)
  const struct gkyl_array *h_ij; // (Can-pb quantitiy) metric tensor (covariant components)
  const struct gkyl_array *h_ij_inv; // (Can-pb quantitiy) inverse metric tensor (contravariant components) 
  const struct gkyl_array *det_h; // (Can-pb quantitiy) determinant of the metric tensor 
  const struct gkyl_array *hamil; // (Can-pb quantitiy) Hamiltonian
  enum gkyl_model_id model_id; // Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h)
  enum gkyl_quad_type quad_type; // type of quadrature to use: defaults to Gaussian
  bool use_last_converged; // Boolean for if we are using the results of the iterative scheme
                           // *even if* the scheme fails to converge. 
  bool use_gpu; // bool for gpu usage
  double eps; // tolerance for the iterator
  int max_iter; // number of total iterations
};

// Correction status
struct gkyl_vlasov_lte_correct_status {
  bool iter_converged; // true if iterations converged
  int num_iter; // number of iterations for the correction
  double error[5]; // error in each moment, up to 5 (vdim+2) components
};  

/**
 * Create new updater to correct the LTE (local thermodynamic equlibrium) distribution 
 * function (Maxwellian for non-relativistic/Maxwell-Juttner for relativistic)
 * so that its moments match desired input moments.
 *
 * @param inp Input parameters defined in gkyl_vlasov_lte_correct_inp struct.
 * @return New updater pointer.
 */
struct gkyl_vlasov_lte_correct* 
gkyl_vlasov_lte_correct_inew(const struct gkyl_vlasov_lte_correct_inp *inp);

/**
 * Fix the LTE (local thermodynamic equlibrium) distribution function
 * (Maxwellian for non-relativistic/Maxwell-Juttner for relativistic)
 * so that *all* its stationary-frame moments (n, V_drift, T/m) match target moments.
 * NOTE: If this algorithm fails, the returns the original distribution function
 * with only the desired stationary-frame density moment corrected.
 *
 * @param up LTE distribution function moment correction updater
 * @param f_lte LTE distribution function to fix (modified in-place)
 * @param moms_target Target stationary-frame moments (n, V_drift, T/m)
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 * @return Status of correction
 */
struct gkyl_vlasov_lte_correct_status gkyl_vlasov_lte_correct_all_moments(gkyl_vlasov_lte_correct *up,
  struct gkyl_array *f_lte, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Host-side wrapper for computing the absolute value of the 
 * difference in cell averages between the target moments and iterative moments.
 */
void gkyl_vlasov_lte_correct_all_moments_abs_diff_cu(const struct gkyl_range *conf_range, 
  int num_comp, int nc, 
  const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *moms_abs_diff);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_vlasov_lte_correct_release(gkyl_vlasov_lte_correct *up);
