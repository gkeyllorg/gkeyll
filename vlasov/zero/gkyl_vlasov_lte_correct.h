#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_vlasov_velocity_map.h>

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
  const struct gkyl_range
    *conf_range_ext; // Extended configuration-space range (for internal memory allocations)
  const struct gkyl_range *vel_range; // velocity space range
  const struct gkyl_range *phase_range; // phase space range
  const struct gkyl_vlasov_velocity_map
    *vel_map; // Velocity-space mapping object. NULL => uniform velocity grid.
  const struct gkyl_range *
    hamil_range; // Range for indexing Hamiltonian (either velocity-space range or full phase-space range).
  const struct gkyl_array *hamil; // (Can-pb quantitiy) Hamiltonian
  enum gkyl_model_id model_id; // Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h)
  enum gkyl_hamil_id
    hamil_id; // Enum for the Hamiltonian representation (sparse/dense velocity-space or phase-space expansion).
  const struct gkyl_array *gamma_inv; // SR quantitiy: 1/gamma = 1/sqrt(1 + p^2)
  const struct gkyl_array *h_ij; // (Can-pb quantitiy) metric tensor (covariant components)
  const struct gkyl_array
    *h_ij_inv; // (Can-pb quantitiy) inverse metric tensor (contravariant components)
  const struct gkyl_array *det_h; // (Can-pb quantitiy) determinant of the metric tensor
  enum gkyl_quad_type quad_type; // type of quadrature to use: defaults to Gaussian
  bool use_last_converged; // In cells that do not converge within max_iter, keep the last
    // iterate when it is closer to the target than the uncorrected projection (true), or
    // always fall back to the uncorrected projection (false).
  bool
    use_extended_hamil_def; // bool to determine if we are using the extended canonical-pb hamiltonian
  const struct gkyl_array
    *background_flows; // Can-pb quantity: Background flows for extended hamiltonians
  const struct gkyl_array *effective_potential; // Can-pb quantity: specified effective potential
  bool use_gpu; // bool for gpu usage
  double eps; // tolerance for the iterator
  int max_iter; // number of total iterations
};

// Correction status
struct gkyl_vlasov_lte_correct_status {
  bool iter_converged; // 0 when every cell converged to the tolerance, 1 otherwise (this
    // is the value the apps write to the correction-status diagnostic).
  int num_iter; // number of iterations for the correction
  int num_cells_dropped; // number of cells that fell back to the uncorrected projection
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
struct gkyl_vlasov_lte_correct *gkyl_vlasov_lte_correct_inew(
  const struct gkyl_vlasov_lte_correct_inp *inp
);

/**
 * Fix the LTE (local thermodynamic equlibrium) distribution function
 * (Maxwellian for non-relativistic/Maxwell-Juttner for relativistic)
 * so that *all* its stationary-frame moments (n, V_drift, T/m) match target moments.
 * The correction is a fixed-point iteration carried out cell by cell. A cell whose iterate
 * loses positivity, or that has not converged after max_iter iterations (see use_last_converged),
 * is dropped from the correction and returned as the projection of its target moments with
 * only the density corrected; the other cells keep their corrected distribution.
 *
 * @param up LTE distribution function moment correction updater
 * @param f_lte LTE distribution function to fix (modified in-place)
 * @param moms_target Target stationary-frame moments (n, V_drift, T/m)
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 * @return Status of correction
 */
struct gkyl_vlasov_lte_correct_status gkyl_vlasov_lte_correct_all_moments(
  gkyl_vlasov_lte_correct *up, struct gkyl_array *f_lte, const struct gkyl_array *moms_target,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local
);

/**
 * Host-side wrapper for computing the per-cell errors of the iterate against the target
 * moments and dropping cells from the correction (see the private header).
 */
void gkyl_vlasov_lte_correct_cell_errors_cu(
  const struct gkyl_range *conf_range, int num_comp, int nc, double tol, int drop_mode,
  const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter,
  const struct gkyl_array *abs_diff_init, struct gkyl_array *corr_mask,
  struct gkyl_array *abs_diff_moms
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_vlasov_lte_correct_release(gkyl_vlasov_lte_correct *up);
