#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dist_type.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

// Object type
typedef struct gkyl_vlasov_dist_correct gkyl_vlasov_dist_correct;

// input packaged as a struct
struct gkyl_vlasov_dist_correct_inp {
  const struct gkyl_dist_proj_type *dist_proj; // Distribution projection object
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid on which to compute moments
  const struct gkyl_rect_grid *vel_grid; // Velocity-space grid
  const struct gkyl_basis *conf_basis; // Configuration space basis functions
  const struct gkyl_basis *vel_basis; // Velocity-space basis functions
  const struct gkyl_basis *phase_basis; // Phase-space basis functions
  const struct gkyl_range *conf_range; // Configuration-space range
  const struct gkyl_range *conf_range_ext; // Extended configuration-space range (for internal memory allocations)
  const struct gkyl_range *vel_range; // velocity space range
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  const struct gkyl_range *phase_range; // phase space range
  enum gkyl_model_id model_id; // Enum identifier for model type (e.g., SR, see gkyl_eqn_type.h)
  double mass;
  enum gkyl_quad_type quad_type; // type of quadrature to use: defaults to Gaussian
  bool use_last_converged; // Boolean for if we are using the results of the iterative scheme
    // *even if* the scheme fails to converge.
  bool use_gpu; // bool for gpu usage
  double eps; // tolerance for the iterator
  int max_iter; // number of total iterations
};

// Correction status
struct gkyl_vlasov_dist_correct_status {
  bool iter_converged; // true if iterations converged
  int num_iter; // number of iterations for the correction
  double *error; // error in each moment
};

/**
 * Create new updater to correct a projected distribution function
 * so that its moments match desired input moments.
 *
 * @param inp Input parameters defined in gkyl_vlasov_dist_correct_inp struct.
 * @return New updater pointer.
 */
struct gkyl_vlasov_dist_correct *gkyl_vlasov_dist_correct_inew(
  const struct gkyl_vlasov_dist_correct_inp *inp
);

/**
 * Correct the distribution function so that its moments match
 * the target distribution moments.
 *
 * @param up Distribution function moment correction updater
 * @param f_dist Distribution function to correct (modified in-place)
 * @param moms_target Target distribution moments
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 * @return Status of correction
 */
struct gkyl_vlasov_dist_correct_status gkyl_vlasov_dist_correct_all_moments(
  gkyl_vlasov_dist_correct *up, struct gkyl_array *f_dist, const struct gkyl_array *moms_target,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local
);

/**
 * Host-side wrapper for computing the absolute value of the 
 * difference in cell averages between the target moments and iterative moments.
 */
void
gkyl_vlasov_dist_correct_all_moments_abs_diff_cu(
  const struct gkyl_range *conf_range, int num_mom, int nc, int vdim, const struct gkyl_array *moms_target,
  const struct gkyl_array *moms_iter, struct gkyl_array *moms_abs_diff
);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void
gkyl_vlasov_dist_correct_release(gkyl_vlasov_dist_correct *up);
