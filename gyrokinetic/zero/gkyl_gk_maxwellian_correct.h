#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_velocity_map.h>

// Object type
typedef struct gkyl_gk_maxwellian_correct gkyl_gk_maxwellian_correct;

// input packaged as a struct
struct gkyl_gk_maxwellian_correct_inp {
  const struct gkyl_rect_grid *phase_grid; // Phase-space grid on which to compute moments
  const struct gkyl_basis *conf_basis; // Configuration-space basis functions
  const struct gkyl_basis *phase_basis; // Phase-space basis functions
  const struct gkyl_range *conf_range; // Configuration-space range
  const struct gkyl_range *conf_range_ext; // Extended configuration-space range (for internal memory allocations)
  const struct gkyl_range *vel_range; // velocity space range
  const struct gk_geometry *gk_geom; // Geometry object.
  const struct gkyl_velocity_map *vel_map; // Velocity space mapping object.
  double mass; // Mass factor
  bool bimaxwellian; // Bool for whether we are projecting a bi-Maxwellian instead of a Maxwellian.
  bool divide_jacobgeo; // Bool for whether to divide out the conf-space Jacobian from density.
  bool use_last_converged; // Boolean for if we are using the results of the iterative scheme
                           // *even if* the scheme fails to converge. 
  bool use_gpu; // Bool for gpu useage.
  double eps; // Tolerance for the iterator.
  int max_iter; // Number of total iterations.
};

// Correction status
struct gkyl_gk_maxwellian_correct_status {
  bool iter_converged; // true if iterations converged
  int num_iter; // number of iterations for the correction
  double error[4]; // error in each moment (n, u_par, T/m) or (n, u_par, T_par/m, T_perp/m)
};  

/**
 * Create new updater to correct the gyrokinetic Maxwellian (or bi-Maxwellian) 
 * distribution function so that its moments match desired input moments.
 *
 * @param inp Input parameters defined in gkyl_gk_maxwellian_correct_inp struct.
 * @return New updater pointer.
 */
struct gkyl_gk_maxwellian_correct* 
gkyl_gk_maxwellian_correct_inew(const struct gkyl_gk_maxwellian_correct_inp *inp);

/**
 * Fix the gyrokinetic Maxwellian (or bi-Maxwellian) distribution function
 * so that *all* its primitive moments (n, upar, T/m) or (n, upar, Tpar/m, Tperp/m) 
 * match target moments.
 * NOTE: If this algorithm fails, the returns the original distribution function
 * with only the desired density moment corrected.
 *
 * @param up gyrokinetic Maxwellian distribution function moment correction updater
 * @param f_max gyrokinetic Maxwellian (or bi-Maxwellian) distribution function to fix (modified in-place)
 * @param moms_target Target primitive moments (n, upar, T/m) or (n, upar, Tpar/m, Tperp/m)
 * @param phase_local Local phase-space range
 * @param conf_local Local configuration space range
 * @return Status of correction
 */
struct gkyl_gk_maxwellian_correct_status 
gkyl_gk_maxwellian_correct_all_moments(gkyl_gk_maxwellian_correct *up,
  struct gkyl_array *f_max, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local);

/**
 * Host-side wrapper for computing the absolute value of the 
 * difference in cell averages between the target moments and iterative moments
 * for correcting the Maxwellian or bi-Maxwellian.
 * 
 * @param conf_range Configuration-space range
 * @param num_comp Number of components being corrected (3 for Maxwellian, 4 for bi-Maxwellian)
 * @param nc Number of configuration space basis functions
 * @param moms_target Target primitive moments (n, u_par, T/m)
 * @param moms_iter Iterative moments used in fixed-point iteration
 * @param moms_abs_diff Absolute value of the difference between the cell averages of moms_iter and moms_target
 */
void gkyl_gk_maxwellian_correct_all_moments_abs_diff_cu(const struct gkyl_range *conf_range, 
  int num_comp, int nc, 
  const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *moms_abs_diff);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_gk_maxwellian_correct_release(gkyl_gk_maxwellian_correct *up);
