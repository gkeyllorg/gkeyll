// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_maxwellian_moments.h>
#include <gkyl_gk_maxwellian_proj_on_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct gkyl_gk_maxwellian_correct
{
  int num_conf_basis; // Number of configuration-space basis functions
  int num_comp; // Number of components being corrected 
                // 3 for a Maxwellian (n, upar, T/m), 4 for a bi-Maxwellian (n, upar, Tpar/m, Tperp/m)
  bool bimaxwellian; // Bool for whether we are correcting a bi-Maxwellian's moments.

  struct gkyl_velocity_map *vel_map; // Velocity space mapping object.

  struct gkyl_array *moms_iter;
  struct gkyl_array *d_moms;
  struct gkyl_array *dd_moms;

  struct gkyl_gk_maxwellian_moments *moments_up;
  struct gkyl_gk_maxwellian_proj_on_basis *proj_max;

  // error estimate (n, u_par, T/m), or (n, u_par, T_par/m, T_perp/m) 0 - success., num. picard iterations
  double *error; // absolute value of difference in cell averages between iteration and target
  double eps; // tolerance for the iterator
  int max_iter; // number of total iterations
  bool use_last_converged; // Boolean for if we are using the results of the iterative scheme
                           // *even if* the scheme fails to converge. 

  bool use_gpu; // Boolean if we are performing projection on device.
  double *error_cu; // error on device if using GPUs 
  struct gkyl_array *abs_diff_moms;
};
