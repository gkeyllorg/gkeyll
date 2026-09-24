// Private header: not for direct use
#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

struct gkyl_vlasov_lte_correct
{
  int num_conf_basis; // Number of configuration-space basis functions
  int num_comp; // Number of components being corrected vdim+2 (n, V_drift, T/m)

  struct gkyl_velocity_map *vel_map; // Velocity space mapping object.

  struct gkyl_array *moms_iter;
  struct gkyl_array *d_moms;
  struct gkyl_array *dd_moms;

  struct gkyl_vlasov_lte_moments *moments_up;
  struct gkyl_vlasov_lte_proj_on_basis *proj_lte;

  // error estimate (n, V_drift, T/m), 0 - success., num. picard iterations
  double *error; // absolute value of difference in cell averages between iteration and target
  double eps; // tolerance for the iterator
  int max_iter; // number of total iterations
  bool use_last_converged; // Boolean for if we are using the results of the iterative scheme
                           // *even if* the scheme fails to converge. 

  bool use_gpu; // Boolean if we are performing projection on device.
  double *error_cu; // error on device if using GPUs 
  struct gkyl_array *abs_diff_moms;
};
