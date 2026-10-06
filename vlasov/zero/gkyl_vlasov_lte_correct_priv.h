// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Which cells to drop from the correction when evaluating the per-cell errors. Cells whose
// iterate has non-positive density or temperature are dropped in every mode but NONE.
enum vlasov_lte_correct_drop_mode {
  VLASOV_LTE_CORRECT_DROP_NONE = 0, // Report the errors only; drop no cell.
  VLASOV_LTE_CORRECT_DROP_NONPOSITIVE, // Drop cells that lost positivity.
  VLASOV_LTE_CORRECT_DROP_ABOVE_TOL, // Also drop cells whose error exceeds the tolerance.
  VLASOV_LTE_CORRECT_DROP_NOT_IMPROVED, // Also drop cells no closer to the target than the uncorrected projection.
};

struct gkyl_vlasov_lte_correct {
  int num_conf_basis; // Number of configuration-space basis functions
  int num_comp; // Number of components being corrected vdim+2 (n, V_drift, T/m)

  struct gkyl_velocity_map *vel_map; // Velocity space mapping object.

  struct gkyl_array *moms_iter; // Moments of the iterate, then the parameters of the next projection.
  struct gkyl_array *d_moms; // Accumulated correction to the target moments.
  struct gkyl_array *dd_moms; // Mismatch between the target and the iterate.
  struct gkyl_array
    *corr_mask; // Per-cell scalar: 1 while the cell is being corrected, 0 once dropped.
  struct gkyl_array *abs_diff_moms; // Per-cell error of each moment (0 in dropped cells).
  struct gkyl_array *abs_diff_init; // Per-cell error of the uncorrected projection.

  struct gkyl_vlasov_lte_moments *moments_up;
  struct gkyl_vlasov_lte_proj_on_basis *proj_lte;

  double *error; // Maximum error of each moment over the cells being corrected.
  double eps; // tolerance for the iterator
  int max_iter; // number of total iterations
  bool use_last_converged; // Keep the last iterate in cells that did not converge but improved.

  bool use_gpu; // Boolean if we are performing projection on device.
  double *error_cu; // error on device if using GPUs
  double *mask_sum_cu; // Sum of the mask on device if using GPUs.
};

// Error of the iterate against the target in one cell: relative error of the cell averages of
// the density and temperature, absolute error of drift components below one (the drift may
// vanish) and relative error above. Drops the cell according to drop_mode by zeroing its mask;
// dropped cells report zero error so that they do not enter the convergence test.
GKYL_CU_DH static inline void
vlasov_lte_correct_cell_errors(
  int num_comp, int nc, double tol, int drop_mode, const double *moms_target,
  const double *moms_iter, const double *abs_diff_init, double *mask, double *abs_diff
)
{
  double err[GKYL_MAX_VDIM + 2];
  double err_max = 0.0, init_max = 0.0;
  int T_idx = num_comp - 1; // T/m is always the last component
  for (int c = 0; c < num_comp; ++c) {
    double diff = fabs(moms_iter[c * nc] - moms_target[c * nc]);
    bool relative = (c == 0) || (c == T_idx) || (fabs(moms_target[c * nc]) >= 1.0);
    err[c] = relative ? diff / fabs(moms_target[c * nc]) : diff;
    err_max = fmax(err_max, err[c]);
    init_max = fmax(init_max, abs_diff_init[c]);
  }

  if (drop_mode != VLASOV_LTE_CORRECT_DROP_NONE) {
    // A NaN error or moment fails these comparisons and drops the cell.
    bool keep = (moms_iter[0] > 0.0) && (moms_iter[T_idx * nc] > 0.0) && (err_max == err_max);
    if (drop_mode == VLASOV_LTE_CORRECT_DROP_ABOVE_TOL) {
      keep = keep && (err_max <= tol);
    } else if (drop_mode == VLASOV_LTE_CORRECT_DROP_NOT_IMPROVED) {
      keep = keep && (err_max <= init_max);
    }
    if (!keep) {
      mask[0] = 0.0;
    }
  }

  for (int c = 0; c < num_comp; ++c) {
    abs_diff[c] = mask[0] > 0.0 ? err[c] : 0.0;
  }
}
