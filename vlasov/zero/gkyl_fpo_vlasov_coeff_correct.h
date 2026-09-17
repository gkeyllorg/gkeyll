#pragma once

#include <gkyl_fpo_vlasov_coeff_correct_priv.h>

// Object type
typedef struct gkyl_fpo_coeff_correct gkyl_fpo_coeff_correct;

/**
 * Create new updater to compute the necessary corrections to drag and diffusion 
 * coefficients to enforce momentum and energy conservation in the FPO.
 *
*/ 
gkyl_fpo_coeff_correct*
gkyl_fpo_coeff_correct_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_rng, bool use_gpu);

gkyl_fpo_coeff_correct*
gkyl_fpo_coeff_correct_cu_dev_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_range *conf_rng);

/**
 * Solve linear system to calculate corrections to drag and diffusion coefficients,
 * and accumulate them onto the existing drag and diffusion coefficient arrays.
 *
*/
void gkyl_fpo_coeff_correct_advance(gkyl_fpo_coeff_correct *up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  const struct gkyl_array *fpo_moms, const struct gkyl_array *boundary_corrections,
  const struct gkyl_array *moms, struct gkyl_array *drag_diff_coeff_corrs,
  struct gkyl_array *drag_coeff, struct gkyl_array *drag_coeff_surf,
  struct gkyl_array *diff_coeff, struct gkyl_array *diff_coeff_surf, bool use_gpu);

void gkyl_fpo_coeff_correct_advance_cu(gkyl_fpo_coeff_correct *up,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  const struct gkyl_array *fpo_moms, const struct gkyl_array *boundary_corrections,
  const struct gkyl_array *moms, struct gkyl_array *drag_diff_coeff_corrs,
  struct gkyl_array *drag_coeff, struct gkyl_array *drag_coeff_surf,
  struct gkyl_array *diff_coeff, struct gkyl_array *diff_coeff_surf);

/**
 * Delete updater.
 *
 * @param up Updater to delete
 */
void gkyl_fpo_vlasov_coeff_correct_release(gkyl_fpo_coeff_correct *up);
