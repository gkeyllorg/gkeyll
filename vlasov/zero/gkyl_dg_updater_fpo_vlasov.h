#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_collisions gkyl_dg_updater_collisions;

// return type for drag and diffusion timers
struct gkyl_dg_updater_fpo_vlasov_tm {
  double drag_tm, diff_tm; // time for drag and diffusion updates
};

/**
 * Create new updater to update FPO equations using hyper dg and hyper dg gen stencil (for diffusion).
 *
 * @param grid Grid object
 * @param pbasis Phase-space basis function
 * @param phase_range Phase space range
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 * @return New fpo updater object
 */
struct gkyl_dg_updater_collisions* 
gkyl_dg_updater_fpo_vlasov_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *pbasis, const struct gkyl_range *phase_range, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param fpo fpo updater object
 * @param update_rng Range on which to compute.
 * @param drag_coeff Drag coefficient vector calculated from first Rosenbluth potential
 * @param diff_coeff Diffusion tensor calculated from second Rosenbluth potential
 * @param drag_coeff_surf Surface expansion of recovered drag coeff at lower cell boundary
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_fpo_vlasov_advance(struct gkyl_dg_updater_collisions *fpo,
  long offsets[36],
  const struct gkyl_range *update_rng,
  const struct gkyl_array *drag_coeff, 
  const struct gkyl_array *drag_coeff_surf,
  const struct gkyl_array *sgn_drag_coeff_surf,
  const struct gkyl_array *const_sgn_drag_coeff_surf,
  const struct gkyl_array *diff_coeff,
  const struct gkyl_array *diff_coeff_surf,
  const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in drag and diffusion terms
 *
 * @param lbo Updater object
 * @return timers
 */
struct gkyl_dg_updater_fpo_vlasov_tm gkyl_dg_updater_fpo_vlasov_get_tm(const struct gkyl_dg_updater_collisions *coll);

/**
 * Delete updater.
 *
 * @param lbo Updater to delete.
 */
void gkyl_dg_updater_fpo_vlasov_release(struct gkyl_dg_updater_collisions* coll);
