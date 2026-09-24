#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_diffusion_gyrokinetic gkyl_dg_updater_diffusion_gyrokinetic;

// return type for drag and diffusion timers
struct gkyl_dg_updater_diffusion_gyrokinetic_tm {
  double diffusion_tm; // time for diffusion updates
};

/**
 * Create new updater to update gyrokinetic diffusion equations using hyper dg.
 *
 * @param grid Grid object.
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param is_diff_constant If diffusion coefficient constant or spatially constant.
 * @param diff_in_dir Whether to apply diffusion in each direction.
 * @param diff_order Diffusion order.
 * @param diff_range Range object to index the diffusion coefficient.
 * @param is_zero_flux_bc True for boundaries with zero flux BCs.
 * @param coeff Diffusion coefficient/tensor.
 * @param jacobgeo_inv Reciprocal of the configuration space jacobian.
 * @param use_gpu Whether to run on host or device.
 * @return New diff updater object
 */
struct gkyl_dg_updater_diffusion_gyrokinetic* gkyl_dg_updater_diffusion_gyrokinetic_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis, bool is_diff_const, const bool *diff_in_dir,
  int diff_order, const struct gkyl_range *diff_range, const bool *is_zero_flux_bc,
  const struct gkyl_array *coeff, const struct gkyl_array *jacobgeo_inv, bool use_gpu);

/**
 * Acquire gyrokinetic diffusion equation object.
 *
 * @param up gyrokinetic diffusion updater object.
 * 
 * @return gyrokinetic diffusion equation object.
 */
struct gkyl_dg_eqn* 
gkyl_dg_updater_diffusion_gyrokinetic_acquire_eqn(const gkyl_dg_updater_diffusion_gyrokinetic* up);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up diffusion updater object.
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater.
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_diffusion_gyrokinetic_advance(struct gkyl_dg_updater_diffusion_gyrokinetic *up,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in diffusion terms
 *
 * @param diffusion Updater object
 * @return timers
 */
struct gkyl_dg_updater_diffusion_gyrokinetic_tm gkyl_dg_updater_diffusion_gyrokinetic_get_tm(const struct gkyl_dg_updater_diffusion_gyrokinetic *up);

/**
 * Delete updater.
 *
 * @param diffusion Updater to delete.
 */
void gkyl_dg_updater_diffusion_gyrokinetic_release(struct gkyl_dg_updater_diffusion_gyrokinetic *up);
