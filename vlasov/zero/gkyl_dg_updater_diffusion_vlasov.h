#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_diffusion_vlasov gkyl_dg_updater_diffusion_vlasov;

// return type for drag and diffusion timers
struct gkyl_dg_updater_diffusion_vlasov_tm {
  double diffusion_tm; // time for diffusion updates
};

/**
 * Create new updater to update vlasov diffusion equations using hyper dg.
 *
 * @param grid Grid object.
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param is_diff_constant If diffusion coefficient constant or spatially constant.
 * @param diff_in_dir Whether to apply diffusion in each direction.
 * @param diff_order Diffusion order.
 * @param diff_range Range object to index the diffusion coefficient.
 * @param is_zero_flux_dir True in directions with (lower and upper) zero flux BCs.
 * @param use_gpu Whether to run on host or device.
 * @return New diff updater object
 */
struct gkyl_dg_updater_diffusion_vlasov* gkyl_dg_updater_diffusion_vlasov_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis, bool is_diff_const, const bool *diff_in_dir,
  int diff_order, const struct gkyl_range *diff_range, const bool *is_zero_flux_dir, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up diffusion updater object.
 * @param update_rng Range on which to compute.
 * @param coeff Diffusion coefficient/tensor.
 * @param fIn Input to updater.
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_diffusion_vlasov_advance(struct gkyl_dg_updater_diffusion_vlasov *up,
  const struct gkyl_range *update_rng, const struct gkyl_array *coeff,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in diffusion terms
 *
 * @param diffusion Updater object
 * @return timers
 */
struct gkyl_dg_updater_diffusion_vlasov_tm gkyl_dg_updater_diffusion_vlasov_get_tm(const struct gkyl_dg_updater_diffusion_vlasov *up);

/**
 * Delete updater.
 *
 * @param diffusion Updater to delete.
 */
void gkyl_dg_updater_diffusion_vlasov_release(struct gkyl_dg_updater_diffusion_vlasov *up);
