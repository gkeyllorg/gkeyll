#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_lbo_vlasov_drag.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_rad_vlasov gkyl_dg_updater_rad_vlasov;

// return type for drag and diffusion timers
struct gkyl_dg_updater_rad_vlasov_tm {
  double drag_tm; // time for drag and diffusion updates
};

/**
 * Create new updater to update radiation operator in Vlasov equation using hyper dg.
 *
 * @param phase_grid Phase space grid object
 * @param conf_basis Configuration space basis functions
 * @param phase_basis Phase space basis function
 * @param conf_range Configuration space range
 * @param drag_inp Input struct to vlasov drag operator (uses the gkyl_dg_lbo_vlasov_drag.h auxiliary struct) 
 * @param use_gpu Bool for whether updater is on host or device
 * @return New radiation updater object
 */
struct gkyl_dg_updater_rad_vlasov* 
gkyl_dg_updater_rad_vlasov_new(const struct gkyl_rect_grid *phase_grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_range, struct gkyl_dg_lbo_vlasov_drag_auxfields *drag_inp, 
  bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param rad Radiation updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_rad_vlasov_advance(struct gkyl_dg_updater_rad_vlasov *rad,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in drag terms
 *
 * @param rad Updater object
 * @return timers
 */
struct gkyl_dg_updater_rad_vlasov_tm gkyl_dg_updater_rad_vlasov_get_tm(const struct gkyl_dg_updater_rad_vlasov *rad);

/**
 * Delete updater.
 *
 * @param rad Updater to delete.
 */
void gkyl_dg_updater_rad_vlasov_release(struct gkyl_dg_updater_rad_vlasov *rad);
