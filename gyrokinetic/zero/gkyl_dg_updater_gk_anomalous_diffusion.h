#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_gk_bc_type.h>

// Object type
typedef struct gkyl_dg_updater_gk_anomalous_diffusion gkyl_dg_updater_gk_anomalous_diffusion;

// return type for drag and diffusion timers
struct gkyl_dg_updater_gk_anomalous_diffusion_tm {
  double diffusion_tm; // time for diffusion updates
};

/**
 * Create new updater to apply the anomalous gyrokinetic diffusion operator using hyper dg.
 *
 * @param grid Grid object.
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param conf_range Conf-space range object.
 * @param bc_x_lower Boundary condition at lower x boundary.
 * @param bc_x_upper Boundary condition at upper x boundary.
 * @param nu Diffusivity.
 * @param jacobgeo_inv Reciprocal of the configuration space jacobian.
 * @param use_gpu Whether to run on host or device.
 * @return New diff updater object
 */
struct gkyl_dg_updater_gk_anomalous_diffusion* gkyl_dg_updater_gk_anomalous_diffusion_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range,
  enum gkyl_gyrokinetic_bc_type bc_x_lower, enum gkyl_gyrokinetic_bc_type bc_x_upper,
  const struct gkyl_array *nu, const struct gkyl_array *jacobgeo_inv, bool use_gpu);

/**
 * Acquire anomalous gyrokinetic diffusion equation object.
 *
 * @param up anomalous gyrokinetic diffusion updater object.
 * 
 * @return anomalous gyrokinetic diffusion equation object.
 */
struct gkyl_dg_eqn* 
gkyl_dg_updater_gk_anomalous_diffusion_acquire_eqn(const gkyl_dg_updater_gk_anomalous_diffusion* up);

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
void gkyl_dg_updater_gk_anomalous_diffusion_advance(struct gkyl_dg_updater_gk_anomalous_diffusion *up,
  const struct gkyl_range *update_rng, const struct gkyl_array* GKYL_RESTRICT fIn,
  struct gkyl_array* GKYL_RESTRICT cflrate, struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in anomalous diffusion terms.
 *
 * @param diffusion Updater object
 * @return timers
 */
struct gkyl_dg_updater_gk_anomalous_diffusion_tm gkyl_dg_updater_gk_anomalous_diffusion_get_tm(const struct gkyl_dg_updater_gk_anomalous_diffusion *up);

/**
 * Delete updater.
 *
 * @param diffusion Updater to delete.
 */
void gkyl_dg_updater_gk_anomalous_diffusion_release(struct gkyl_dg_updater_gk_anomalous_diffusion *up);
