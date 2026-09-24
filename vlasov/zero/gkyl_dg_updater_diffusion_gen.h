#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_dg_updater_diffusion_gen gkyl_dg_updater_diffusion_gen;

// return type for drag and diffusion timers
struct gkyl_dg_updater_diffusion_gen_tm {
  double diffusion_tm; // time for diffusion updates
};

/**
 * Create new updater to update a fluid diffusion equation with
 * tensorial diffusion.
 *
 * @param grid Grid object.
 * @param basis Basis functions of the equation system.
 * @param diff_range Range object to index the diffusion coefficient.
 * @param use_gpu Whether to run on host or device.
 * @return New diff updater object
 */
struct gkyl_dg_updater_diffusion_gen* gkyl_dg_updater_diffusion_gen_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_range *diff_range, bool use_gpu);

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
void gkyl_dg_updater_diffusion_gen_advance(struct gkyl_dg_updater_diffusion_gen *up,
  const struct gkyl_range *update_rng, const struct gkyl_array *coeff,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Return total time spent in diffusion terms
 *
 * @param diffusion Updater object
 * @return timers
 */
struct gkyl_dg_updater_diffusion_gen_tm gkyl_dg_updater_diffusion_gen_get_tm(const struct gkyl_dg_updater_diffusion_gen *up);

/**
 * Delete updater.
 *
 * @param diffusion Updater to delete.
 */
void gkyl_dg_updater_diffusion_gen_release(struct gkyl_dg_updater_diffusion_gen *up);
