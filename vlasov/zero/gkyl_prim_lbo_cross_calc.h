#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_prim_lbo_cross_calc gkyl_prim_lbo_cross_calc;

/**
 * Create new updater to compute cross-primitive moments of 
 * distribution function. Free using gkyl_prim_vlasov_cross_calc_release.
 *
 * @param grid Grid object
 * @param prim Pointer to primitive moment type object
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_prim_lbo_cross_calc *gkyl_prim_lbo_cross_calc_new(
  const struct gkyl_rect_grid *grid, struct gkyl_prim_lbo_type *prim, bool use_gpu);

/**
 * Compute cross-primitive moments of distribution function. The conf_rng
 * MUST be a sub-range of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param calc Primitive moment calculator updater to run.
 * @param conf_rng Config-space range.
 * @param alpha_E Morse's alpha_E factor.
 * @param self_m Mass of the species.
 * @param self_moms Moments of distribution function (Zeroth, First, and Second).
 * @param self_prim_moms Drift velocity & thermal speed squared of this species.
 * @param other_m Mass of the colliding species.
 * @param other_moms Moments of distribution function (Zeroth, First, and Second).
 * @param other_prim_moms Drift velocity & thermal speed squared of the colliding species.
 * @param boundary_corrections Momentum and Energy boundary corrections.
 * @param nu Collision frequency.
 * @param prim_moms_out Output drift velocity and thermal speed squared.
 */
void gkyl_prim_lbo_cross_calc_advance(struct gkyl_prim_lbo_cross_calc *calc,
  const struct gkyl_range *conf_rng, const struct gkyl_array *alphaE, double self_m,
  const struct gkyl_array *self_moms, const struct gkyl_array *self_prim_moms, double other_m,
  const struct gkyl_array *other_moms, const struct gkyl_array *other_prim_moms,
  const struct gkyl_array *boundary_corrections, const struct gkyl_array *nu,
  struct gkyl_array *prim_moms_out);

/**
 * Delete pointer to primitive moment calculator updater.
 *
 * @param calc Updater to delete.
 */
void gkyl_prim_lbo_cross_calc_release(gkyl_prim_lbo_cross_calc *calc);

/**
 * Return pointer to primitive moment type structure.

 * @param calc Updater pointer.
 */
const struct gkyl_prim_lbo_type *gkyl_prim_lbo_cross_calc_get_prim(gkyl_prim_lbo_cross_calc *calc);

// "derived" class constructors
struct gkyl_prim_lbo_cross_calc *gkyl_prim_lbo_vlasov_cross_calc_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_rng, bool use_gpu);

struct gkyl_prim_lbo_cross_calc *gkyl_prim_lbo_gyrokinetic_cross_calc_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_rng, bool use_gpu);
