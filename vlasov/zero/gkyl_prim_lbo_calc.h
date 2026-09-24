#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_prim_lbo_type.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_prim_lbo_calc gkyl_prim_lbo_calc;

/**
 * Create new updater to compute primitive moments of distribution
 * function. Free using gkyl_prim_vlasov_calc_release.
 *
 * @param grid Grid object
 * @param prim Pointer to primitive moment type object
 * @param use_gpu bool to determine if on GPU
 * @return New updater pointer.
 */
struct gkyl_prim_lbo_calc* 
gkyl_prim_lbo_calc_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim, bool use_gpu);

/**
 * Compute primitive moments of distribution function. The phase_rng and conf_rng
 * MUST be a sub-ranges of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param calc Primitive moment calculator updater to run
 * @param conf_rng Config-space range
 * @param moms Moments of distribution function (Zeroth, First, and Second)
 * @param boundary_corrections Momentum and Energy boundary corrections
 * @param nu Collision frequency.
 * @param prim_moms_out Output drift velocity and thermal speed squared.
 */
void gkyl_prim_lbo_calc_advance(struct gkyl_prim_lbo_calc* calc, 
  const struct gkyl_range *conf_rng, const struct gkyl_array *moms,
  const struct gkyl_array *boundary_corrections, const struct gkyl_array *nu,
  struct gkyl_array *prim_moms_out);

/**
 * Delete pointer to primitive moment calculator updater.
 *
 * @param calc Updater to delete.
 */
void gkyl_prim_lbo_calc_release(struct gkyl_prim_lbo_calc* calc);

/**
 * Return pointer to primitive moment type structure.

 * @param calc Updater pointer.
 */
const struct gkyl_prim_lbo_type* gkyl_prim_lbo_calc_get_prim(gkyl_prim_lbo_calc* calc);

// "derived" class constructors
struct gkyl_prim_lbo_calc* 
gkyl_prim_lbo_vlasov_calc_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_rng, bool use_gpu);

struct gkyl_prim_lbo_calc* 
gkyl_prim_lbo_pkpm_calc_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_rng, bool use_gpu);

struct gkyl_prim_lbo_calc* 
gkyl_prim_lbo_gyrokinetic_calc_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_rng, bool use_gpu);
