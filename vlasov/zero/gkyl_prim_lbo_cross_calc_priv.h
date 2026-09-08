#pragma once

#include <gkyl_prim_lbo_type.h>
#include <gkyl_rect_grid.h>

struct gkyl_prim_lbo_cross_calc {
  struct gkyl_rect_grid grid;
  const struct gkyl_prim_lbo_type *prim;

  bool is_first; // flag to indicate first call to update
  struct gkyl_nmat *As, *xs; // matrices for LHS and RHS
  gkyl_nmat_mem *mem; // memory for use in batched linear solve
  int nspecies; // number of species to cross-collide with

  uint32_t flags;
  struct gkyl_prim_lbo_cross_calc *on_dev; // pointer to itself or device data
};

#ifdef GKYL_HAVE_CUDA
/**
 * Create new updater to compute cross-primitive moments of 
 * distribution function on NV-GPU. See new() method for documentation.
 */
struct gkyl_prim_lbo_cross_calc *gkyl_prim_lbo_cross_calc_cu_dev_new(
  const struct gkyl_rect_grid *grid, struct gkyl_prim_lbo_type *prim);

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
void gkyl_prim_lbo_cross_calc_advance_cu(struct gkyl_prim_lbo_cross_calc *calc,
  const struct gkyl_range *conf_rng, const struct gkyl_array *alpha_E, double self_m,
  const struct gkyl_array *self_moms, const struct gkyl_array *self_prim_moms, double other_m,
  const struct gkyl_array *other_moms, const struct gkyl_array *other_prim_moms,
  const struct gkyl_array *boundary_corrections, const struct gkyl_array *nu,
  struct gkyl_array *prim_moms_out);
#endif
