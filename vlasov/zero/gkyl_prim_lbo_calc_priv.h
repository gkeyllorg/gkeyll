#pragma once

#include <gkyl_prim_lbo_type.h>
#include <gkyl_rect_grid.h>

struct gkyl_prim_lbo_calc {
  struct gkyl_rect_grid grid;
  const struct gkyl_prim_lbo_type *prim;

  bool is_first; // flag to indicate first call to update
  struct gkyl_nmat *As, *xs; // matrices for LHS and RHS
  gkyl_nmat_mem *mem; // memory for use in batched linear solve

  uint32_t flags;
  struct gkyl_prim_lbo_calc *on_dev; // pointer to itself or device data
};

#ifdef GKYL_HAVE_CUDA
/**
 * Create new updater to compute primitive moments of distribution function on
 * NV-GPU. See new() method for documentation.
 */
struct gkyl_prim_lbo_calc* 
gkyl_prim_lbo_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim);

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
 * @param prim_moms_out Output drift velocity and thermal speed squared.
 */
void gkyl_prim_lbo_calc_advance_cu(struct gkyl_prim_lbo_calc* calc, 
  const struct gkyl_range *conf_rng, const struct gkyl_array *moms,
  const struct gkyl_array *boundary_corrections, const struct gkyl_array *nu,
  struct gkyl_array* prim_moms_out);
#endif
