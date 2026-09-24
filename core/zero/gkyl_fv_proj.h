#pragma once

#include <gkyl_proj_on_basis.h>

// Object type
typedef struct gkyl_proj_on_basis gkyl_fv_proj;

/**
 * Create new updater to compute cell-average of function on a
 * grid. Free using gkyl_fv_proj_release method.
 *
 * @param grid Grid object
 * @param num_quad Number of quadrature nodes
 * @param num_ret_vals Number of values 'eval' sets
 * @param eval Function to project (See gkyl_proj_on_basis for signature).
 * @param ctx Context for function evaluation. Can be NULL.
 * @return New updater pointer.
 */
gkyl_fv_proj* gkyl_fv_proj_new(const struct gkyl_rect_grid *grid,
  int num_quad, int num_ret_vals, evalf_t eval, void *ctx);

/**
 * Compute cell averages. The update_rng MUST be a sub-range of
 * the range on which the array is defined. That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param pob Project on basis updater to run
 * @param tm Time at which projection must be computed
 * @param update_rng Range on which to run projection.
 * @param out Output array
 */
void gkyl_fv_proj_advance(const gkyl_fv_proj *pob,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_fv_proj_release(gkyl_fv_proj* pob);
