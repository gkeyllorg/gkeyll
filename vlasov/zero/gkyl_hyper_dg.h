#pragma once

#include <gkyl_util.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_hyper_dg gkyl_hyper_dg;

/**
 * Create new updater to update equations using DG algorithm.
 *
 * @param grid Grid object
 * @param basis Basis functions
 * @param equation Equation object
 * @param num_up_dirs Number of directions to update
 * @param update_dirs List of directions to update (size 'num_up_dirs')
 * @param zero_flux_flags[2*GKYL_MAX_DIM] Flags to indicate if boundary has zero-flux BCs
 * @param update_vol_term Set to 0 to skip volume update
 * @param use_gpu bool to determine if on GPU
 */
gkyl_hyper_dg *gkyl_hyper_dg_new(
  const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  const struct gkyl_dg_eqn *equation, int num_up_dirs, int update_dirs[GKYL_MAX_DIM],
  int zero_flux_flags[2 * GKYL_MAX_DIM], int update_vol_term, bool use_gpu
);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Hyper DG updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_hyper_dg_advance(
  gkyl_hyper_dg *up, const struct gkyl_range *update_rng, const struct gkyl_array *fIn,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
);

// Number of relative offsets used by the generic stencil update: a 3-cell stencil for
// each of (up to) three update directions and a 9-cell stencil for each of (up to)
// three pairs of update directions.
#define GKYL_HYPER_DG_GEN_NUM_OFFSETS 36

/**
 * Construct the relative (linear index) offsets of the stencils used by the generic
 * stencil update. For update directions d1 <= d2 (indices into update_dirs):
 *   - offsets[3*d1 + i], i=0,1,2: lower, center and upper cell along d1 (diagonal term).
 *   - offsets[9*(d1+d2) + 3*i + j], i,j=0,1,2: cell at position i along d1 and
 *     position j along d2 (cross terms d1,d2 and d2,d1).
 *
 * @param num_up_dirs Number of update directions (at most 3).
 * @param update_dirs Update directions, in increasing order.
 * @param range Range the updated arrays are defined on (or a sub-range of it).
 * @param offsets On output, relative offsets of the stencils.
 */
void gkyl_hyper_dg_gen_stencil_offsets(
  int num_up_dirs, const int *update_dirs, const struct gkyl_range *range,
  long offsets[GKYL_HYPER_DG_GEN_NUM_OFFSETS]
);

/**
 * Compute RHS of DG generic stencil update, for second-order operators with mixed
 * partial derivatives (e.g. a full diffusion tensor). The term for each pair of update
 * directions (d1,d2) is updated with its own stencil: 3 cells along d1 if d1=d2, and
 * the 9 cells in the (d1,d2) plane otherwise. The equation's gen_surf_term is passed a
 * kernel index that encodes whether the cell is in the interior or at a lower/upper
 * edge of update_rng in the stencil's directions (assuming update_rng spans the whole
 * domain in the update directions, starting at index 1), so equations with special
 * kernels at domain boundaries can select them; equations that rely on ghost cells
 * instead can ignore it.
 * The update_rng MUST be a sub-range of the range on which the array is defined.
 * That is, it must be either the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Hyper DG generic stencil updater object
 * @param offsets Relative offsets of the stencils (see gkyl_hyper_dg_gen_stencil_offsets).
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_hyper_dg_gen_stencil_advance(
  gkyl_hyper_dg *up, long offsets[GKYL_HYPER_DG_GEN_NUM_OFFSETS],
  const struct gkyl_range *update_rng, const struct gkyl_array *fIn, struct gkyl_array *cflrate,
  struct gkyl_array *rhs
);

/**
 * Set if volume term should be computed or not.
 *
 * @param up Hyper DG updater object
 * @param update_vol_term Set to 1 to update vol term, 0 otherwise
 */
void gkyl_hyper_dg_set_update_vol(gkyl_hyper_dg *up, int update_vol_term);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_hyper_dg_release(gkyl_hyper_dg *up);
