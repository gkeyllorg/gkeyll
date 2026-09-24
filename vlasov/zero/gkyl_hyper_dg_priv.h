#pragma once

#include <gkyl_dg_eqn.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

struct gkyl_hyper_dg {
  struct gkyl_rect_grid grid; // grid object
  int ndim; // number of dimensions
  int num_basis; // number of basis functions
  int num_up_dirs; // number of update directions
  int update_dirs[GKYL_MAX_DIM]; // directions to update
  // zero_flux_flags[d] == 1 means zero-flux BC in 'd'
  int zero_flux_flags[2*GKYL_MAX_DIM];
  int update_vol_term; // should we update volume term?
  const struct gkyl_dg_eqn *equation; // equation object

  uint32_t flags;
  struct gkyl_hyper_dg *on_dev; // pointer to itself or device data
  bool use_gpu; // Whether to run on the gpu.
};

#ifdef GKYL_HAVE_CUDA

/**
 * Create new updater on CUDA device to update equations using DG algorithm.
 *
 * @param grid_cu Grid object (on device)
 * @param basis Basis functions
 * @param equation Equation object
 * @param num_up_dirs Number of directions to update
 * @param update_dirs List of directions to update (size 'num_up_dirs')
 * @param zero_flux_flags[2*GKYL_MAX_DIM] Flags to indicate if boundary has zero-flux BCs
 * @param update_vol_term Set to 0 to skip volume update
 */
gkyl_hyper_dg* gkyl_hyper_dg_cu_dev_new(const struct gkyl_rect_grid *grid_cu,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation_cu,
  int num_up_dirs, int update_dirs[GKYL_MAX_DIM], int zero_flux_flags[2*GKYL_MAX_DIM],
  int update_vol_term);

/**
 * Compute RHS of DG update on the device. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param hdg Hyper DG updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_hyper_dg_advance_cu(gkyl_hyper_dg* up, const struct gkyl_range *update_range,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs);

/**
 * Set if volume term should be computed or not.
 *
 * @param up Hyper DG updater object
 * @param update_vol_term Set to 1 to update vol term, 0 otherwise
 */
void
gkyl_hyper_dg_set_update_vol_cu(gkyl_hyper_dg *up, int update_vol_term);

#endif
