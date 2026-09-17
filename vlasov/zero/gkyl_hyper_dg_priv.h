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

GKYL_CU_DH
static int 
idx_to_inloup_ker(int dim, const int *idx, const int *dirs, const int *num_cells) {
  int iout = 0;

  for (int d=0; d<dim; ++d) {
    if (idx[dirs[d]] == 1) {
      iout = 2*iout+(int)(pow(3,d)+0.5);
    } else if (idx[dirs[d]] == num_cells[dirs[d]]) {
      iout = 2*iout+(int)(pow(3,d)+0.5)+1;
    }
  }
  return iout;
}

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
 * Compute RHS of DG generic stencil update on device.
 * In this case, generic stencil means all neighbor values are accessed and stored
 * in memory (i.e., in 2D, 9 cells are stored to update the center cell and in 3D
 * 27 cells are stored to update the center cell)
 * The update_rng MUST be a sub-range of the range on which the array is defined. 
 * That is, it must be either the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Hyper DG generic stencil updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_hyper_dg_gen_stencil_advance_cu(gkyl_hyper_dg* hdg, long offsets[36], const struct gkyl_range *update_rng,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, 
  struct gkyl_array *rhs);

/**
 * Set if volume term should be computed or not.
 *
 * @param up Hyper DG updater object
 * @param update_vol_term Set to 1 to update vol term, 0 otherwise
 */
void
gkyl_hyper_dg_set_update_vol_cu(gkyl_hyper_dg *up, int update_vol_term);

#endif
