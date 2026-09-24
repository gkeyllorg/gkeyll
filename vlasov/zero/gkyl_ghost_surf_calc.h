#pragma once

#include <gkyl_util.h>
#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_ghost_surf_calc gkyl_ghost_surf_calc;

/**
 * Create new updater to update equations in the ghost cells using DG algorithm.
 *
 * @param grid_cu Grid object (on device)
 * @param equation Equation object
 * @param use_gpu bool to determine if on GPU
 */
struct gkyl_ghost_surf_calc* gkyl_ghost_surf_calc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_dg_eqn *equation, int cdim, bool use_gpu);

/**
 * Create new updater on CUDA device to update equations in the ghost cells using DG algorithm.
 *
 * @param grid_cu Grid object (on device)
 * @param equation Equation object
 */
struct gkyl_ghost_surf_calc* gkyl_ghost_surf_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_dg_eqn *equation, int cdim);

/**
 * Compute RHS of DG update in the ghost cells. The update_rng MUST be the external range
 * on which the array is defined or it will update the edge cells instead.
 *
 * @param gcalc Ghost surface updater object
 * @param phase_rng Phase space local_ext range.
 * @param conf_rng Configuration space range.
 * @param fIn Input to updater
 * @param rhs RHS output
 */
void gkyl_ghost_surf_calc_advance(gkyl_ghost_surf_calc *gcalc,
  const struct gkyl_range *phase_rng,
  const struct gkyl_array *fIn, struct gkyl_array *rhs);

// CUDA call
void gkyl_ghost_surf_calc_advance_cu(gkyl_ghost_surf_calc *gcalc,
  const struct gkyl_range *phase_rng,
  const struct gkyl_array *fIn, struct gkyl_array *rhs);
  
/**
 * Delete updater.
 *
 * @param gcalc Updater to delete.
 */
void gkyl_ghost_surf_calc_release(gkyl_ghost_surf_calc* gcalc);
