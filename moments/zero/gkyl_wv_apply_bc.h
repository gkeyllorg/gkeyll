#pragma once

#include <gkyl_array.h>
#include <gkyl_evalf_def.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>

// Object type
typedef struct gkyl_wv_apply_bc gkyl_wv_apply_bc;

/**
 * Create new updater to apply set of boundary conditions.
 *
 * @param grid Grid object
 * @param eqn Equation object to use
 * @param geom Wave geometry to use in boundary conditions
 * @param dir Direction to apply BC 
 * @param edge Edge to apply BC (0: left-edge; 1: right-edge)
 * @param nghost Number of ghost cells in each direction
 * @param Boundary condition function to apply
 * @param ctx Context to pass to bcfunc.
 * @return New updater pointer.
 */
gkyl_wv_apply_bc* gkyl_wv_apply_bc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_wv_eqn *eqn, const struct gkyl_wave_geom *geom,
  int dir, enum gkyl_edge_loc edge, const int *nghost,
  wv_bc_func_t bcfunc, void *ctx);

/**
 * Apply boundary condition on specified field. If the update_rng does
 * not touch the edge in specified dir then nothing is done.
 *
 * @param bc BC updater to run
 * @param tm Time at which BC is applied
 * @param update_rng Range on which BC is applied. See note above.
 * @param out Output array
 */
void gkyl_wv_apply_bc_advance(const gkyl_wv_apply_bc *bc, double tm,
  const struct gkyl_range *update_rng, struct gkyl_array *out);

/**
 * Apply boundary condition on specified field, however, copying the
 * output to a buffer instead of the field itself. If the update_rng
 * does not touch the edge in specified dir then nothing is done.
 *
 * @param bc BC updater to run
 * @param tm Time at which BC is applied
 * @param update_rng Range on which BC is applied. See note above.
 * @param inp Input array
 * @param buffer Output buffer in which BCs are copied
 */
void gkyl_wv_apply_bc_to_buff(const gkyl_wv_apply_bc *bc, double tm,
  const struct gkyl_range *update_rng, const struct gkyl_array *inp, double *buffer);

/**
 * Delete updater.
 *
 * @param bc Updater to delete.
 */
void gkyl_wv_apply_bc_release(gkyl_wv_apply_bc* bc);
