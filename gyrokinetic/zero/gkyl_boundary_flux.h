#pragma once

#include <gkyl_range.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_util.h>

// Object type
typedef struct gkyl_boundary_flux gkyl_boundary_flux;

/**
 * Create new updater to compute fluxes through the boundary surface due to the
 * surface terms in one or more equation systems (terms).
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param grid_cu Grid object (on device)
 * @param skin_r Skin range.
 * @param ghost_r Ghost range.
 * @param num_eqns Number of equation objects.
 * @param eqns Equation objects.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_boundary_flux *gkyl_boundary_flux_new(int dir, enum gkyl_edge_loc edge,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *skin_r,
  const struct gkyl_range *ghost_r, int num_equations, const struct gkyl_dg_eqn **eqns,
  bool use_gpu);

/**
 * Compute the boundary flux.
 *
 * @param up Boundary flux updater.
 * @param fIn Input distribution function.
 * @param rhs Output flux.
 */
void gkyl_boundary_flux_advance(
  gkyl_boundary_flux *up, const struct gkyl_array *fIn, struct gkyl_array *fluxOut);

/**
 * Free memory associated with this updater.
 *
 * @param up Updater to delete.
 */
void gkyl_boundary_flux_release(gkyl_boundary_flux *up);
