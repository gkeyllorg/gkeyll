#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type definition.
typedef struct gkyl_flash_interp_vars gkyl_flash_interp_vars;

/**
 * Initialize new updater for interpolating:
 * 1. Momentum density (rho u_x, rho u_y rho u_z) from finite volume to DG grid.
 * 2. Pressure tensor (p_ij) from DG grid to finite volume grid. 
 * 
 * @param basis DG basis.
 * @param grid_kinetic 
 * @param grid_fluid 
 * @param range_fluid_ext Extended fluid range.
 * @return struct gkyl_flash_interp_vars*
 */
struct gkyl_flash_interp_vars*
gkyl_flash_interp_vars_new(struct gkyl_basis *basis, 
  struct gkyl_rect_grid *grid_kinetic, struct gkyl_rect_grid *grid_fluid,
  struct gkyl_range *range_fluid_ext);

/**
 * Interpolate momentum density (rho u_x, rho u_y rho u_z) from kinetic DG grid
 * to coarse DG grid. Note that the function expects the input array to include
 * *all five* Euler variables, [rho, rho ux, rho uy, rho uz, p]
 * 
 * @param up 
 * @param range_kinetic 
 * @param range_fluid 
 * @param fv_fluid_in 
 * @param dg_rhou_out 
 */
void
gkyl_flash_interp_vars_rhou_to_dg(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *fv_fluid_in, struct gkyl_array *dg_rhou_out);

/**
 * Interpolate fluid velocity vector (u_x, u_y, u_z) from finite volume grid to kinetic DG grid.
 * 
 * @param up 
 * @param range_kinetic 
 * @param range_fluid 
 * @param fv_fluid_vel_in 
 * @param dg_vel_out 
 */
void
gkyl_flash_interp_vars_vel_to_dg(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *fv_fluid_vel_in, struct gkyl_array *dg_vel_out);

/**
 * Interpolate initial condition arrays [n, ux, uy, uz, T/m] onto DG kinetic grid.
 * 
 * @param up 
 * @param range_kinetic 
 * @param range_fluid 
 * @param fv_moms_in 
 * @param dg_moms_out 
 */
void
gkyl_flash_interp_vars_init_dist_moms_to_dg(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *fv_moms_in, struct gkyl_array *dg_moms_out);

/**
 * Interpolate [div(P), P : grad(u)] onto finite volume grid.
 * 
 * @param up 
 * @param range_kinetic 
 * @param range_fluid 
 * @param dg_coupling_fluxes 
 * @param fv_coupling_fluxes 
 */
void
gkyl_flash_interp_vars_coupling_fluxes_to_fv(struct gkyl_flash_interp_vars *up,
  struct gkyl_range *range_kinetic, struct gkyl_range *range_fluid,
  const struct gkyl_array *dg_coupling_fluxes, struct gkyl_array *fv_coupling_fluxes);

/**
 * Release memory allocated for updater. 
 * 
 * @param up 
 */
void
gkyl_flash_interp_vars_release(struct gkyl_flash_interp_vars *up);