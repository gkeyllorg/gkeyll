// Header containing functions required for coupling with FLASH
//
#pragma once

#include <gkyl_flash_multib.h>
#include <gkyl_flash.h>

/**
 * Used for coupling to an external fluid solver.
 * Store the moments [n, T/m] used to initialize the distribution function as a
 * Maxwellian. The input is finite volume data on the fluid grid (2x the
 * kinetic configuration-space resolution), interpolated onto the kinetic grid.
 * The species must have been created with init_dist_from_arrays = true.
 *
 * @param app App object.
 * @param init_dist_moms Finite volume moments [n, T/m].
 */
void
gkyl_flash_app_set_init_dist_moms_from_fv(gkyl_flash_app *app, struct gkyl_array *init_dist_moms);

/**
 * Used for coupling to an external fluid solver.
 * Initialize the distribution function as a Maxwellian matching finite volume
 * moments [n, T/m] on the fluid grid, then apply boundary conditions.
 *
 * @param app App object.
 * @param init_dist_moms Finite volume moments [n, T/m].
 */
void
gkyl_flash_app_init_maxwellian_from_moms_from_fv(gkyl_flash_app *app, struct gkyl_array *init_dist_moms);

/**
 * Pass array of fluid velocity [ux, uy, uz] to set the velocity of the moving
 * reference frame for the kinetic solver.
 * The input array is finite volume data living on a grid with 2x the configuration space
 * resolution of the kinetic grid in the app.
 * 
 * @param app 
 * @param fluid_mom_in 
 */
void
gkyl_flash_app_set_fluid_vel_from_fv(struct gkyl_flash_app *app, struct gkyl_array *fluid_vel_in);

/**
 * Used for coupling to an external fluid solver.
 * Set the fluid velocity that defines the moving frame for the transformed
 * kinetic equation on all local blocks, then sync it at inter-block
 * boundaries. fluid_vel_in[b] is finite volume data [ux, uy, uz] for local
 * block b (global block gkyl_flash_multib_app_local_block_id(app, b)).
 * Must be called once on all ranks.
 *
 * @param app Multiblock app object.
 * @param fluid_vel_in Fluid velocity arrays, one per local block.
 */
void
gkyl_flash_multib_app_set_fluid_vel_from_fv(gkyl_flash_multib_app *app, struct gkyl_array *fluid_vel_in[]);

/**
 * Used for coupling to an external fluid solver.
 * Compute the momentum and energy flux contributions for coupling 
 * the kinetic equation a fluid system, storing the result
 * in the coupling_fluxes array and interpolating to finite volume.
 * 
 * @param app 
 * @param coupling_fluxes 
 */
void
gkyl_flash_app_get_kinetic_fluxes_to_fv(struct gkyl_flash_app *app, struct gkyl_array *coupling_fluxes);

/**
 * Used for coupling to an external fluid solver.
 * Compute the momentum and energy flux contributions for coupling 
 * the kinetic equation a fluid system, storing the result
 * in the coupling_fluxes array and interpolating to finite volume.
 * 
 * @param app 
 * @param blockID 
 * @param coupling_fluxes 
 */
void
gkyl_flash_multib_app_get_kinetic_fluxes_to_fv(gkyl_flash_multib_app *app, int blockID, struct gkyl_array *coupling_fluxes);
