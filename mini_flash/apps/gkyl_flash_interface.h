// Header containing functions required for coupling with FLASH
//
#pragma once

#include <gkyl_flash_multib.h>
#include <gkyl_flash.h>

/**
 * Used for coupling to an external fluid solver.
 * Compute the momentum and energy flux contributions for coupling 
 * the kinetic equation a fluid system, storing the result
 * in the coupling_fluxes array.
 * 
 * @param app 
 * @param coupling_fluxes 
 */
void
gkyl_flash_app_get_kinetic_fluxes(struct gkyl_flash_app *app, struct gkyl_array *coupling_fluxes);

void
gkyl_flash_multib_app_get_kinetic_fluxes(gkyl_flash_multib_app *app, int blockID, struct gkyl_array *coupling_fluxes);


/**
 * Used for coupling to an external fluid solver.
 * Set the fluid velocity that defines the moving frame for the
 * transformed kinetic equation.
 * 
 * @param app 
 * @param fluid_mom_in 
 */
void
gkyl_flash_app_set_fluid_vel(struct gkyl_flash_app *app, struct gkyl_array *fluid_mom_in);

/**
 * @brief 
 * 
 * @param app 
 * @param blockID 
 * @param fluid_vel_in 
 */
void
gkyl_flash_multib_app_set_fluid_vel(gkyl_flash_multib_app *app, int blockID, struct gkyl_array *fluid_vel_in);

/**
 * Used for coupling to an external fluid solver.
 * Set the initial conditions for the distribution function as a Maxwellian
 * with moments given by the input array [n, ux, uy, uz, T/m].
 * 
 * @param app 
 * @param init_dist_moms 
 */
void
gkyl_flash_app_set_init_dist_moms(gkyl_flash_app *app, struct gkyl_array *init_dist_moms);

/**
 * Initialize the distribution function as a Maxwellian matching the moment array
 * of finite volume data [n, ux, uy, uz, T/m]. 
 * 
 * @param app 
 * @param init_dist_moms 
 */
void
gkyl_flash_app_init_maxwellian_from_moms(gkyl_flash_app *app, struct gkyl_array *init_dist_moms);

/**
 * Initialize the distribution function as a Maxwellian matching the moment array
 * of finite volume data [n, ux, uy, uz, T/m]. 
 * 
 * @param app 
 * @param init_dist_moms 
 */
void
gkyl_flash_multib_app_init_maxwellian_from_moms_block(gkyl_flash_multib_app *app, int blockID, struct gkyl_array *init_dist_moms);

/**
 * Update velocity space extents for each species' internal phase space and velocity space grids to 
 * the input lower and upper extents.
 * 
 * @param app 
 * @param lower_vel 
 * @param upper_vel 
 */
void 
gkyl_flash_app_set_vel_space_extents(gkyl_flash_app *app, double lower_vel[GKYL_MAX_VDIM], double upper_vel[GKYL_MAX_VDIM]);

/**
 * Update velocity space extents for each species' internal phase space and velocity space grids to 
 * the input lower and upper extents.
 * 
 * @param app 
 * @param lower_vel 
 * @param upper_vel 
 */
void 
gkyl_flash_multib_app_set_vel_space_extents(gkyl_flash_multib_app *app, int blockID, double lower_vel[GKYL_MAX_VDIM], double upper_vel[GKYL_MAX_VDIM]);


/**
 * Return the configuration space range over the fluid grid.
 * 
 * @param app 
 * @return struct gkyl_range 
 */
struct gkyl_range
gkyl_flash_app_get_fluid_conf_ext_range(gkyl_flash_app *app);

/**
 * Return the configuration space range over the fluid grid.
 * 
 * @param app 
 * @return struct gkyl_range 
 */
struct gkyl_range
gkyl_flash_app_get_fluid_conf_range(gkyl_flash_app *app);

struct gkyl_range
gkyl_flash_multib_app_get_fluid_conf_ext_range(gkyl_flash_multib_app *app, int blockID);

struct gkyl_range
gkyl_flash_multib_app_get_fluid_conf_range(gkyl_flash_multib_app *app, int blockID);
