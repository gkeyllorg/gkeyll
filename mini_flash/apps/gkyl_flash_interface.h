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
void gkyl_flash_app_set_init_dist_moms_from_fv(
  gkyl_flash_app *app, struct gkyl_array *init_dist_moms
);

/**
 * Used for coupling to an external fluid solver.
 * Initialize the distribution function as a Maxwellian matching finite volume
 * moments [n, T/m] on the fluid grid, then apply boundary conditions.
 *
 * @param app App object.
 * @param init_dist_moms Finite volume moments [n, T/m].
 */
void gkyl_flash_app_init_maxwellian_from_moms_from_fv(
  gkyl_flash_app *app, struct gkyl_array *init_dist_moms
);

/**
 * Used for coupling to an external fluid solver.
 * Initialize the distribution functions on all local blocks as Maxwellians
 * matching finite volume moments [n, T/m] on each block's fluid grid, then
 * apply physical and inter-block boundary conditions. init_dist_moms[b] holds
 * the moments for local block b (global block gkyl_flash_multib_app_local_block_id(app, b))
 * on its extended fluid range (gkyl_flash_multib_app_get_fluid_conf_ext_range),
 * including ghost cells. Must be called once on all ranks.
 *
 * @param app Multiblock app object.
 * @param init_dist_moms Finite volume moments [n, T/m], one per local block.
 */
void gkyl_flash_multib_app_init_maxwellian_from_moms_from_fv(gkyl_flash_multib_app *app);

/**
 * Used for coupling to an external fluid solver.
 * Set the array of moments used to initialize a Maxwellian for the given block with
 * local block ID on the caller's rank, `blockId_local`.
 * 
 * @param app 
 * @param blockId_local 
 * @param init_dist_moms 
 */
void gkyl_flash_multib_app_set_init_dist_moms_from_fv_block(
  gkyl_flash_multib_app *app, int blockId_local, struct gkyl_array *init_dist_moms
);

/**
 * Pass array of fluid velocity [ux, uy, uz] to set the velocity of the moving
 * reference frame for the kinetic solver.
 * The input array is finite volume data living on a grid with 2x the configuration space
 * resolution of the kinetic grid in the app.
 * 
 * @param app 
 * @param fluid_mom_in 
 */
void gkyl_flash_app_set_fluid_vel_from_fv(
  struct gkyl_flash_app *app, struct gkyl_array *fluid_vel_in
);

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
void gkyl_flash_multib_app_set_fluid_vel_from_fv(
  gkyl_flash_multib_app *app, struct gkyl_array *fluid_vel_in[]
);

/**
 * Used for coupling to an external fluid solver.
 * Set the fluid velocity [ux, uy, uz] (finite volume, on the extended fluid
 * range) of local block blockId_local, without communication. After setting
 * all local blocks, call gkyl_flash_multib_app_sync_fluid_vel once on all ranks
 * before the next update.
 *
 * @param app Multiblock app object.
 * @param blockId_local Local block index.
 * @param fluid_vel_in Finite volume fluid velocity for this block.
 */
void gkyl_flash_multib_app_set_fluid_vel_from_fv_block(
  gkyl_flash_multib_app *app, int blockId_local, struct gkyl_array *fluid_vel_in
);

/**
 * Used for coupling to an external fluid solver.
 * Sync the fluid velocity set with gkyl_flash_multib_app_set_fluid_vel_from_fv_block
 * at inter-block boundaries. Must be called once on all ranks, after all local
 * blocks have been set.
 *
 * @param app Multiblock app object.
 */
void gkyl_flash_multib_app_sync_fluid_vel(gkyl_flash_multib_app *app);

/**
 * Used for coupling to an external fluid solver.
 * Compute the maximum stable kinetic timestep of local block blockId_local,
 * without communication. The timestep passed to
 * gkyl_flash_multib_app_kinetic_update_block must be the same on all blocks
 * and ranks, so the caller must reduce these (minimum over all local blocks
 * and all ranks) first.
 *
 * @param app Multiblock app object.
 * @param blockId_local Local block index.
 * @return Maximum stable timestep of this block.
 */
double gkyl_flash_multib_app_compute_dt_max_block(gkyl_flash_multib_app *app, int blockId_local);

/**
 * Used for coupling to an external fluid solver.
 * Take a single forward Euler step of ONLY the kinetic equation on local
 * block blockId_local with timestep dt, and apply its physical boundary
 * conditions, without communication. dt must not exceed the global stable
 * timestep (gkyl_flash_multib_app_compute_dt_max); the app exits if this
 * block can only take a smaller step, which would leave the blocks out of
 * sync. After updating all local blocks, call
 * gkyl_flash_multib_app_kinetic_update_sync once on all ranks.
 *
 * @param app Multiblock app object.
 * @param blockId_local Local block index.
 * @param dt Timestep.
 */
void gkyl_flash_multib_app_kinetic_update_block(
  gkyl_flash_multib_app *app, int blockId_local, double dt
);

/**
 * Used for coupling to an external fluid solver.
 * Complete a kinetic step taken with gkyl_flash_multib_app_kinetic_update_block:
 * advance the app time by dt and sync the distribution functions at
 * inter-block boundaries. Must be called once on all ranks after all local
 * blocks have been updated.
 *
 * @param app Multiblock app object.
 * @param dt Timestep taken by the blocks.
 */
void gkyl_flash_multib_app_kinetic_update_sync(gkyl_flash_multib_app *app, double dt);

/**
 * Used for coupling to an external fluid solver.
 * Compute the momentum and energy flux contributions for coupling 
 * the kinetic equation a fluid system, storing the result
 * in the coupling_fluxes array and interpolating to finite volume.
 * 
 * @param app 
 * @param coupling_fluxes 
 */
void gkyl_flash_app_get_kinetic_fluxes_to_fv(
  struct gkyl_flash_app *app, struct gkyl_array *coupling_fluxes
);

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
void gkyl_flash_multib_app_get_kinetic_fluxes_to_fv_block(
  gkyl_flash_multib_app *app, int blockId_local, struct gkyl_array *coupling_fluxes
);

// Initialization from an external caller.
// Self-contained helpers meant to be copied into the external code.

// Maximum number of blocks on a single rank. The external code may define it
// (e.g. as FLASH's MAXBLOCKS) before including this header.
#ifndef GKYL_FLASH_MAX_LOCAL_BLOCKS
#define GKYL_FLASH_MAX_LOCAL_BLOCKS 1000
#endif

// Context struct for creating a multi-block app from an external caller.
// Arrays indexed by b hold data for local block b (global ID local_block_ids[b]).
// struct gkyl_flash_coupling_ctx {
//   char name[128]; // App name (defaults to "gkyl_app" if empty).
//   int cdim, vdim; // Configuration and velocity space dimensionality.

//   double mass; // Species mass.
//   double nu; // Constant self-collision frequency.
//   double epsilon0, mu0; // Vacuum permittivity and permeability.
//   double B0[3]; // Uniform static magnetic field [Bx, By, Bz].

//   int num_blocks; // Number of blocks in the global domain.
//   int num_local_blocks; // Number of blocks on the current rank.
//   int local_block_ids[GKYL_FLASH_MAX_LOCAL_BLOCKS]; // 0-based global IDs of local blocks, ascending.

//   int cells_block[GKYL_MAX_CDIM]; // Fluid (finite volume) cells per block, identical for all blocks.
//   // Configuration space bounds of each local block.
//   double lower_block[GKYL_FLASH_MAX_LOCAL_BLOCKS][GKYL_MAX_CDIM];
//   double upper_block[GKYL_FLASH_MAX_LOCAL_BLOCKS][GKYL_MAX_CDIM];
//   // Global IDs of the neighbours of each local block at the (lower, upper) edge
//   // in each direction, or -1 at a physical boundary. Periodic directions are
//   // given by the wrapped neighbour IDs.
//   int neigh_block_ids[GKYL_FLASH_MAX_LOCAL_BLOCKS][GKYL_MAX_CDIM][2];

//   // Species BCs at the (lower, upper) physical domain boundary in each
//   // direction. Only used at edges with neighbour ID -1.
//   enum gkyl_species_bc_type bcs[GKYL_MAX_CDIM][2];

//   int cells_vel[3]; // Velocity space cells.
//   double lower_vel[3], upper_vel[3]; // Velocity space bounds.

//   double cfl_frac; // CFL scaling factor for the kinetic time step (defaults to 1.0).
//   bool write_pkpm_diags; // Flag to write PKPM diagnostic moments.
// };

// /**
//  * Create a multi-block app with a single kinetic species and no fluid
//  * species from a coupling context. The distribution function is initialized
//  * later from moment arrays (e.g. gkyl_flash_multib_app_init_maxwellian_from_moms),
//  * and the static field holds the uniform magnetic field B0.
//  * Uses the linear decomposition (block b on rank b), so the local block IDs
//  * in ctx must match the ones the app assigns; the app exits otherwise.
//  * The app keeps its own references to comm, which the caller may release
//  * at any time. Must be called on all ranks.
//  *
//  * @param ctx Coupling context.
//  * @param comm Communicator over all ranks.
//  * @return New multi-block app.
//  */
// gkyl_flash_multib_app*
// gkyl_flash_multib_app_new_from_ctx(const struct gkyl_flash_coupling_ctx *ctx, struct gkyl_comm *comm);
