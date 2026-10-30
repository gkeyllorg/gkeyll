#pragma once

#include <gkyl_array.h>
#include <gkyl_gyrokinetic.h>

/**
 * @file gkyl_gk_field_priv.h
 * @brief Private header for gyrokinetic field solver functions.
 *
 * This header contains forward declarations for internal field-related
 * functions used by the gyrokinetic application. These functions handle
 * FEM projections, field solves, boundary conditions, energy diagnostics,
 * FLR corrections and polarization potentials.
 */

/** FEM Initialization Functions **/

/**
 * Initialize the FEM object for the Boltzmann field solver.
 * Sets up the finite element machinery needed for Boltzmann electron response.
 *
 * @param app Gyrokinetic application object.
 * @param gkf Field object to initialize.
 */
void gk_field_fem_new_boltzmann(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf);

/**
 * Initialize the FEM object for 1D (1x) field solves.
 * Used for simulations with one configuration space dimension.
 *
 * @param app Gyrokinetic application object.
 * @param gkf Field object to initialize.
 */
void gk_field_fem_new_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf);

/**
 * Initialize the FEM object for 2D or 3D (2x/3x) field solves.
 * Used for simulations with two or three configuration space dimensions.
 *
 * @param app Gyrokinetic application object.
 * @param gkf Field object to initialize.
 */
void gk_field_fem_new_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf);

/** FEM Projection Functions **/

/**
 * Create the projection onto the continuous parallel (z) FEM basis, field->fem_parproj,
 * optionally weighted. With a weight it solves the 1x field equation.
 *
 * @param app Gyrokinetic application object.
 * @param gkf Field object.
 * @param weight Weight of the projection (NULL for none), defined on the local range.
 */
void gk_field_fem_parproj_new(
  struct gkyl_gyrokinetic_app *app, struct gk_field *gkf, const struct gkyl_array *weight
);

/**
 * Project a DG field onto the parallel FEM basis to enforce continuity
 * along the parallel (z) direction. This is essential for field solves
 * that require C0 continuity in the parallel direction.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing solver state.
 * @param arr_dg Input DG array to be projected.
 * @param arr_fem Output FEM array (projected result with parallel continuity).
 */
void gk_field_fem_projection_par(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *arr_dg,
  struct gkyl_array *arr_fem
);

/** Charge Density Accumulation Functions **/

/**
 * Accumulate charge density for a full Poisson solve.
 * Computes the charge density from all kinetic species for use
 * in solving the gyrokinetic Poisson equation.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object for storing charge density.
 * @param gks Species object providing density moments.
 * @param bflux Boundary flux array (for boundary contributions).
 */
void gk_field_accumulate_rho_c_poisson(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gk_species *gks,
  struct gkyl_array **bflux
);

/**
 * No-op function for parallel boundary conditions (disabled variant).
 * Used when parallel BCs are not needed or handled elsewhere.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object.
 * @param finout Array (unchanged).
 */
void gk_field_enforce_parallel_bc_disabled(
  const gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *finout
);

/** Adiabatic Species Functions **/

/**
 * Set up the response of the adiabatic species: allocate their background charge density and
 * K, and select the methods for the dimensionality (a weight in the 1x parallel projection, or
 * a flux-surface average and the Woodbury solver in 2x/3x). Call after the field solver is
 * initialized. K and the dimension-specific objects are computed by gk_field_adiabatic_init
 * once the initial conditions are applied.
 *
 * @param app Gyrokinetic application object.
 * @param gkf Field object.
 */
void gk_field_adiabatic_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf);

/**
 * Release the resources of the adiabatic species response.
 *
 * @param app Gyrokinetic application object.
 * @param gkf Field object.
 */
void gk_field_adiabatic_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf);

/** Finite Larmor Radius (FLR) Correction Functions **/

/**
 * Initialize FLR correction objects for the field solver.
 * Sets up machinery for gyroaveraging and FLR effects in the potential.
 *
 * @param app Gyrokinetic application object.
 * @param gkf Field object to initialize FLR corrections for.
 */
void gk_field_flr_new(struct gkyl_gyrokinetic_app *app, struct gk_field *gkf);

/**
 * Release FLR correction resources.
 *
 * @param app Gyrokinetic application object.
 * @param gkf Field object whose FLR resources are to be released.
 */
void gk_field_flr_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *gkf);

/**
 * Retrieve the electrostatic potential phi from the modified potential Phi_0 with the
 * simplified term phi = Phi_0 + (rho_i^2/eps_pol)*rho_c/J, valid for constant rho_i
 * and eps_pol. Also conserves Phi_0 for the field energy computation.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing solver state.
 * @param phi Electrostatic potential array (input/output).
 */
void gk_field_invert_flr_const(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi
);

/**
 * Same as gk_field_invert_flr_const but applying A = 1 - rho_i^2*nabla_perp^2
 * to Phi_0 with the FEM perpendicular operator.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing solver state.
 * @param phi Electrostatic potential array (input/output).
 */
void gk_field_invert_flr_op(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi
);

/**
 * No-op function for FLR inversion (when FLR corrections are disabled).
 * Used when long-wavelength approximation is valid and FLR effects are negligible.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object.
 * @param phi Electrostatic potential array (unchanged).
 */
void gk_field_invert_flr_none(
  gkyl_gyrokinetic_app *app, struct gk_field *field, struct gkyl_array *phi
);
