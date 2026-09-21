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
 * FLR corrections, polarization potentials, and biased wall configurations.
 */

/** FEM Initialization Functions **/

/**
 * Initialize the FEM object for the Boltzmann field solver.
 * Sets up the finite element machinery needed for Boltzmann electron response.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize.
 */
void
gk_field_fem_new_boltzmann(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize the FEM object for 1D (1x) field solves.
 * Used for simulations with one configuration space dimension.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize.
 */
void
gk_field_fem_new_1x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize the FEM object for 2D or 3D (2x/3x) field solves.
 * Used for simulations with two or three configuration space dimensions.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize.
 */
void
gk_field_fem_new_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/** FEM Projection Functions **/

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
void
gk_field_fem_projection_par(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *arr_dg, struct gkyl_array *arr_fem);

/** Charge Density Accumulation Functions **/

/**
 * Accumulate charge density for an adiabatic species response.
 * Computes the charge density contribution assuming an adiabatic
 * (flux-surface averaged) response for the species.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object for storing charge density.
 * @param s Species object providing density moments.
 * @param bflux Boundary flux array (for boundary contributions).
 */
void
gk_field_accumulate_rho_c_adiabatic(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gk_species *s, struct gkyl_array **bflux);

/**
 * Accumulate charge density for a full Poisson solve.
 * Computes the charge density from all kinetic species for use
 * in solving the gyrokinetic Poisson equation.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object for storing charge density.
 * @param s Species object providing density moments.
 * @param bflux Boundary flux array (for boundary contributions).
 */
void
gk_field_accumulate_rho_c_poisson(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gk_species *s, struct gkyl_array **bflux);

/**
 * No-op function for parallel boundary conditions (disabled variant).
 * Used when parallel BCs are not needed or handled elsewhere.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object.
 * @param finout Array (unchanged).
 */
void
gk_field_enforce_parallel_bc_disabled(const gkyl_gyrokinetic_app *app,
  struct gk_field *field, struct gkyl_array *finout);

/** Adiabatic Electron (2x/3x) Functions **/

/**
 * Set the background (electron) density times the Jacobian, n0*J, and the
 * background charge density q_e*n0*J, from electron_density or the optional
 * electron_density_profile. Used by all dimensionalities.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object.
 */
void
gk_field_adiabatic_density_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Set the coefficients of the adiabatic electron response, K = (e^2 n0/Te) J
 * and the Helmholtz coefficient kSq = -K. Call before creating the perpendicular
 * Helmholtz solver.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object.
 */
void
gk_field_adiabatic_coefs_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Initialize the flux-surface average, the 1D->cdim extension and the zonal
 * (Woodbury) system for the adiabatic electron response. Call after the
 * Helmholtz solver and the parallel smoothers exist.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object.
 */
void
gk_field_adiabatic_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Compute the flux-surface average <phi> = int J phi dy dz / int J dy dz
 * into f->adiab.psi (1D in x).
 *
 * @param app Gyrokinetic application object.
 * @param f Field object.
 * @param phi Potential to average.
 */
void
gk_field_adiabatic_fsa(gkyl_gyrokinetic_app *app, const struct gk_field *f, const struct gkyl_array *phi);

/**
 * Extend a 1D function of x to a cdim field constant along the other directions.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object.
 * @param psi 1D field in x.
 * @param out cdim output field.
 */
void
gk_field_adiabatic_inflate(gkyl_gyrokinetic_app *app, const struct gk_field *f, const struct gkyl_array *psi, struct gkyl_array *out);

/**
 * Solve the 2x/3x field equation with adiabatic electrons,
 * -div(eps grad_perp phi) + K (phi - <phi>) = rho, via the Woodbury identity.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object.
 */
void
gk_field_adiabatic_rhs_phi_2x3x(struct gkyl_gyrokinetic_app *app, struct gk_field *field);

/**
 * Add factor*(1/2) int K (phi-<phi>)^2 over the local range to out.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object.
 * @param factor Factor multiplying the integral.
 * @param out Scalar to accumulate into (device memory if use_gpu).
 */
void
gk_field_adiabatic_energy_accumulate(gkyl_gyrokinetic_app *app, const struct gk_field *f, double factor, double *out);

/**
 * Release the adiabatic electron response resources.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object.
 */
void
gk_field_adiabatic_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/** Finite Larmor Radius (FLR) Correction Functions **/

/**
 * Initialize FLR correction objects for the field solver.
 * Sets up machinery for gyroaveraging and FLR effects in the potential.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize FLR corrections for.
 */
void
gk_field_flr_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release FLR correction resources.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object whose FLR resources are to be released.
 */
void
gk_field_flr_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Apply FLR inversion to obtain the gyroaveraged potential.
 * Solves the FLR equation to compute the potential seen by particles.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing solver state.
 * @param phi Electrostatic potential array (input/output).
 */
void
gk_field_invert_flr(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *phi);

/**
 * No-op function for FLR inversion (when FLR corrections are disabled).
 * Used when long-wavelength approximation is valid and FLR effects are negligible.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object.
 * @param phi Electrostatic potential array (unchanged).
 */
void
gk_field_invert_flr_none(gkyl_gyrokinetic_app *app, struct gk_field *field,
  struct gkyl_array *phi);

/** Biased Wall Functions **/

/**
 * Initialize biased wall boundary condition objects.
 * Sets up machinery for applying time-dependent wall potentials.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object to initialize biased wall for.
 */
void
gk_field_biased_wall_new(struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Release biased wall resources.
 *
 * @param app Gyrokinetic application object.
 * @param f Field object whose biased wall resources are to be released.
 */
void
gk_field_biased_wall_release(const struct gkyl_gyrokinetic_app *app, struct gk_field *f);

/**
 * Calculate and apply the wall potential at the current simulation time.
 * Updates phi at the wall boundaries based on the biased wall configuration.
 *
 * @param app Gyrokinetic application object.
 * @param field Field object containing the potential.
 * @param tm Current simulation time.
 */
void
gk_field_calc_phi_wall(gkyl_gyrokinetic_app *app, struct gk_field *field, double tm);
