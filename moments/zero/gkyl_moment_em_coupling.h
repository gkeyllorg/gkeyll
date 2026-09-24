#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_eqn_type.h>
#include <gkyl_fv_proj.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_gr_spacetime.h>

struct gkyl_moment_em_coupling_data {
  enum gkyl_eqn_type type; // Equation type.

  double charge; // Species charge.
  double mass; // Species mass.

  double k0; // Closure parameter (for 10-moment equations only; defaults to 0.0).
};

struct gkyl_moment_em_coupling_inp {
  const struct gkyl_rect_grid *grid; // Grid over which the equations are solved.
  int nfluids; // Number of fluid species.
  struct gkyl_moment_em_coupling_data param[GKYL_MAX_SPECIES]; // Data for each fluid species.
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.

  bool static_field; // Is the plasma field static? If true, only J is updated to new time step. 
  double t_ramp_E; // Ramp-up time for the linear ramp function for initializing external electric fields.
  double t_ramp_curr; // Ramp-up time for the linear ramp function for initializing applied currents.

  bool has_collision; // Run with collisions switched on.
  bool use_rel; // Assume special relativistic fluid species.
  
  // Matrix of scaling factors for collision frequencies. Should be symmetric (i.e. nu_base_sr = nu_base_rs).
  // These are defined such that nu_sr = nu_base_sr / rho_s, and nu_rs = nu_base_rs / rho_r.
  double nu_base[GKYL_MAX_SPECIES][GKYL_MAX_SPECIES];

  bool use_explicit_em_coupling; // Use the explicit source-solver for handling moment-EM coupling (not operational yet).

  bool has_nT_sources; // Run with number density and temperature sources.

  bool has_frictional_sources; // Run with frictional sources.
  bool use_explicit_friction; // Use an explicit (SSP-RK3) solver for integrating frictional sources.
  double friction_Z; // Ionization number for frictional sources.
  double friction_T_elc; // Electron temperature for frictional sources.
  double friction_Lambda_ee; // Electron-electron collisional term for frictional sources.

  bool has_volume_sources; // Run with volume-based geometrical sources.
  double volume_gas_gamma; // Adiabatic index for volume-based geometrical sources.
  double volume_U0; // Initial comoving plasma velocity for volume-based geometrical sources.
  double volume_R0; // Initial radial distance from expansion/contraction center for volume-based geometrical sources.

  bool has_reactive_sources; // Run with reactive sources.
  double reactivity_gas_gamma; // Adiabatic index for reactive sources.
  double reactivity_specific_heat_capacity; // Specific heat capacity for reactive sources.
  double reactivity_energy_of_formation; // Energy of formation for reactive sources.
  double reactivity_ignition_temperature; // Ignition temperature for reactive sources.
  double reactivity_reaction_rate; // Reaction rate for reactive sources.

  bool has_einstein_medium_sources; // Run with coupled fluid-Einstein sources in plane-symmetric spacetimes.
  double medium_gas_gamma; // Adiabatic index for coupled fluid-Einstein sources in plane-symmetric spacetimes.
  double medium_kappa; // Stress-energy prefactor for coupled fluid-Einstein sources in plane-symmetric spacetimes.

  bool has_gr_ultra_rel_sources; // Run with general relativistic source terms (Euler equations, ultra-relativistic equation of state).
  double gr_ultra_rel_gas_gamma; // Adiabatic index for general relativistic Euler equations (ultra-relativistic equation of state).

  bool has_gr_euler_sources; // Run with general relativistic source terms (Euler equations, general equation of state).
  double gr_euler_gas_gamma; // Adiabatic index for general relativistic Euler equations (general equation of state).

  bool has_gr_twofluid_sources; // Run with general relativistic two-fluid source terms.
  double gr_twofluid_mass_elc; // Electron mass for general relativistic two-fluid equations.
  double gr_twofluid_mass_ion; // Ion mass for general relativistic two-fluid equations.
  double gr_twofluid_charge_elc; // Electron charge for general relativistic two-fluid equations.
  double gr_twofluid_charge_ion; // Ion charge for general relativistic two-fluid equations.
  double gr_twofluid_gas_gamma_elc; // Adiabatic index for electrons in general relativistic two-fluid equations.
  double gr_twofluid_gas_gamma_ion; // Adiabatic index for ions in general relativistic two-fluid equations.
  double gr_twofluid_e_fact; // Electric field divergence error propagation speed for general relativistic two-fluid equations.

  bool has_vacuum_einstein_sources; // Run with vacuum Einstein sources using the Bona-Masso formalism.
  double vacuum_einstein_excision_threshold; // Excision threshold (lapse) for vacuum Einstein equations using the Bona-Masso formalism.
  enum gkyl_spacetime_slicing vacuum_einstein_spacetime_slicing; // Spacetime slicing condition for vacuum Einstein equations using the Bona-Masso formalism.
  enum gkyl_spacetime_evolution vacuum_einstein_spacetime_evolution; // Spacetime evolution system for vacuum Einstein equations using the Bona-Masso formalism.

  bool has_vacuum_einstein_conformal_sources; // Run with vacuum Einstein sources using the conformal Bona-Masso formalism.
  double vacuum_einstein_conformal_excision_threshold; // Excision threshold (lapse) for vacuum Einstein equations using the conformal Bona-Masso formalism.
  enum gkyl_spacetime_slicing vacuum_einstein_conformal_spacetime_slicing; // Spacetime slicing condition for vacuum Einstein equations using the conformal Bona-Masso formalism.
  enum gkyl_spacetime_evolution vacuum_einstein_conformal_spacetime_evolution; // Spacetime evolution system for vacuum Einstein equations using the conformal Bona-Masso formalism.
  
  bool has_gr_mhd_sources; // Run with general relativistic source terms (general relativistic magnetohydrodynamics equations).
  double gr_mhd_gas_gamma; // Adiabatic index for general relativistic magnetohydrodynamics equations.
};

// Moment-EM coupling object.
typedef struct gkyl_moment_em_coupling gkyl_moment_em_coupling;

/**
* Create a new moment-EM coupling object, used for integrating the electromagnetic source terms that appear in the multi-fluid equations.
*
* @param inp Input parameters for the moment-EM coupling object.
* @return Moment-EM coupling object.
*/
gkyl_moment_em_coupling*
gkyl_moment_em_coupling_new(struct gkyl_moment_em_coupling_inp inp);

/**
* Integrate the electromagnetic source terms in the multi-fluid equation system using an implicit forcing solver (specifically the time-centered
* Crank-Nicolson/implicit Runge-Kutta method). The gkyl_range object to be updated should be a (non-strict) subrange of the range over which the
* fluid variable array is defined (i.e. either the fluid variable array gkyl_range object itself, or a gkyl_range object initialized using the
* gkyl_sub_range_init method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param update_range Range object over which to integrate the electromagnetic sources using an implicit time-centered method.
* @param fluid Array of fluid variables (array size = nfluids).
* @param app_accel Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param p_rhs Array of RHS/source terms to be applied to the fluid variables  
*              (e.g., Braginskii transport for Euler/Isothermal Euler; gradient-based closure for ten-moment).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param nT_sources Array of number density and temperature source terms.
*/
void
gkyl_moment_em_coupling_implicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array *p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array* app_current, const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES]);

/**
* Integrate the electromagnetic source terms in the multi-fluid equation system using an explicit forcing solver (specifically either the strong
* stability-preserving third-order Runge-Kutta method or the simple first-order forward-Euler method). The gkyl_range object to be updated should be a
* (non-strict) subrange of the range over which the fluid variable array is defined (i.e. either the fluid variable array gkyl_range object itself, or
* a gkyl_range object initialized using the gkyl_sub_range_init method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param update_range Range object over which to integrate the electromagnetic sources using an explicit time integration method.
* @param fluid Array of fluid variables (array size = nfluids).
* @param app_accel Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param p_rhs Array of RHS/source terms to be applied to the pressure tensor (for the case of 10-moment gradient-based closure only).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param app_current1 Array of stage-1 current terms to be applied to the fluid equations (for stage-1 of external current driving).
* @param app_current2 Array of stage-2 current terms to be applied to the fluid equations (for stage-2 of external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param nT_sources Array of number density and temperature source terms.
* @param proj_app_curr The finite-volume projection routine for the external current.
* @param nstrang Indicator of which step in the Strang splitting we are currently considering.
*/
void
gkyl_moment_em_coupling_explicit_advance(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, const struct gkyl_range* update_range,
  struct gkyl_array* fluid[GKYL_MAX_SPECIES], const struct gkyl_array* app_accel[GKYL_MAX_SPECIES], const struct gkyl_array* p_rhs[GKYL_MAX_SPECIES],
  struct gkyl_array* em, const struct gkyl_array *app_current, const struct gkyl_array* app_current1, const struct gkyl_array* app_current2,
  const struct gkyl_array* ext_em, const struct gkyl_array* nT_sources[GKYL_MAX_SPECIES], gkyl_fv_proj* proj_app_curr, int nstrang);

/**
* Delete moment-EM coupling object.
*
* @param mom_em Moment-EM coupling object to delete.
*/
void
gkyl_moment_em_coupling_release(gkyl_moment_em_coupling* mom_em);