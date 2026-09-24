#pragma once

#include <string.h>
#include <gkyl_moment_em_coupling.h>
#include <gkyl_sources_implicit_priv.h>
#include <gkyl_sources_explicit_priv.h>
#include <gkyl_gr_spacetime.h>

struct gkyl_moment_em_coupling {
  struct gkyl_rect_grid grid; // Grid over which the equations are solved.
  int ndim; // Number of grid dimensions.
  int nfluids; // Number of fluid species.
  struct gkyl_moment_em_coupling_data param[GKYL_MAX_SPECIES]; // Data for each fluid species.

  bool is_charged_species; // Does the fluid species carry a charge?
  double epsilon0; // Permittivity of free space.
  double mu0; // Permeability of free space.

  bool static_field; // Is the plasma field static? If true, only J is updated to new time step. 
  bool ramp_app_E; // Use a linear ramp function for initializing external electric fields.
  double t_ramp_E; // Ramp-up time for the linear ramp function for initializing external electric fields.
  bool ramp_app_curr; // Use a linear ramp function for initializing applied currents.
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
  double friction_Lambda_ee; // Electron-electron collisional terms for frictional sources.

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
