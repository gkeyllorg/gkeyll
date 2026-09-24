#pragma once
#include <gkyl_gr_spacetime.h>

// Forward-declaration of the private gkyl_moment_em_coupling object type.
typedef struct gkyl_moment_em_coupling gkyl_moment_em_coupling;

/**
* Integrate the number density and temperature source terms for a single fluid species within a single cell, using an explicit forcing solver
* (specifically a simple first-order forward-Euler method).
*
* @param mass Mass of the fluid species.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables for a single fluid species (before source update).
* @param fluid_new Array of new fluid variables for a single fluid species (after source update).
* @param nT_sources Array of number density and temperature source terms for a single fluid species.
*/
void
explicit_nT_source_update_euler(const double mass, const double dt, double* fluid_old, double* fluid_new, const double* nT_sources);

/**
* Integrate the number density and temperature source terms in the multi-fluid equation system within a single cell, using an explicit forcing
* solver (specifically a simple first-order forward-Euler method)
*
* @param mom_em Moment-EM coupling object.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param nT_sources_s Array of number density and temperature source terms.
*/
void
explicit_nT_source_update(const gkyl_moment_em_coupling* mom_em, const double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* nT_sources_s[GKYL_MAX_SPECIES]);

/**
* Integrate the frictional source terms in the multi-fluid equation system within a single cell, using an explicit forcing solver (specifically a
* simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param Z Ionization number.
* @param T_elc Electron temperature.
* @param Lambda_ee Electron-electron collisional term.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param f_elc_old Array of old electron fluid variables (before source update).
* @param f_ion_old Array of old ion fluid variables (before source update).
* @param f_elc_new Array of new electron fluid variables (after source update).
* @param f_ion_new Array of new ion fluid variables (after source update).
*/
void
explicit_frictional_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double Z, const double T_elc, const double Lambda_ee,
  double t_curr, const double dt, double* f_elc_old, double* f_ion_old, double* f_elc_new, double* f_ion_new);

/**
* Integrate the frictional source terms in the multi-fluid equation system within a single cell, using an explicit forcing solver (specifically a
* strong stability-preserving third-order Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
*/
void
explicit_frictional_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES]);

/**
* Integrate the volume-based geometrical source terms (e.g. for expanding/contracting box formalism) in the multi-fluid equation system for a
* single 5-moment fluid, within a single cell, using an explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param gas_gamma Adiabatic index.
* @param U0 (Initial) comoving plasma velocity.
* @param R0 (Initial) radial distance from expansion/contraction center.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old single-species 5-moment fluid variables (before source update).
* @param fluid_new Array of new single-species 5-moment fluid variables (after source update).
*/
void
explicit_volume_source_5m_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, const double U0, const double R0,
  double t_curr, const double dt, double* fluid_old, double* fluid_new);

/**
* Integrate the volume-based geometrical source terms (e.g. for expanding/contracting box formalism) in the multi-fluid equation system for a
* single 10-moment fluid, within a single cell, using an explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param U0 (Initial) comoving plasma velocity.
* @param R0 (Initial) radial distance from expansion/contraction center.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old single-species 10-moment fluid variables (before source update).
* @param fluid_new Array of new single-species 10-moment fluid variables (after source update).
*/
void
explicit_volume_source_10m_update_euler(const gkyl_moment_em_coupling* mom_em, const double U0, const double R0, double t_curr, const double dt,
  double* fluid_old, double* fluid_new);

/**
* Integrate the volume-based geometrical source terms (e.g. for expanding/contracting box formalism) in the multi-fluid equation system for a
* single fluid, within a single cell, using an explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param U0 (Initial) comoving plasma velocity.
* @param R0 (Initial) radial distance from expansion/contraction center.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param em_old Array of old electromagnetic variables (before source update).
* @param em_new Array of new electromagnetic variables (after source update).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
*/
void
explicit_volume_source_maxwell_update_euler(const gkyl_moment_em_coupling* mom_em, const double U0, const double R0, double t_curr,
  const double dt, double* em_old, double* em_new, const double* ext_em);

/**
* Integrate the volume-based geometrical source terms (e.g. for expanding/contracting box formalism) in the multi-fluid equation system within
* a single cell, using an explicit forcing solver (specifically a strong stability-preserving third-order Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param em Array of electromagnetic variables.
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
*/
void
explicit_volume_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES],
  double* em, const double* ext_em);

/**
* Integrate the reactive source terms in the multi-fluid equation system within a single cell, using an explicit forcing solver (specifically a 
* simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param gas_gamma Adiabatic index.
* @param specific_heat_capacity Specific heat capacity.
* @param energy_of_formation Energy of formation.
* @param ignition_temperature Ignition temperature.
* @param reaction_rate Reaction rate.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_reactive_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, const double specific_heat_capacity,
  const double energy_of_formation, const double ignition_temperature, const double reaction_rate, double t_curr, const double dt,
  double* fluid_old, double* fluid_new);

/**
* Integrate the reactive source terms in the multi-fluid equation system within a single cell, using an explicit forcing solver (specifically a
* strong stability-preserving third-order Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
*/
void
explicit_reactive_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES]);

/**
* Integrate the coupled fluid-Einstein source terms in plane-symmetric spacetimes in the multi-fluid equation system within a single cell, using an
* explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param gas_gamma Adiabatic index.
* @param kappa Stress-energy prefactor in the Einstein field equations.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_medium_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, const double kappa, double t_curr,
  const double dt, double* fluid_old, double* fluid_new);

/**
* Integrate the coupled fluid-Einstein source terms in plane-symmetric spacetimes in the multi-fluid equation system within a single cell, using an
* explicit forcing solver (specifically a strong stability-preserving third-order Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
*/
void
explicit_medium_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES]);

/**
* Integrate the general relativistic source terms (Euler equations, ultra-relativistic equation of state) in the multi-fluid equation system within a
* single cell, using an explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param gas_gamma Adiabatic index.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_gr_ultra_rel_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, double t_curr, const double dt,
  double* fluid_old, double* fluid_new);

/**
* Integrate the general relativistic source terms (Euler equations, ultra-relativistic equation of state) in the multi-fluid equation system within a
* single cell, using an explicit forcing solver (specifically a strong stability-preserving third-order Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
*/
void
explicit_gr_ultra_rel_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES]);

/**
* Integrate the general relativistic source terms (Euler equations, general equation of state) in the multi-fluid equation system within a
* single cell, using an explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param gas_gamma Adiabatic index.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_gr_euler_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, double t_curr, const double dt,
  double* fluid_old, double* fluid_new);

/**
* Integrate the general relativistic source terms (Euler equations, general equation of state) in the multi-fluid equation system within a
* single cell, using an explicit forcing solver (specifically a strong stability-preserving third-order Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
*/
void
explicit_gr_euler_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES]);

/**
* Integrate the electron coupling source terms in the general relativistic two-fluid equation system within a single cell, using an explicit
* forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param gas_gamma_elc Adiabatic index (electrons).
* @param mass_elc Electron mass.
* @param charge_elc Electron charge.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_gr_twofluid_source_update_elc_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma_elc, const double mass_elc,
  const double charge_elc, double t_curr, const double dt, double* fluid_old, double* fluid_new);

/**
* Integrate the ion coupling source terms in the general relativistic two-fluid equation system within a single cell, using an explicit
* forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param gas_gamma_ion Adiabatic index (ions).
* @param mass_ion Ion mass.
* @param charge_ion Ion charge.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_gr_twofluid_source_update_ion_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma_ion, const double mass_ion,
  const double charge_ion, double t_curr, const double dt, double* fluid_old, double* fluid_new);

/**
* Integrate the electromagnetic coupling source terms in the general relativistic two-fluid equation system within a single cell, using an
* explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param gas_gamma_elc Adiabatic index (electrons).
* @param gas_gamma_ion Adiabatic index (ions).
* @param mass_elc Electron mass.
* @param charge_elc Electron charge.
* @param mass_ion Ion mass.
* @param charge_ion Ion charge.
* @param e_fact Electric field divergence error propagation speed.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_gr_twofluid_source_update_em_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma_elc, const double gas_gamma_ion,
  const double mass_elc, const double charge_elc, const double mass_ion, const double charge_ion, double e_fact, double t_curr, const double dt,
  double* fluid_old, double* fluid_new);

/**
* Integrate the curved spacetime coupling source terms for the electrons in the general relativistic two-fluid equation system within a single
* cell, using an explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param gas_gamma_elc Adiabatic index (electrons).
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_gr_twofluid_source_update_elc_spacetime_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma_elc, double t_curr,
  const double dt, double* fluid_old, double* fluid_new);

/**
* Integrate the curved spacetime coupling source terms for the ions in the general relativistic two-fluid equation system within a single
* cell, using an explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param gas_gamma_ion Adiabatic index (ions).
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_gr_twofluid_source_update_ion_spacetime_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma_ion, double t_curr,
  const double dt, double* fluid_old, double* fluid_new);

/**
* Integrate all electron, ion, and electromagnetic coupling source terms in the general relativistic two-fluid equation system within a
* single cell, using an explicit forcing solver (specifically a strong stability-preserving third-order Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
*/
void
explicit_gr_twofluid_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES]);

/**
* Integrate the vacuum Einstein source terms (Bona-Masso formalism) in the dynamic spacetime equation system within a single cell, using an
* explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
@ @param excision_threshold Excision threshold (lapse).
* @param spacetime_slicing Spacetime slicing condition.
* @param spacetime_evolution Spacetime evolution system.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_vacuum_einstein_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double excision_threshold,
  const enum gkyl_spacetime_slicing spacetime_slicing, const enum gkyl_spacetime_evolution spacetime_evolution, double t_curr, const double dt,
  double* fluid_old, double* fluid_new);

/**
* Integrate the vacuum Einstein source terms (Bona-Masso formalism) in the dynamic spacetime equation system within a single cell, using an
* explicit forcing solver (specifically a strong stability-preserving third-order Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
*/
void
explicit_vacuum_einstein_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES]);

/**
* Integrate the vacuum Einstein source terms (conformal Bona-Masso formalism) in the dynamic spacetime equation system within a single cell, using an
* explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
@ @param excision_threshold Excision threshold (lapse).
* @param spacetime_slicing Spacetime slicing condition.
* @param spacetime_evolution Spacetime evolution system.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_vacuum_einstein_conformal_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double excision_threshold,
  const enum gkyl_spacetime_slicing spacetime_slicing, const enum gkyl_spacetime_evolution spacetime_evolution, double t_curr, const double dt,
  double* fluid_old, double* fluid_new);

/**
* Integrate the vacuum Einstein source terms (conformal Bona-Masso formalism) in the dynamic spacetime equation system within a single cell, using an
* explicit forcing solver (specifically a strong stability-preserving third-order Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
*/
void
explicit_vacuum_einstein_conformal_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES]);

/**
* Integrate the general relativistic source terms (general relativistic magnetohydrodynamics equations) in the multi-fluid equation system within a
single cell, using an explicit forcing solver (specifically a simple first-order forward-Euler method).
*
* @param mom_em Moment-EM coupling object.
* @param gas_gamma Adiabatic index.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_old Array of old fluid variables (before source update).
* @param fluid_new Array of new fluid variables (after source update).
*/
void
explicit_gr_mhd_source_update_euler(const gkyl_moment_em_coupling* mom_em, const double gas_gamma, double t_curr, const double dt,
  double* fluid_old, double* fluid_new);

/**
* Integrate the general relativistic source terms (general relativistic magnetohydrodynamics) in the multi-fluid equation system within a
* single cell, using an explicit forcing solver (specifically a strong stability-preserving third-order Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
*/
void
explicit_gr_mhd_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, const double dt, double* fluid_s[GKYL_MAX_SPECIES]);

/**
* Integrate the electric field source terms in the multi-field equation system within a single cell, using an explicit forcing solver (specifically
* a simple first-order forward-Euler method), assuming a cold relativistic fluid.
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param e_field_old Array of old electric field variables (before source update).
* @param e_field_new Array of new electric field variables (after source update).
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
*/
void
explicit_e_field_source_update_euler(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double e_field_old[3], double* e_field_new,
  double* fluid_s[GKYL_MAX_SPECIES], const double* app_current);

/**
* Integrate the electric field source terms in the multi-fluid equation system within a single cell, using an explicit forcing solver (specifically
* a strong stability-preserving third-order Runge-Kutta method), assuming a cold relativistic fluid.
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param app_current1 Array of stage-1 current terms to be applied to the fluid equations (for stage-1 of external current driving).
* @param app_current2 Array of stage-2 current terms to be applied to the fluid equations (for stage-2 of external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
*/
void
explicit_e_field_source_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  double* em, const double* app_current, const double* app_current1, const double* app_current2, const double* ext_em);

/**
* Perform a Higuera-Cary particle push of the fluid momentum in the multi-fluid equation system within a single cell, assuming a cold relativistic
* fluid.
*
* @param vel Fluid species velocity vector.
* @param q Charge of fluid species.
* @param m Mass of fluid species.
* @param dt Current stable time-step.
* @param e_field Array of electric field variables.
* @param b_field Array of magnetic field variables.
*/
void
explicit_higuera_cary_push(double* vel, const double q, const double m, const double dt, const double c, const double e_field[3],
  const double b_field[3]);

/**
* Integrate the momentum source terms in the multi-fluid equation system within a single cell, using a Higuera-Cary particle push, assuming a cold
* relativistic fluid.
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param app_accel_s Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param em Array of electromagnetic variables.
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
*/
void
explicit_higuera_cary_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], double* em, const double* ext_em);

/**
* Integrate the electromagnetic source terms in the multi-fluid equation system within each cell, using an explicit forcing solver (specifically
* a strong stability-preserving third-order Runge-Kutta method combined with a Higuera-Cary update), assuming a cold relativistic fluid.
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation tine.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param app_accel_s Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param app_current1 Array of stage-1 current terms to be applied to the fluid equations (for stage-1 of external current driving).
* @param app_current2 Array of stage-2 current terms to be applied to the fluid equations (for stage-2 of external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param nstrang Indictator of which step in the Strang splitting we are currently considering.
*/
void
explicit_source_coupling_update(const gkyl_moment_em_coupling* mom_em, double t_curr, double dt, double* fluid_s[GKYL_MAX_SPECIES],
  const double* app_accel_s[GKYL_MAX_SPECIES], double* em, const double* app_current, const double* app_current1, const double* app_current2,
  const double* ext_em, int nstrang);