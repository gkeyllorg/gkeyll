#pragma once

#include <gkyl_wv_eqn.h>
#include <gkyl_gr_spacetime.h>

// Type of Riemann-solver to use:
enum gkyl_wv_gr_twofluid_rp {
  WV_GR_TWOFLUID_RP_HLL, // Default (HLL fluxes).
  WV_GR_TWOFLUID_RP_LAX
};

// Input context, packaged as a struct.
struct gkyl_wv_gr_twofluid_inp {
  double mass_elc; // Electron mass.
  double mass_ion; // Ion mass;
  double charge_elc; // Electron charge.
  double charge_ion; // Ion charge.
  double gas_gamma_elc; // Adiabatic index (electrons).
  double gas_gamma_ion; // Adiabatic index (ions).
  double light_speed; // Speed of light.
  double e_fact; // Factor of speed of light for electric field correction.
  double b_fact; // Factor of speed of light for magnetic field correction.

  enum gkyl_spacetime_gauge spacetime_gauge; // Spacetime gauge choice.
  int reinit_freq; // Spacetime reinitialization frequency.
  struct gkyl_gr_spacetime *spacetime; // Pointer to base spacetime object.

  enum gkyl_wv_gr_twofluid_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new general relativistic two-fluid equations object with ideal gas equation of state.
*
* @param mass_elc Electron mass.
* @param mass_ion Ion mass.
* @param charge_elc Electron charge.
* @param charge_ion Ion charge.
* @param gas_gamma_elc Adiabatic index (electrons).
* @param gas_gamma_ion Adiabatic index (ions).
* @param light_speed Speed of light.
* @param e_fact Factor of speed of light for electric field correction.
* @param b_fact Factor of speed of light for magnetic field correction.
* @param spacetime_gauge Spacetime gauge choice.
* @param reinit_freq Spacetime reinitialization frequency.
* @param spacetime Pointer to base spacetime object.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the general relativistic two-fluid equations object with ideal gas equation of state.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_twofluid_new(double mass_elc, double mass_ion, double charge_elc, double charge_ion, double gas_gamma_elc, double gas_gamma_ion,
  double light_speed, double e_fact, double b_fact, enum gkyl_spacetime_gauge spacetime_gauge, int reinit_freq, struct gkyl_gr_spacetime* spacetime, bool use_gpu);

/**
* Create a new general relativistic two-fluid equations object with ideal gas equation of state, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the general relativistic two-fluid equations object with ideal gas equation of state.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_twofluid_inew(const struct gkyl_wv_gr_twofluid_inp* inp);

/**
* Get electron mass.
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Electron mass.
*/
double
gkyl_wv_gr_twofluid_mass_elc(const struct gkyl_wv_eqn* eqn);

/**
* Get ion mass.
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Ion mass.
*/
double
gkyl_wv_gr_twofluid_mass_ion(const struct gkyl_wv_eqn* eqn);

/**
* Get electron charge.
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Electron charge.
*/
double
gkyl_wv_gr_twofluid_charge_elc(const struct gkyl_wv_eqn* eqn);

/**
* Get ion charge.
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Ion charge.
*/
double
gkyl_wv_gr_twofluid_charge_ion(const struct gkyl_wv_eqn* eqn);

/**
* Get adiabatic index (electrons).
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Adiabatic index (electrons).
*/
double
gkyl_wv_gr_twofluid_gas_gamma_elc(const struct gkyl_wv_eqn* eqn);

/**
* Get adiabatic index (ions).
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Adiabatic index (ions).
*/
double
gkyl_wv_gr_twofluid_gas_gamma_ion(const struct gkyl_wv_eqn* eqn);

/**
* Get speed of light.
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Speed of light.
*/
double
gkyl_wv_gr_twofluid_light_speed(const struct gkyl_wv_eqn* eqn);

/**
* Get factor of speed of light for electric field correction.
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Factor of speed of light for electric field correction.
*/
double
gkyl_wv_gr_twofluid_e_fact(const struct gkyl_wv_eqn* eqn);

/**
* Get factor of speed of light for magnetic field correction.
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Factor of speed of light for magnetic field correction.
*/
double
gkyl_wv_gr_twofluid_b_fact(const struct gkyl_wv_eqn* eqn);

/**
* Get spacetime gauge choice.
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Spacetime gauge choice.
*/
enum gkyl_spacetime_gauge
gkyl_wv_gr_twofluid_spacetime_gauge(const struct gkyl_wv_eqn* eqn);

/**
* Get spacetime reinitialization frequency.
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Spacetime reinitialization frequency.
*/
int
gkyl_wv_gr_twofluid_reinit_freq(const struct gkyl_wv_eqn* eqn);

/**
* Get base spacetime object.
*
* @param eqn General relativistic two-fluid equations object with ideal gas equation of state.
* @return Pointer to the base spacetime object.
*/
struct gkyl_gr_spacetime*
gkyl_wv_gr_twofluid_spacetime(const struct gkyl_wv_eqn* eqn);