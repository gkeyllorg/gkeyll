#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_reactive_euler_rp {
  WV_REACTIVE_EULER_RP_LAX = 0, // Default (Lax fluxes).
  WV_REACTIVE_EULER_RP_ROE,
};

// Input context, packaged as a struct.
struct gkyl_wv_reactive_euler_inp {
  double gas_gamma; // Adiabatic index.
  double specific_heat_capacity; // Specific heat capacity.
  double energy_of_formation; // Energy of formation.
  double ignition_temperature; // Ignition temperature.
  double reaction_rate; // Reaction rate.

  enum gkyl_wv_reactive_euler_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new reactive Euler equations object.
*
* @param gas_gamma Adiabatic index.
* @param specific_heat_capacity Specific heat capacity.
* @param energy_of_formation Energy of formation.
* @param ignition_temperature Ignition temperature.
* @param reaction_rate Reaction rate.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the reactive Euler equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_reactive_euler_new(double gas_gamma, double specific_heat_capacity, double energy_of_formation, double ignition_temperature,
  double reaction_rate, bool use_gpu);

/**
* Create a new reactive Euler equations object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the reactive Euler equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_reactive_euler_inew(const struct gkyl_wv_reactive_euler_inp* inp);

/**
* Get adiabatic index.
*
* @param wv Reactive Euler equations object.
* @return Adiabatic index.
*/
double
gkyl_wv_reactive_euler_gas_gamma(const struct gkyl_wv_eqn* wv);

/**
* Get specific heat capacity.
*
* @param wv Reactive Euler equations object.
* @return Specific heat capacity.
*/
double
gkyl_wv_reactive_euler_specific_heat_capacity(const struct gkyl_wv_eqn* wv);

/**
* Get energy of formation.
*
* @param wv Reactive Euler equations object.
* @return Energy of formation.
*/
double
gkyl_wv_reactive_euler_energy_of_formation(const struct gkyl_wv_eqn* wv);

/**
* Get ignition temperature.
*
* @param wv Reactive Euler equations object.
* @return Ignition temperature.
*/
double
gkyl_wv_reactive_euler_ignition_temperature(const struct gkyl_wv_eqn* wv);

/**
* Get reaction rate.
*
* @param wv Reactive Euler equations object.
* @return Reaction rate.
*/
double
gkyl_wv_reactive_euler_reaction_rate(const struct gkyl_wv_eqn* wv);