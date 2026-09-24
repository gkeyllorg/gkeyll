#pragma once

#include <gkyl_wv_eqn.h>
#include <gkyl_gr_spacetime.h>

// Type of Riemann-solver to use:
enum gkyl_wv_gr_maxwell_rp {
  WV_GR_MAXWELL_RP_ROE = 0, // Default (Roe fluxes).
  WV_GR_MAXWELL_RP_LAX,
};

// Input context, packaged as a struct.
struct gkyl_wv_gr_maxwell_inp {
  double light_speed; // Speed of light.
  double e_fact; // Factor of speed of light for electric field correction.
  double b_fact; // Factor of speed of light for magnetic field correction.

  enum gkyl_spacetime_gauge spacetime_gauge; // Spacetime gauge choice.
  int reinit_freq; // Spacetime reinitialization frequency.
  struct gkyl_gr_spacetime *spacetime; // Pointer to base spacetime object.

  enum gkyl_wv_gr_maxwell_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new general relativistic Maxwell equations object.
*
* @param light_speed Speed of light.
* @param e_fact Factor of speed of light for electric field correction.
* @param b_fact Factor of speed of light for magnetic field correction.
* @param spacetime_gauge Spacetime gauge choice.
* @param reinit_freq Spacetime reinitialization frequency.
* @param spacetime Pointer to base spacetime object.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the general relativistic Maxwell equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_maxwell_new(double light_speed, double e_fact, double b_fact, enum gkyl_spacetime_gauge spacetime_gauge, int reinit_freq, struct gkyl_gr_spacetime* spacetime, bool use_gpu);

/**
* Create a new general relativistic Maxwell equations object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the general relativistic Maxwell equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_gr_maxwell_inew(const struct gkyl_wv_gr_maxwell_inp* inp);

/**
* Get speed of light.
*
* @param eqn General relativistic Maxwell equations object.
* @return Speed of light.
*/
double
gkyl_wv_gr_maxwell_light_speed(const struct gkyl_wv_eqn* eqn);

/**
* Get factor of speed of light for electric field correction.
*
* @param eqn General relativistic Maxwell equations object.
* @return Factor of speed of light for electric field correction.
*/
double
gkyl_wv_gr_maxwell_e_fact(const struct gkyl_wv_eqn* eqn);

/**
* Get factor of speed of light for magnetic field correction.
*
* @param eqn General relativistic Maxwell equations object.
* @return Factor of speed of light for magnetic field correction.
*/
double
gkyl_wv_gr_maxwell_b_fact(const struct gkyl_wv_eqn* eqn);

/**
* Get spacetime gauge choice.
*
* @param eqn General relativistic Maxwell equations object.
* @return Spacetime gauge choice.
*/
enum gkyl_spacetime_gauge
gkyl_wv_gr_maxwell_spacetime_gauge(const struct gkyl_wv_eqn* eqn);

/**
* Get spacetime reinitialization frequency.
*
* @param eqn General relativistic Maxwell equations object.
* @return Spacetime reinitialization frequency.
*/
int
gkyl_wv_gr_maxwell_reinit_freq(const struct gkyl_wv_eqn* eqn);

/**
* Get base spacetime object.
*
* @param eqn General relativistic Maxwell equations object.
* @return Pointer to the base spacetime object.
*/
struct gkyl_gr_spacetime*
gkyl_wv_gr_maxwell_spacetime(const struct gkyl_wv_eqn* eqn);