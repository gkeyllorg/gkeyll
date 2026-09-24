#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_euler_rgfm_rp {
  WV_EULER_RGFM_RP_LAX = 0, // Default (Lax fluxes).
};

// Input context, packaged as a struct.
struct gkyl_wv_euler_rgfm_inp {
  int num_species; // Number of distinct species in the domain.
  double* gas_gamma_s; // Adiabatic indices for each species in the domain.
  int reinit_freq; // Reinitialization frequency for the level set.

  enum gkyl_wv_euler_rgfm_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new Euler Riemann ghost fluid equations object.
*
* @param num_species Number of distinct species in the domain.
* @param gas_gamma_s Adiabatic indices for each species in the domain.
* @param reinit_freq Reinitialization frequency for the level set.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the Euler Riemann ghost fluid equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_euler_rgfm_new(int num_species, double* gas_gamma_s, int reinit_freq, bool use_gpu);

/**
* Create a new Euler Riemann ghost fluid equations object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the Euler Riemann ghost fluid equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_euler_rgfm_inew(const struct gkyl_wv_euler_rgfm_inp* inp);

/**
* Get number of distinct species in the domain.
*
* @param wv Euler Riemann ghost fluid equations object.
* @return Number of distinct species in the domain.
*/
int
gkyl_wv_euler_rgfm_num_species(const struct gkyl_wv_eqn* wv);

/**
* Get adiabatic indices for each species in the domain.
*
* @param wv Euler Riemann ghost fluid equations object.
* @return Adiabatic indices for each species in the domain.
*/
double*
gkyl_wv_euler_rgfm_gas_gamma_s(const struct gkyl_wv_eqn* wv);

/**
* Get reinitialization frequency for the level set.
*
* @param wv Euler Riemann ghost fluid equations object.
* @return Reinitialization frequency for the level set.
*/
int
gkyl_wv_euler_rgfm_reinit_freq(const struct gkyl_wv_eqn* wv);