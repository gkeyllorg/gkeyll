#pragma once

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_euler_mixture_rp {
  WV_EULER_MIXTURE_RP_LAX = 0, // Default (Lax fluxes).
  WV_EULER_MIXTURE_RP_ROE,
};

// Input context, packaged as a struct.
struct gkyl_wv_euler_mixture_inp {
  int num_species; // Number of distinct species in mixture.
  double* gas_gamma_s; // Adiabatic indices for each species in mixture.

  enum gkyl_wv_euler_mixture_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new Euler mixture equations object.
*
* @param num_species Number of distinct species in mixture.
* @param gas_gamma_s Adiabatic indices for each species in mixture.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the Euler mixture equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_euler_mixture_new(int num_species, double* gas_gamma_s, bool use_gpu);

/**
* Create a new Euler mixture equations object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the Euler mixture equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_euler_mixture_inew(const struct gkyl_wv_euler_mixture_inp* inp);

/**
* Get number of distinct species in mixture.
*
* @param wv Euler mixture equations object.
* @return Number of distinct species in mixture.
*/
int
gkyl_wv_euler_mixture_num_species(const struct gkyl_wv_eqn* wv);

/**
* Get adiabatic indices for each species in mixture.
*
* @param wv Euler mixture equations object.
* @return Adiabatic indices for each species in mixture.
*/
double*
gkyl_wv_euler_mixture_gas_gamma_s(const struct gkyl_wv_eqn* wv);