#pragma once

// ** Partially formally-verified solvers for the isothermal Euler equations **
// ** Lax-Friedrichs Solver: **
// ** Proof of hyperbolicity preservation (rho and mom_x components): ../proofs/finite_volume/proof_isothermal_euler_mom_x_lax_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (rho and mom_x components): ../proofs/finite_volume/proof_isothermal_euler_mom_x_lax_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_strict_hyperbolicity.rkt **
// ** Proof of CFL stability (rho and mom_x components): ../proofs/finite_volume/proof_isothermal_euler_mom_x_lax_cfl_stability.rkt **
// ** Proof of CFL stability (mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_cfl_stability.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function (rho and mom_x components): NOT PROVEN **
// ** Proof of local Lipschitz continuity of discrete flux function (mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_lax_local_lipschitz.rkt **
// ** Roe Solver: **
// ** Proof of hyperbolicity preservation (rho and mom_x components): NOT PROVEN **
// ** Proof of hyperbolicity preservation (mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_roe_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (rho and mom_x components): NOT PROVEN **
// ** Proof of strict hyperbolicity preservation (mom_y and mom_z components): NOT PROVEN **
// ** Proof of flux conservation (jump continuity, rho and mom_x components): NOT PROVEN **
// ** Proof of flux conservation (jump continuity, mom_y and mom_z components): ../proofs/finite_volume/proof_isothermal_euler_mom_yz_roe_flux_conservation.rkt **

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_iso_euler_rp {
  WV_ISO_EULER_RP_LAX = 0, // Default (Lax fluxes).
  WV_ISO_EULER_RP_ROE,
};

// Input context, packaged as a struct.
struct gkyl_wv_iso_euler_inp {
  double vt; // Thermal velocity.

  enum gkyl_wv_iso_euler_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new isothermal Euler equations object.
*
* @param vt Thermal velocity.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the isothermal Euler equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_iso_euler_new(double vt, bool use_gpu);

/**
* Create a new isothermal Euler equations object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the isothermal Euler equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_iso_euler_inew(const struct gkyl_wv_iso_euler_inp* inp);

/**
 * Get thermal velocity.
 * 
 * @param wv Isothermal Euler equations object.
 * @return Thermal velocity.
 */
 double gkyl_wv_iso_euler_vt(const struct gkyl_wv_eqn* wv);