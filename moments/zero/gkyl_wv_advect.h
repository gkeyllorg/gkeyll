#pragma once

// ** Fully formally-verified solvers for the linear advection equation **
// ** Lax-Friedrichs Solver: **
// ** Proof of hyperbolicity preservation: ../proofs/finite_volume/proof_linear_advection_lax_hyperbolicity.rkt **
// ** Proof of CFL stability: ../proofs/finite_volume/proof_linear_advection_lax_cfl_stability.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function: ../proofs/finite_volume/proof_linear_advection_lax_local_lipschitz.rkt **
// ** Roe Solver: **
// ** Proof of hyperbolicity preservation: ../proofs/finite_volume/proof_linear_advection_roe_hyperbolicity.rkt **
// ** Proof of flux conservation (jump continuity): ../proofs/finite_volume/proof_linear_advection_roe_flux_conservation.rkt **

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_advect_rp {
  WV_ADVECT_RP_ROE = 0, // Default (Roe fluxes).
  WV_ADVECT_RP_LAX
};

// Input context, packaged as a struct.
struct gkyl_wv_advect_inp {
  double a; // Advection speed.

  enum gkyl_wv_advect_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new linear advection equation object.
*
* @param a Advection speed.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the linear advection equation object.
*/
struct gkyl_wv_eqn*
gkyl_wv_advect_new(double a, bool use_gpu);

/**
* Create a new linear advection equation object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the linear advection equation object.
*/
struct gkyl_wv_eqn*
gkyl_wv_advect_inew(const struct gkyl_wv_advect_inp* inp);
