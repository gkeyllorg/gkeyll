
#pragma once

// ** Fully formally-verified solvers for the inviscid Burgers' equation **
// ** Lax-Friedrichs Solver: **
// ** Proof of hyperbolicity preservation: ../proofs/finite_volume/proof_inviscid_burgers_lax_hyperbolicity.rkt **
// ** Proof of CFL stability: ../proofs/finite_volume/proof_inviscid_burgers_lax_cfl_stability.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function: ../proofs/finite_volume/proof_inviscid_burgers_lax_local_lipschitz.rkt **
// ** Roe Solver: **
// ** Proof of hyperbolicity preservation: ../proofs/finite_volume/proof_inviscid_burgers_roe_hyperbolicity.rkt **
// ** Proof of flux conservation (jump continuity): ../proofs/finite_volume/proof_inviscid_burgers_roe_flux_conservation.rkt **

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_burgers_rp {
  WV_BURGERS_RP_ROE = 0, // Default (Roe fluxes).
  WV_BURGERS_RP_LAX
};

// Input context, packaged as a struct.
struct gkyl_wv_burgers_inp {
  enum gkyl_wv_burgers_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new inviscid Burgers' equation object.
*
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the inviscid Burgers' equation object.
*/
struct gkyl_wv_eqn*
gkyl_wv_burgers_new(bool use_gpu);

/**
* Create a new inviscid Burgers' equation object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the inviscid Burgers' equation object.
*/
struct gkyl_wv_eqn*
gkyl_wv_burgers_inew(const struct gkyl_wv_burgers_inp* inp);
