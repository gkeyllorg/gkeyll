#pragma once

// ** Fully formally-verified solvers for the perfectly hyperbolic Maxwell equations **
// ** Lax-Friedrichs Solver: **
// ** Proof of hyperbolicity preservation (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_strict_hyperbolicity.rkt **
// ** Proof of CFL stability (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_cfl_stability.rkt **
// ** Proof of CFL stability (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_cfl_stability.rkt **
// ** Proof of CFL stability (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_cfl_stability.rkt **
// ** Proof of CFL stability (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_cfl_stability.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_lax_local_lipschitz.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_lax_local_lipschitz.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_lax_local_lipschitz.rkt **
// ** Proof of local Lipschitz continuity of discrete flux function (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_lax_local_lipschitz.rkt **
// ** Roe Solver: **
// ** Proof of hyperbolicity preservation (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_hyperbolicity.rkt **
// ** Proof of hyperbolicity preservation (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_strict_hyperbolicity.rkt **
// ** Proof of strict hyperbolicity preservation (Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_strict_hyperbolicity.rkt **
// ** Proof of flux conservation (jump continuity, Ex and phi components): ../proofs/finite_volume/proof_maxwell_1d_Ex_phi_roe_flux_conservation.rkt **
// ** Proof of flux conservation (jump continuity, Ey and Bz components): ../proofs/finite_volume/proof_maxwell_1d_Ey_Bz_roe_flux_conservation.rkt **
// ** Proof of flux conservation (jump continuity, Ez and By components): ../proofs/finite_volume/proof_maxwell_1d_Ez_By_roe_flux_conservation.rkt **
// ** Proof of flux conservation (jump continuity, Bx and psi components): ../proofs/finite_volume/proof_maxwell_1d_Bx_psi_roe_flux_conservation.rkt **

#include <gkyl_wv_eqn.h>

// Type of Riemann-solver to use:
enum gkyl_wv_maxwell_rp {
  WV_MAXWELL_RP_ROE = 0, // Default (Roe fluxes).
  WV_MAXWELL_RP_LAX
};

// Input context, packaged as a struct.
struct gkyl_wv_maxwell_inp {
  double c; // Speed of light.
  double e_fact; // Factor of speed of light for electric field correction.
  double b_fact; // Factor of speed of light for magnetic field correction.

  struct gkyl_wv_embed_geo *embed_geo; // embedded geometry
  enum gkyl_wv_maxwell_rp rp_type; // Type of Riemann-solver to use.
  bool use_gpu; // Whether the wave equation object is on the host (false) or the device (true).
};

/**
* Create a new Maxwell equations object.
*
* @param c Speed of light.
* @param e_fact Factor of speed of light for electric field correction.
* @param b_fact Factor of speed of light for magnetic field correction.
* @param use_gpu Whether the wave equation object is on the host (false) or the device (true).
* @return Pointer to the Maxwell equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_maxwell_new(double c, double e_fact, double b_fact, bool use_gpu);

/**
 * Create a new Maxwell equation object that lives on NV-GPU.
 * see new() method above for documentation.
 */
 struct gkyl_wv_eqn* gkyl_wv_maxwell_cu_dev_new(double c, double e_fact, double b_fact);

/**
* Create a new Maxwell equations object, from an input context struct.
*
* @param inp Input context struct.
* @return Pointer to the Maxwell equations object.
*/
struct gkyl_wv_eqn*
gkyl_wv_maxwell_inew(const struct gkyl_wv_maxwell_inp* inp);
