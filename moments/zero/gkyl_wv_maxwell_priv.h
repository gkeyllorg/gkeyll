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

// Private header, not for direct use in user-facing code.

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_embed_geo.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>

struct wv_maxwell {
  struct gkyl_wv_eqn eqn; // Base equation object.
  double c; // Speed of light.
  double e_fact; // Factor of speed of light for electric field correction.
  double b_fact; // Factor of speed of light for magnetic field correction.
};

/**
* Free Maxwell equations object.
*
* @param ref Reference counter for Maxwell equations.
*/
void
gkyl_wv_maxwell_free(const struct gkyl_ref_count* ref);

/**
* Compute maximum absolute wave speed.
*
* @param c Speed of light.
* @param e_fact Factor of speed of light for electric field correction.
* @param b_fact Factor of speed of light for magnetic field correction.
* @param q Conserved variable vector.
* @return Maximum absolute wave speed for a given q.
*/
GKYL_CU_DH
static inline double
gkyl_maxwell_max_abs_speed(double c, double e_fact, double b_fact, const double* q)
{
  return fmax(fmax(fabs((c * e_fact)), fabs(c)), fabs((b_fact * c)));
}

/**
* Compute flux vector. Assumes rotation to local coordinate system.
*
* @param c Speed of light.
* @param e_fact Factor of speed of light for electric field correction.
* @param b_fact Factor of speed of light for magnetic field correction.
* @param q Conserved variable vector.
* @param flux Flux vector in direction 'dir' (output).
*/
GKYL_CU_DH
static inline void
gkyl_maxwell_flux(double c, double e_fact, double b_fact, const double* q, double* flux)
{
  flux[0] = (e_fact * ((c * c) * q[6]));
  flux[1] = ((c * c) * q[5]);
  flux[2] = (-1.0 * ((c * c) * q[4]));
  flux[3] = (b_fact * q[7]);
  flux[4] = (-1.0 * q[2]);
  flux[5] = q[1];
  flux[6] = (e_fact * q[0]);
  flux[7] = (b_fact * ((c * c) * q[3]));
}

/**
* Compute eigenvalues of the flux Jacobian. Assumes rotation to local coordinate system.
*
* @param c Speed of light.
* @param e_fact Factor of speed of light for electric field correction.
* @param b_fact Factor of speed of light for magnetic field correction.
* @param q Conserved variable vector.
* @param flux_deriv Flux Jacobian eigenvalues in direction 'dir' (output).
*/
GKYL_CU_DH
static inline void
gkyl_maxwell_flux_deriv(double c, double e_fact, double b_fact, const double* q, double* flux_deriv)
{
  flux_deriv[0] = (-1.0 * (c * e_fact));
  flux_deriv[1] = (-1.0 * c);
  flux_deriv[2] = c;
  flux_deriv[3] = (-1.0 * (b_fact * c));
  flux_deriv[4] = (-1.0 * c);
  flux_deriv[5] = c;
  flux_deriv[6] = (c * e_fact);
  flux_deriv[7] = (b_fact * c);
}

/**
* Compute Riemann variables given the conserved variables.
*
* @param eqn Base equation object.
* @param qstate Current state vector.
* @param qin Conserved variable vector (input).
* @param wout Riemann variable vector (output).
*/
GKYL_CU_DH
static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 8; i++) {
    wout[i] = qin[i];
  }
}

/**
* Compute conserved variables given the Riemann variables.
*
* @param eqn Base equation object.
* @param qstate Current state vector.
* @param win Riemann variable vector (input).
* @param qout Conserved variable vector (output).
*/
GKYL_CU_DH
static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double *qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 8; i++) {
    qout[i] = win[i];
  }
}

/**
* Boundary condition function for applying wall boundary conditions for the Maxwell equations.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply wall boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
GKYL_CU_DH
static void
maxwell_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  // Zero tangent for electric field.
  ghost[0] = skin[0];
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];

  // Zero normal for magnetic field.
  ghost[3] = -skin[3];
  ghost[4] = skin[4];
  ghost[5] = skin[5];

  // Correction potentials.
  ghost[6] = -skin[6];
  ghost[7] = skin[7];
}

/**
* Boundary condition function for applying no-slip boundary conditions for the Maxwell equations.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply no-slip boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
GKYL_CU_DH
static void
maxwell_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  // Zero tangent for the electric field.
  ghost[0] = skin[0];
  ghost[1] = -skin[1];
  ghost[2] = -skin[2];

  // Zero normal for the magnetic field.
  ghost[3] = -skin[3];
  ghost[4] = skin[4];
  ghost[5] = skin[5];

  // Correction potentials.
  ghost[6] = -skin[6];
  ghost[7] = skin[7];
}

/**
* Rotate state vector from global to local coordinate frame.
*
* @param eqn Base equation object.
* @param tau1 First tangent vector of the coordinate frame.
* @param tau2 Second tangent vector of the coordinate frame.
* @param norm Normal vector of the coordinate frame.
* @param qglobal State vector in global coordinate frame (input).
* @param qlocal State vector in local coordinate frame (output).
*/
GKYL_CU_DH
static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal)
{
  // Rotate electric field vector to local coordinates.
  qlocal[0] = (qglobal[0] * norm[0]) + (qglobal[1] * norm[1]) + (qglobal[2] * norm[2]);
  qlocal[1] = (qglobal[0] * tau1[0]) + (qglobal[1] * tau1[1]) + (qglobal[2] * tau1[2]);
  qlocal[2] = (qglobal[0] * tau2[0]) + (qglobal[1] * tau2[1]) + (qglobal[2] * tau2[2]);
  
  // Rotate magnetic field vector to local coordinates.
  qlocal[3] = (qglobal[3] * norm[0]) + (qglobal[4] * norm[1]) + (qglobal[5] * norm[2]);
  qlocal[4] = (qglobal[3] * tau1[0]) + (qglobal[4] * tau1[1]) + (qglobal[5] * tau1[2]);
  qlocal[5] = (qglobal[3] * tau2[0]) + (qglobal[4] * tau2[1]) + (qglobal[5] * tau2[2]);
  
  // Correction potentials are scalars (so remain unchanged).
  qlocal[6] = qglobal[6];
  qlocal[7] = qglobal[7];
}

/**
* Rotate state vector from local to global coordinate frame.
*
* @param eqn Base equation object.
* @param tau1 First tangent vector of the coordinate frame.
* @param tau2 Second tangent vector of the coordinate frame.
* @param norm Normal vector of the coordinate frame.
* @param qlocal State vector in local coordinate frame (input).
* @param qglobal State vector in global coordinate frame (output).
*/
GKYL_CU_DH
static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal)
{
  // Rotate electric field vector to global coordinates.
  qglobal[0] = (qlocal[0] * norm[0]) + (qlocal[1] * tau1[0]) + (qlocal[2] * tau2[0]);
  qglobal[1] = (qlocal[0] * norm[1]) + (qlocal[1] * tau1[1]) + (qlocal[2] * tau2[1]);
  qglobal[2] = (qlocal[0] * norm[2]) + (qlocal[1] * tau1[2]) + (qlocal[2] * tau2[2]);

  // Rotate magnetic field vector to global coordinates.
  qglobal[3] = (qlocal[3] * norm[0]) + (qlocal[4] * tau1[0]) + (qlocal[5] * tau2[0]);
  qglobal[4] = (qlocal[3] * norm[1]) + (qlocal[4] * tau1[1]) + (qlocal[5] * tau2[1]);
  qglobal[5] = (qlocal[3] * norm[2]) + (qlocal[4] * tau1[2]) + (qlocal[5] * tau2[2]);

  // Correction potentials are scalars (so remain unchanged).
  qglobal[6] = qlocal[6];
  qglobal[7] = qlocal[7];
}

/**
* Compute waves and speeds using Lax fluxes.
*
* @param eqn Base equation object.
* @param delta Jump across interface to split.
* @param ql Conserved variables on the left of the interface.
* @param qr Conserved variables on the right of the interface.
* @param waves Waves (output).
* @param s Wave speeds (output).
* @return Maximum wave speed.
*/
GKYL_CU_DH
static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_maxwell *maxwell = container_of(eqn, struct wv_maxwell, eqn);
  double c = maxwell->c; // Speed of light.
  double e_fact = maxwell->e_fact; // Factor of speed of light for electric field correction.
  double b_fact = maxwell->b_fact; // Factor of speed of light for magnetic field correction.

  double sl = gkyl_maxwell_max_abs_speed(c, e_fact, b_fact, ql);
  double sr = gkyl_maxwell_max_abs_speed(c, e_fact, b_fact, qr);
  double amax = fmax(sl, sr);

  double fl[8], fr[8];

  fl[0] = (e_fact * ((c * c) * ql[6]));
  fl[1] = ((c * c) * ql[5]);
  fl[2] = (-1.0 * ((c * c) * ql[4]));
  fl[3] = (b_fact * ql[7]);
  fl[4] = (-1.0 * ql[2]);
  fl[5] = ql[1];
  fl[6] = (e_fact * ql[0]);
  fl[7] = (b_fact * ((c * c) * ql[3]));

  fr[0] = (e_fact * ((c * c) * qr[6]));
  fr[1] = ((c * c) * qr[5]);
  fr[2] = (-1.0 * ((c * c) * qr[4]));
  fr[3] = (b_fact * qr[7]);
  fr[4] = (-1.0 * qr[2]);
  fr[5] = qr[1];
  fr[6] = (e_fact * qr[0]);
  fr[7] = (b_fact * ((c * c) * qr[3]));

  double *w0 = &waves[0], *w1 = &waves[8];
  for (int i = 0; i < 8; i++) {
    w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
    w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
  }

  s[0] = -amax;
  s[1] = amax;


  return s[1];
}

/**
* Compute fluctuations using Lax fluxes.
*
* @param eqn Base equation object.
* @param ql Conserved variable vector on the left of the interface.
* @param qr Conserved variable vector on the right of the interface.
* @param waves Waves (input).
* @param s Wave speeds (input).
* @param amdq Left-moving fluctuations (output).
* @param apdq Right-moving fluctuations (output).
*/
GKYL_CU_DH
static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0], *w1 = &waves[8];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 8; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

GKYL_CU_DH
static void
wave_embed_copy_B(const double *q, double *qphi, double *delta, void *ctx)
{
  qphi[0] = 0.0;
  qphi[1] = 0.0;
  qphi[2] = 0.0;
  qphi[3] = q[3];
  qphi[4] = q[4];
  qphi[5] = q[5];
  qphi[6] = q[6];
  qphi[7] = q[7];
}

GKYL_CU_DH
static void
wave_embed_pec(const double *q, double *qphi, double *delta, void *ctx)
{
  qphi[0] = q[0];
  qphi[1] = -q[1];
  qphi[2] = -q[2];
  qphi[3] = -q[3];
  qphi[4] = q[4];
  qphi[5] = q[5];
  qphi[6] = -q[6];
  qphi[7] = q[7];
}

GKYL_CU_DH
static double
wave_embedded(const struct gkyl_wv_eqn *eqn, const double *delta,
  const double *ql, const double *qr, const double phil, const double phir,
  double *waves, double *s)
{
  double qphi[8] = {0.0};
  double deltaphi[8] = {0.0};

  if ((phil < 0.0) && (phir > 0.0)) {
    eqn->embed_geo->embed_func(qr, qphi, deltaphi, eqn->embed_geo->ctx);

    for (int d=0; d<8; ++d) deltaphi[d] = qr[d] - qphi[d];
  }
  if ((phil > 0.0) && (phir < 0.0)) {
    eqn->embed_geo->embed_func(ql, qphi, deltaphi, eqn->embed_geo->ctx);

    for (int d=0; d<8; ++d) deltaphi[d] = qphi[d] - ql[d];
  }

  return wave_lax(eqn, deltaphi, ql, qr, waves, s);
}

/**
* Compute waves and speeds using Lax fluxes (with potential fallback).
*
* @param eqn Base equation object.
* @param type Type of Riemann-solver flux to use.
* @param delta Jump across interface to split.
* @param ql Conserved variables on the left of the interface.
* @param qr Conserved variables on the right of the interface.
* @param waves Waves (output).
* @param s Wave speeds (output).
* @return Maximum wave speed.
*/
GKYL_CU_DH
static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, double phil, double phir, double* waves, double* s)
{
  if ((phil < 0.0) || (phir < 0.0))
    return wave_embedded(eqn, delta, ql, qr, phil, phir, waves, s);
  else
    return wave_lax(eqn, delta, ql, qr, waves, s);
}

/**
* Compute fluctuations using Lax fluxes (with potential fallback),
*
* @param eqn Base equation object.
* @param type Type of Riemann-solver flux to use.
* @param ql Conserved variable vector on the left of the interface.
* @param qr Conserved variable vector on the right of the interface.
* @param waves Waves (input).
* @param s Wave speeds (input).
* @param amdq Left-moving fluctuations (output).
* @param apdq Right-moving fluctuations (output).
*/
GKYL_CU_DH
static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, double phil, double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

/**
* Compute waves and speeds using Roe fluxes.
*
* @param eqn Base equation object.
* @param delta Jump across interface to split.
* @param ql Conserved variables on the left of the interface.
* @param qr Conserved variables on the right of the interface.
* @param waves Waves (output).
* @param s Wave speeds (output).
* @return Maximum wave speed.
*/
GKYL_CU_DH
static double
wave_roe(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_maxwell *maxwell = container_of(eqn, struct wv_maxwell, eqn);
  double c = maxwell->c; // Speed of light.
  double e_fact = maxwell->e_fact; // Factor of speed of light for electric field correction.
  double b_fact = maxwell->b_fact; // Factor of speed of light for magnetic field correction.

  double a1 = 0.5 * (delta[3] - ((1.0 / c) * delta[7]));
  double a2 = 0.5 * (delta[3] + ((1.0 / c) * delta[7]));
  double a3 = 0.5 * (delta[0] - (c * delta[6]));
  double a4 = 0.5 * (delta[0] + (c * delta[6]));
  double a5 = 0.5 * (delta[1] - (c * delta[5]));
  double a6 = 0.5 * (delta[1] + (c * delta[5]));
  double a7 = 0.5 * (delta[2] - (c * delta[4]));
  double a8 = 0.5 * (delta[2] + (c * delta[4]));

  double *w0 = &waves[0 * 8], *w1 = &waves[1 * 8], *w2 = &waves[2 * 8], *w3 = &waves[3 * 8], *w4 = &waves[4 * 8], *w5 = &waves[5 * 8];
  for (int i = 0; i < 8; i++) {
    w0[i] = 0.0; w1[i] = 0.0; w2[i] = 0.0; w3[i] = 0.0; w4[i] = 0.0; w5[i] = 0.0;
  }

  w0[3] = a1; w0[7] = -c * a1;
  s[0] = -c * b_fact;

  w1[3] = a2; w1[7] = c * a2;
  s[1] = c * b_fact;

  w2[0] = a3; w2[6] = -(1.0 / c) * a3;
  s[2] = -c * e_fact;

  w3[0] = a4; w3[6] = (1.0 / c) * a4;
  s[3] = c * e_fact;

  w4[1] = a5; w4[2] = a8;
  w4[4] = (1.0 / c) * a8; w4[5] = -(1.0 / c) * a5;
  s[4] = -c;

  w5[1] = a6; w5[2] = a7;
  w5[4] = -(1.0 / c) * a7; w5[5] = (1.0 / c) * a6;
  s[5] = c;

  return c;
}

/**
* Compute fluctuations using Roe fluxes.
*
* @param eqn Base equation object.
* @param ql Conserved variable vector on the left of the interface.
* @param qr Conserved variable vector on the right of the interface.
* @param waves Waves (input).
* @param s Wave speeds (input).
* @param amdq Left-moving fluctuations (output).
* @param apdq Right-moving fluctuations (output).
*/
GKYL_CU_DH
static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const double *w0 = &waves[0 * 8], *w1 = &waves[1 * 8], *w2 = &waves[2 * 8], *w3 = &waves[3 * 8], *w4 = &waves[4 * 8], *w5 = &waves[5 * 8];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]), s2m = fmin(0.0, s[2]), s3m = fmin(0.0, s[3]), s4m = fmin(0.0, s[4]), s5m = fmin(0.0, s[5]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]), s2p = fmax(0.0, s[2]), s3p = fmax(0.0, s[3]), s4p = fmax(0.0, s[4]), s5p = fmax(0.0, s[5]);

  for (int i = 0; i < 8; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]) + (s2m * w2[i]) + (s3m * w3[i]) + (s4m * w4[i]) + (s5m * w5[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]) + (s2p * w2[i]) + (s3p * w3[i]) + (s4p * w4[i]) + (s5p * w5[i]);
  }
}

/**
* Compute waves and speeds using Roe fluxes (with potential fallback).
*
* @param eqn Base equation object.
* @param type Type of Riemann-solver flux to use.
* @param delta Jump across interface to split.
* @param ql Conserved variables on the left of the interface.
* @param qr Conserved variables on the right of the interface.
* @param waves Waves (output).
* @param s Wave speeds (output).
* @return Maximum wave speed.
*/
GKYL_CU_DH
static double
wave(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, double phil, double phir, double* waves, double* s)
{
  // clear waves and wave speeds
  int mwaves = (type == GKYL_WV_HIGH_ORDER_FLUX) ? eqn->num_waves : 2;
  int meqn = eqn->num_equations;
  for (int i=0; i<mwaves; ++i) {
    double *w = &waves[i*meqn];
    for (int j=0; j<meqn; ++j)
      w[j] = 0.0;
    s[i] = 0.0;
  }
  // Check if inside an embedded surface
  if ((phil < 0.0) || (phir < 0.0))
    return wave_embedded(eqn, delta, ql, qr, phil, phir, waves, s);
  else {
    if (type == GKYL_WV_HIGH_ORDER_FLUX)
      return wave_roe(eqn, delta, ql, qr, waves, s);
    else
      return wave_lax(eqn, delta, ql, qr, waves, s);
  }

  return 0.0; // Unreachable code.
}

/**
* Compute fluctuations using Roe fluxes (with potential fallback),
*
* @param eqn Base equation object.
* @param type Type of Riemann-solver flux to use.
* @param ql Conserved variable vector on the left of the interface.
* @param qr Conserved variable vector on the right of the interface.
* @param waves Waves (input).
* @param s Wave speeds (input).
* @param amdq Left-moving fluctuations (output).
* @param apdq Right-moving fluctuations (output).
*/
GKYL_CU_DH
static void
qfluct(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, double phil, double phir, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  if (type == GKYL_WV_HIGH_ORDER_FLUX) {
    return qfluct_roe(eqn, ql, qr, waves, s, amdq, apdq);
  }
  else {
    return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
  }
}

/**
* Compute jump in flux given two conserved variable states.
*
* @param eqn Base equation object.
* @param ql Conserved variable vector on the left of the interface (input).
* @param qr Conserved variable vector on the right of the interface (input).
* @param flux_jump Jump in flux vector (output).
* @return Maximum wave speeds for states ql and qr.
*/
GKYL_CU_DH
static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_maxwell *maxwell = container_of(eqn, struct wv_maxwell, eqn);
  double c = maxwell->c; // Speed of light.
  double e_fact = maxwell->e_fact; // Factor of speed of light for electric field correction.
  double b_fact = maxwell->b_fact; // Factor of speed of light for magnetic field correction.

  double fl[8], fr[8];

  fl[0] = (e_fact * ((c * c) * ql[6]));
  fl[1] = ((c * c) * ql[5]);
  fl[2] = (-1.0 * ((c * c) * ql[4]));
  fl[3] = (b_fact * ql[7]);
  fl[4] = (-1.0 * ql[2]);
  fl[5] = ql[1];
  fl[6] = (e_fact * ql[0]);
  fl[7] = (b_fact * ((c * c) * ql[3]));

  fr[0] = (e_fact * ((c * c) * qr[6]));
  fr[1] = ((c * c) * qr[5]);
  fr[2] = (-1.0 * ((c * c) * qr[4]));
  fr[3] = (b_fact * qr[7]);
  fr[4] = (-1.0 * qr[2]);
  fr[5] = qr[1];
  fr[6] = (e_fact * qr[0]);
  fr[7] = (b_fact * ((c * c) * qr[3]));

  for (int i = 0; i < 8; i++) {
    flux_jump[i] = fr[i] - fl[i];
  }

  double amaxl = gkyl_maxwell_max_abs_speed(c, e_fact, b_fact, ql);
  double amaxr = gkyl_maxwell_max_abs_speed(c, e_fact, b_fact, qr);

  return fmax(amaxl, amaxr);
}

/**
* Determine whether invariant domain of the Maxwell equations is satisfied.
*
* @param eqn Base equation object.
* @param q Conserved variable vector.
* @return Whether the invariant domain is satisfied.
*/
GKYL_CU_DH
static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  return true; // All states are assumed to be valid.
}

/**
* Compute maximum wave speed from a conserved variable vector.
*
* @param eqn Base equation object.
* @param q Conserved variable vector.
* @return Maximum absolute wave speed.
*/
GKYL_CU_DH
static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_maxwell *maxwell = container_of(eqn, struct wv_maxwell, eqn);
  double c = maxwell->c; // Speed of light.
  double e_fact = maxwell->e_fact; // Factor of speed of light for electric field correction.
  double b_fact = maxwell->b_fact; // Factor of speed of light for magnetic field correction.

  return gkyl_maxwell_max_abs_speed(c, e_fact, b_fact, q);
}

/**
* Convert conserved variables to diagnostic variables.
*
* @param eqn Base equation object.
* @param qin Conserved variable vector (input).
* @param diag Diagnostic variable vector (output).
*/
GKYL_CU_DH
static inline void
maxwell_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 6; i++) {
    diag[i] = qin[i] * qin[i];
  }
}

/**
* Compute forcing/source term vector from conserved variable vector.
*
* @param eqn Base equation object.
* @param qin Conserved variable vector (input).
* @param sout Forcing/source term vector (output).
*/
GKYL_CU_DH
static inline void
maxwell_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout);
