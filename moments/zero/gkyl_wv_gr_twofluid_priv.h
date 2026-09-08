#pragma once

// Private header, not for direct use in user-facing code.

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct wv_gr_twofluid {
  struct gkyl_wv_eqn eqn; // Base equation object.
  struct gkyl_gr_spacetime *spacetime; // Pointer to base spacetime object.

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
};

/**
* Compute flux vector. Assumes rotation to local coordinate system.
*
* @param gas_gamma_elc Adiabatic index (electrons).
* @param gas_gamma_ion Adiabatic index (ions).
* @param light_speed Speed of light.
* @param e_fact Factor of speed of light for electric field correction.
* @param b_fact Factor of speed of light for magnetic field correction.
* @param q Conserved variable vector.
* @param flux Flux vector in direction 'dir' (output).
*/
GKYL_CU_D void gkyl_gr_twofluid_flux(double gas_gamma_elc, double gas_gamma_ion, double light_speed,
  double e_fact, double b_fact, const double q[84], double flux[84]);

/**
* Compute primitive variables given the conserved variables.
*
* @param gas_gamma_elc Adiabatic index (electrons).
* @param gas_gamma_ion Adiabatic index (ions).
* @param q Conserved variable vector.
* @param v Primitive variable vector (output).
*/
GKYL_CU_D void gkyl_gr_twofluid_prim_vars(
  double gas_gamma_elc, double gas_gama_ion, const double q[84], double v[84]);

/**
* Compute inverse spatial metric tensor (in covariant component form) given the conserved variables.
*
* @param q Conserved variable vector.
* @param inv_spatial_metric Inverse spatial metric tensor (output).
*/
GKYL_CU_D void gkyl_gr_twofluid_inv_spatial_metric(
  const double q[84], double ***inv_spatial_metric);

/**
* Compute perfect fluid stress-energy tensor (in contravariant component form) for the electrons, given the conserved variables.
*
* @param gas_gamma_elc Adiabatic index (electrons).
* @param gas_gamma_ion Adiabatic index (ions).
* @param q Conserved variable vector.
* @param stress_energy Stress-energy tensor (output).
*/
GKYL_CU_D void gkyl_gr_twofluid_stress_energy_tensor_elc(
  double gas_gamma_elc, double gas_gamma_ion, const double q[84], double ***stress_energy_elc);

/**
* Compute perfect fluid stress-energy tensor (in contravariant component form) for the ions, given the conserved variables.
*
* @param gas_gamma_elc Adiabatic index (electrons).
* @param gas_gamma_ion Adiabatic index (ions).
* @param q Conserved variable vector.
* @param stress_energy Stress-energy tensor (output).
*/
GKYL_CU_D void gkyl_gr_twofluid_stress_energy_tensor_ion(
  double gas_gamma_elc, double gas_gamma_ion, const double q[84], double ***stress_energy_ion);

/**
* Compute maximum absolute wave speed.
*
* @param gas_gamma_elc Adiabatic index (electrons).
* @param gas_gamma_ion Adiabatic index (ions).
* @param light_speed Speed of light.
* @param q Conserved variable vector.
* @return Maximum absolute wave speed for a given q.
*/
GKYL_CU_D static inline double gkyl_gr_twofluid_max_abs_speed(
  double gas_gamma_elc, double gas_gamma_ion, double light_speed, const double q[84]);

/**
* Compute Riemann variables given the conserved variables.
*
* @param eqn Base equation object.
* @param qstate Current state vector.
* @param qin Conserved variable vector (input).
* @param wout Riemann variable vector (output).
*/
GKYL_CU_D static inline void cons_to_riem(
  const struct gkyl_wv_eqn *eqn, const double *qstate, const double *qin, double *wout);

/**
* Compute conserved variables given the Riemann variables.
*
* @param eqn Base equation object.
* @param qstate Current state vector.
* @param win Riemann variable vector (input).
* @param qout Conserved variable vector (output).
*/
GKYL_CU_D static inline void riem_to_cons(
  const struct gkyl_wv_eqn *eqn, const double *qstate, const double *win, double *qout);

/**
* Boundary condition function for applying wall boundary conditions for the general relativistic two-fluid equations with ideal gas equation of state.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply wall boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
GKYL_CU_D static void gr_twofluid_wall(const struct gkyl_wv_eqn *eqn, double t, int nc,
  const double *skin, double *GKYL_RESTRICT ghost, void *ctx);

/**
* Boundary condition function for applying no-slip boundary conditions for the general relativistic two-fluid equations with ideal gas equation of state.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply no-slip boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
GKYL_CU_D static void gr_twofluid_no_slip(const struct gkyl_wv_eqn *eqn, double t, int nc,
  const double *skin, double *GKYL_RESTRICT ghost, void *ctx);

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
GKYL_CU_D static inline void rot_to_local(const struct gkyl_wv_eqn *eqn, const double *tau1,
  const double *tau2, const double *norm, const double *GKYL_RESTRICT qglobal,
  double *GKYL_RESTRICT qlocal);

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
GKYL_CU_D static inline void rot_to_global(const struct gkyl_wv_eqn *eqn, const double *tau1,
  const double *tau2, const double *norm, const double *GKYL_RESTRICT qlocal,
  double *GKYL_RESTRICT qglobal);

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
GKYL_CU_D static double wave_lax(const struct gkyl_wv_eqn *eqn, const double *delta,
  const double *ql, const double *qr, double *waves, double *s);

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
GKYL_CU_D static void qfluct_lax(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq);

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
GKYL_CU_D static double wave_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, const double phil, const double phir,
  double *waves, double *s);

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
GKYL_CU_D static void qfluct_lax_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves,
  const double *s, double *amdq, double *apdq);

/**
* Compute waves and speeds using HLL fluxes.
*
* @param eqn Base equation object.
* @param delta Jump across interface to split.
* @param ql Conserved variables on the left of the interface.
* @param qr Conserved variables on the right of the interface.
* @param waves Waves (output).
* @param s Wave speeds (output).
* @return Maximum wave speed.
*/
GKYL_CU_D static double wave_hll(const struct gkyl_wv_eqn *eqn, const double *delta,
  const double *ql, const double *qr, double *waves, double *s);

/**
* Compute fluctuations using HLL fluxes.
*
* @param eqn Base equation object.
* @param ql Conserved variable vector on the left of the interface.
* @param qr Conserved variable vector on the right of the interface.
* @param waves Waves (input).
* @param s Wave speeds (input).
* @param amdq Left-moving fluctuations (output).
* @param apdq Right-moving fluctuations (output).
*/
GKYL_CU_D static void qfluct_hll(const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr,
  const double *waves, const double *s, double *amdq, double *apdq);

/**
* Compute waves and speeds using HLL fluxes (with potential fallback).
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
GKYL_CU_D static double wave_hll_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *delta, const double *ql, const double *qr, const double phil, const double phir,
  double *waves, double *s);

/**
* Compute fluctuations using HLL fluxes (with potential fallback),
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
GKYL_CU_D static void qfluct_hll_l(const struct gkyl_wv_eqn *eqn, enum gkyl_wv_flux_type type,
  const double *ql, const double *qr, const double phil, const double phir, const double *waves,
  const double *s, double *amdq, double *apdq);

/**
* Compute jump in flux given two conserved variable states.
*
* @param eqn Base equation object.
* @param ql Conserved variable vector on the left of the interface (input).
* @param qr Conserved variable vector on the right of the interface (input).
* @param flux_jump Jump in flux vector (output).
* @return Maximum wave speeds for states ql and qr.
*/
GKYL_CU_D static double flux_jump(
  const struct gkyl_wv_eqn *eqn, const double *ql, const double *qr, double *flux_jump);

/**
* Determine whether invariant domain of the general relativistic two-fluid equations with ideal gas equation of state is satisfied.
*
* @param eqn Base equation object.
* @param q Conserved variable vector.
* @return Whether the invariant domain is satisfied.
*/
GKYL_CU_D static bool check_inv(const struct gkyl_wv_eqn *eqn, const double *q);

/**
* Compute maximum wave speed from a conserved variable vector.
*
* @param eqn Base equation object.
* @param q Conserved variable vector.
* @return Maximum absolute wave speed.
*/
GKYL_CU_D static double max_speed(const struct gkyl_wv_eqn *eqn, const double *q);

/**
* Convert conserved variables to diagnostic variables.
*
* @param eqn Base equation object.
* @param qin Conserved variable vector (input).
* @param diag Diagnostic variable vector (output).
*/
GKYL_CU_D static inline void gr_twofluid_cons_to_diag(
  const struct gkyl_wv_eqn *eqn, const double *qin, double *diag);

/**
* Compute forcing/source term vector from conserved variables.
*
* @param eqn Base equation object.
* @param qin Conserved variable vector (input).
* @param sout Forcing/source term vector (output).
*/
GKYL_CU_DH static inline void gr_twofluid_source(
  const struct gkyl_wv_eqn *eqn, const double *qin, double *sout);

/**
* Free general relativistic two-fluid equations object with ideal gas equation of state.
*
* @param ref Reference counter for general relativistic two-fluid equations with ideal gas equation of state.
*/
void gkyl_gr_twofluid_free(const struct gkyl_ref_count *ref);
