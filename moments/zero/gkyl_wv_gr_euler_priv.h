#pragma once

// Private header, not for direct use in user-facing code.

#include <math.h>
#include <gkyl_array.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gr_minkowski.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

struct gkyl_tov; // Frozen TOV equilibrium table, for optional well-balancing.

struct wv_gr_euler {
  struct gkyl_wv_eqn eqn; // Base equation object.
  struct gkyl_gr_spacetime *spacetime; // Pointer to base spacetime object.
  double gas_gamma; // Adiabatic index.

  enum gkyl_spacetime_gauge spacetime_gauge; // Spacetime gauge choice.
  int reinit_freq; // Spacetime reinitialization frequency.

  // Optional well-balancing against a static TOV equilibrium. When active, the flux waves and the geometric source subtract a frozen discrete equilibrium: the t=0 cell-averaged conserved state,
  // projected once and carried in the spare slots q[71]=<sqrt(g) D_eq>, q[72]=<sqrt(g) Etot_eq> (Etot is the Valencia energy variable in slot [4], = E - D; the driver calls it Etot).
  // tov_eq below is the frozen TOV table, but here it is used only as the on/off switch (non-NULL => well-balance). The equilibrium values are not re-looked-up in the table at runtime; they are the discrete cell-averages in q[71],q[72]. 
  const struct gkyl_tov *tov_eq; // NULL -> plain GR Euler, no WB. Used as a presence flag only.

  bool disable_well_balanced; // tov_eq set but WB subtraction switched off everywhere (collapse stage).

  double rho_atm, p_atm; // Atmosphere floors used by primitive recovery.
};

/**
* Compute flux vector. Assumes rotation to local coordinate system.
*
* @param gas_gamma Adiabatic index.
* @param q Conserved variable vector.
* @param flux Flux vector in direction 'dir' (output).
*/
GKYL_CU_D
void
gkyl_gr_euler_flux(double gas_gamma, const double q[71], double flux[71]);

/**
* Build the well-balancing equilibrium reference state q_eq (family or frozen, per the eqn config) at the cell's live metric. Used by the MP scheme's well-balanced reconstruction so the reconstructed face states share the same equilibrium the flux/source subtract.
*
* @param eqn Base equation object (must be GR Euler).
* @param q Conserved variable vector (73 components, incl. metric and WB slots).
* @param q_eq Equilibrium conserved + metric state (output, 0..70).
*/
void
gkyl_gr_euler_equilibrium(const struct gkyl_wv_eqn *eqn, const double q[76], double q_eq[71]);

/**
* Whether the WB subtraction is switched off for this GR-Euler equation (collapse stage). Single source of truth for the three WB sites (reconstruction, flux jump, MoL source).
*/
bool
gkyl_gr_euler_wb_disabled(const struct gkyl_wv_eqn *eqn);

/**
* Runtime toggle of the frozen-equilibrium well-balancing (mid-run dynamical-collapse trigger).
*/
void
gkyl_gr_euler_set_wb_disabled(const struct gkyl_wv_eqn *eqn, bool disabled);


/**
* Specific enthalpy h = 1 + eps + p/rho from (rho, p) for the gamma-law EOS.
*/
double
gkyl_gr_euler_specific_enthalpy(double gas_gamma, double rho, double p);

/**
* Compute primitive variables given the conserved variables.
*
* @param gas_gamma Adiabatic index.
* @param q Conserved variable vector.
* @param v Primitive variable vector (output).
*/
GKYL_CU_D
void
gkyl_gr_euler_prim_vars(double gas_gamma, const double q[71], double v[71]);

/**
* Compute inverse spatial metric tensor (in covariant component form) given the conserved variables.
*
* @param q Conserved variable vector.
* @param inv_spatial_metric Inverse spatial metric tensor (output).
*/
GKYL_CU_D
void
gkyl_gr_euler_inv_spatial_metric(const double q[71], double ***inv_spatial_metric);

/**
* Compute perfect fluid stress-energy tensor (in contravariant component form) given the conserved variables.
*
* @param gas_gamma Adiabatic index.
* @param q Conserved variable vector.
* @param stress_energy Stress-energy tensor (output).
*/
GKYL_CU_D
void
gkyl_gr_euler_stress_energy_tensor(double gas_gamma, const double q[71], double ***stress_energy);

/**
* Compute maximum absolute wave speed.
*
* @param gas_gamma Adiabatic index.
* @param q Conserved variable vector.
* @return Maximum absolute wave speed for a given q.
*/
GKYL_CU_D
static inline double
gkyl_gr_euler_max_abs_speed(double gas_gamma, const double q[71]);

/**
* Compute Riemann variables given the conserved variables.
*
* @param eqn Base equation object.
* @param qstate Current state vector.
* @param qin Conserved variable vector (input).
* @param wout Riemann variable vector (output).
*/
GKYL_CU_D
static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout);

/**
* Compute conserved variables given the Riemann variables.
*
* @param eqn Base equation object.
* @param qstate Current state vector.
* @param win Riemann variable vector (input).
* @param qout Conserved variable vector (output).
*/
GKYL_CU_D
static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double *qout);

/**
* Boundary condition function for applying wall boundary conditions for the general relativistic Euler equations with ideal gas equation of state.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply wall boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
GKYL_CU_D
static void
gr_euler_wall(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx);

/**
* Boundary condition function for applying no-slip boundary conditions for the general relativistic Euler equations with ideal gas equation of state.
*
* @param eqn Base equation object.
* @param t Current simulation time.
* @param nc Number of boundary cells to which to apply no-slip boundary conditions.
* @param skin Skin cells in boundary region (from which values are copied).
* @param ghost Ghost cells in boundary region (to which values are copied).
* @param ctx Context to pass to the function.
*/
GKYL_CU_D
static void
gr_euler_no_slip(const struct gkyl_wv_eqn* eqn, double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx);

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
GKYL_CU_D
static inline void
rot_to_local(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal,
  double* GKYL_RESTRICT qlocal);

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
GKYL_CU_D
static inline void
rot_to_global(const struct gkyl_wv_eqn* eqn, const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal,
  double* GKYL_RESTRICT qglobal);

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
GKYL_CU_D
static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s);

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
GKYL_CU_D
static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq);

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
GKYL_CU_D
static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s);

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
GKYL_CU_D
static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq);

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
GKYL_CU_D
static double
wave_roe(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s);

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
GKYL_CU_D
static void
qfluct_roe(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq);

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
GKYL_CU_D
static double
wave_roe_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s);

/**
* Compute fluctuations using Roe fluxes (with potential fallback).
*
* @param eqn Base equation object.
* @param type Type of Riemann-solver flux to use.
* @param ql Conserved variable vector on the left of the interface.
* @param qr Conserved variable vector on the right of the interface.
* @param waves Waves (input).
* @param s Wave speeds (input).
* @param amdq Left-moving fluctuations (output).
* @param qpdq Right-moving fluctuations (output).
*/
GKYL_CU_D
static void
qfluct_roe_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq);

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
GKYL_CU_D
static double
wave_hll(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s);

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
GKYL_CU_D
static void
qfluct_hll(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq);

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
GKYL_CU_D
static double
wave_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, const double phil, const double phir, double* waves, double* s);

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
GKYL_CU_D
static void
qfluct_hll_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double phil, const double phir, const double* waves, const double* s,
  double* amdq, double* apdq);

/**
* Compute jump in flux given two conserved variable states.
*
* @param eqn Base equation object.
* @param ql Conserved variable vector on the left of the interface (input).
* @param qr Conserved variable vector on the right of the interface (input).
* @param flux_jump Jump in flux vector (output).
* @return Maximum wave speeds for states ql and qr.
*/
GKYL_CU_D
static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump);

/**
* Determine whether invariant domain of the general relativistic Euler equations with ideal gas equation of state is satisfied.
*
* @param eqn Base equation object.
* @param q Conserved variable vector.
* @return Whether the invariant domain is satisfied.
*/
GKYL_CU_D
static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q);

/**
* Compute maximum wave speed from a conserved variable vector.
*
* @param eqn Base equation object.
* @param q Conserved variable vector.
* @return Maximum absolute wave speed.
*/
GKYL_CU_D
static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q);

/**
* Convert conserved variables to diagnostic variables.
*
* @param eqn Base equation object.
* @param qin Conserved variable vector (input).
* @param diag Diagnostic variable vector (output).
*/
GKYL_CU_D
static inline void
gr_euler_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag);

/**
* Compute forcing/source term vector from conserved variables.
*
* @param eqn Base equation object.
* @param qin Conserved variable vector (input).
* @param sout Forcing/source term vector (output).
*/
GKYL_CU_DH
static inline void
gr_euler_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout);

/**
* Free general relativistic Euler equations object with ideal gas equation of state.
*
* @param ref Reference counter for general relativistic Euler equations with ideal gas equation of state.
*/
void gkyl_gr_euler_free(const struct gkyl_ref_count* ref);
