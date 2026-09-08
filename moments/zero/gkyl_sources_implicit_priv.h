#pragma once

// Forward-declaration of the private gkyl_moment_em_coupling object type.
typedef struct gkyl_moment_em_coupling gkyl_moment_em_coupling;

/**
* Rotate the pressure tensor for the 10-moment equations with respect
* to the magnetic field.
*
* @param q_over_m Ratio of species charge to species mass.
* @param dt Current stable time-step.
* @param em Array of electromagnetic variables.
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param p_tensor_old Old components of the pressure tensor (before rotation).
* @param p_tensor_rhs Source terms appearing on the right-hand-side of the pressure tensor equations.
* @param p_tensor_new New components of the pressure tensor (after rotation).
*/
void pressure_tensor_rotate(double q_over_m, double dt, const double *em, const double *ext_em,
  double p_tensor_old[6], double p_tensor_rhs[6], double p_tensor_new[6]);

/**
* Integrate the electromagnetic source terms of a charged multi-fluid equation system within a single cell, using an implicit forcing solver
* (specifically the time-centered Crank-Nicolson/implicit Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_rhs_s Array of input (rhs) fluid variables (array size = nfluids x 4). 
* @param fluid_s Array of output fluid variables (array size = nfluids).
* @param app_accel_s Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
*/
void implicit_em_source_update(const gkyl_moment_em_coupling *mom_em, double t_curr, double dt,
  double fluid_rhs_s[GKYL_MAX_SPECIES][4], double *fluid_s[GKYL_MAX_SPECIES],
  const double *app_accel_s[GKYL_MAX_SPECIES], double *em, const double *app_current,
  const double *ext_em);

/**
* Integrate the momentum source terms of a neutral multi-fluid equation system within a single cell, using an implicit forcing solver (specifically
* the time-centered Crank-Nicolson/implicit Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_rhs_s Array of input (rhs) fluid variables (array size = nfluids x 4).  
* @param fluid_s Array of output fluid variables (array size = nfluids).
* @param app_accel_s Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param p_rhs_s Array of right-hand-side source terms; for neutral ten-moment fluids the pressure-tensor
*   components (indices 4-9) are applied as a direct forward-Euler increment (e.g. a heat-flux closure).
*/
void implicit_neut_source_update(const gkyl_moment_em_coupling *mom_em, double t_curr, double dt,
  double fluid_rhs_s[GKYL_MAX_SPECIES][4], double *fluid_s[GKYL_MAX_SPECIES],
  const double *app_accel_s[GKYL_MAX_SPECIES], const double *p_rhs_s[GKYL_MAX_SPECIES]);

/**
* Integrate the collisional source terms of a multi-fluid equation system within a single cell, using an implicit forcing solver (specifically the
* time-centered Crank-Nicolson/implicit Runge-Kutta method, with a direct matrix inversion).
*
* @param mom_em Moment-EM coupling object.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
*/
void implicit_collision_source_update(
  const gkyl_moment_em_coupling *mom_em, double dt, double *fluid_s[GKYL_MAX_SPECIES]);

/**
* Integrate the frictional source terms in the multi-fluid equation system within a single cell, using an implicit forcing solver (specifically
* the time-centered Crank-Nicolson/implicit Runge-Kutta method, with a direct matrix inversion), over half a stable time-step.
*
* @param mom_em Moment-EM coupling object.
* @param Z Ionization number.
* @param T_elc Electron temperature.
* @param Lambda_ee Electron-electron collisional term.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param f_elc_old Array of old electron fluid variables (before source update).
* @param f_ion_old Array of old ion fluid variables (before source update).
* @param f_elc_new Array of new electron fluid variables (after source update).
* @param f_ion_new Array of new ion fluid variables (after source update).
* @param app_accel_s Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param em_old Array of old electromagnetic variables (before source update).
* @param em_new Array of new electromagnetic variables (after source update).
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
*/
void implicit_frictional_source_update_half(const gkyl_moment_em_coupling *mom_em, const double Z,
  const double T_elc, const double Lambda_ee, double t_curr, const double dt, double *f_elc_old,
  double *f_ion_old, double *f_elc_new, double *f_ion_new,
  const double *app_accel_s[GKYL_MAX_SPECIES], double *em_old, double *em_new,
  const double *app_current, const double *ext_em);

/**
* Integrate the frictional source terms in the multi-fluid equation system within a single cell, using an implicit forcing solver (specifically
* the time-centered Crank-Nicolson/implicit Runge-Kutta method, with a direct matrix inversion).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid varaibles (array size = nfluids).
* @param app_accel_s Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param em Array of electromagnetic variables
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
*/
void implicit_frictional_source_update(const gkyl_moment_em_coupling *mom_em, double t_curr,
  double dt, double *fluid_s[GKYL_MAX_SPECIES], const double *app_accel_s[GKYL_MAX_SPECIES],
  double *em, const double *app_current, const double *ext_em);

/**
* Integrate the electromagnetic source terms in the multi-fluid equation system within each cell, using an implicit forcing solver (specifically
* the time-centered Crank-Nicolson/implicit Runge-Kutta method).
*
* @param mom_em Moment-EM coupling object.
* @param t_curr Current simulation time.
* @param dt Current stable time-step.
* @param fluid_s Array of fluid variables (array size = nfluids).
* @param app_accel_s Array of acceleration terms to be applied to the fluid equations (for external forces).
* @param p_rhs_s Array of RHS/source terms to be applied to the pressure tensor (for the case of 10-moment gradient-based closure only).
* @param em Array of electromagnetic variables.
* @param app_current Array of current terms to be applied to the fluid equations (for external current driving).
* @param ext_em External electromagnetic variables (for EM fields coming from external sources, e.g. coils, capacitors, etc.).
* @param nT_sources_s Array of number density and temperature source terms.
*/
void implicit_source_coupling_update(const gkyl_moment_em_coupling *mom_em, double t_curr,
  double dt, double *fluid_s[GKYL_MAX_SPECIES], const double *app_accel_s[GKYL_MAX_SPECIES],
  const double *p_rhs_s[GKYL_MAX_SPECIES], double *em, const double *app_current,
  const double *ext_em, const double *nT_sources_s[GKYL_MAX_SPECIES]);
