#pragma once 
#include <math.h> 
#include <gkyl_wave_geom.h> 
#include <gkyl_wv_eqn.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

// Rotate the state [rho, rho ux, rho uy, rho uz] to the frame of an interface (normal, tangents). 
GKYL_CU_DH static inline void 
iso_euler_rot_to_local(const struct gkyl_wv_eqn *eqn, const double *tau1, const double *tau2, 
  const double *norm, const double *GKYL_RESTRICT qglobal, double *GKYL_RESTRICT qlocal) 
{ 
  qlocal[0] = qglobal[0]; 
  qlocal[1] = qglobal[1]*norm[0] + qglobal[2]*norm[1] + qglobal[3]*norm[2]; 
  qlocal[2] = qglobal[1]*tau1[0] + qglobal[2]*tau1[1] + qglobal[3]*tau1[2]; 
  qlocal[3] = qglobal[1]*tau2[0] + qglobal[2]*tau2[1] + qglobal[3]*tau2[2]; 
} 

// Rotate the state back to the global frame. 
GKYL_CU_DH static inline void 
iso_euler_rot_to_global(const struct gkyl_wv_eqn *eqn, const double *tau1, const double *tau2, 
  const double *norm, const double *GKYL_RESTRICT qlocal, double *GKYL_RESTRICT qglobal) 
{ 
  qglobal[0] = qlocal[0]; 
  qglobal[1] = qlocal[1]*norm[0] + qlocal[2]*tau1[0] + qlocal[3]*tau2[0]; 
  qglobal[2] = qlocal[1]*norm[1] + qlocal[2]*tau1[1] + qlocal[3]*tau2[1]; 
  qglobal[3] = qlocal[1]*norm[2] + qlocal[2]*tau1[2] + qlocal[3]*tau2[2]; 
} 

// Flux along the normal of the state in the frame of the interface. 
GKYL_CU_DH static inline void 
iso_euler_flux(double vt, const double *q, double *flux) 
{ 
  const double u = q[1]/q[0]; 
  flux[0] = q[1]; 
  flux[1] = q[1]*u + q[0]*vt*vt; 
  flux[2] = q[2]*u; 
  flux[3] = q[3]*u; 
} 

GKYL_CU_DH double iso_euler_vol_1x_ser_p1(const double *w, const double *dxv, double vt, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_pressure_1x_ser_p1(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_iso_source_1x_ser_p1(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_iso_em_source_1x_ser_p1(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH int fluid_vars_iso_positivity_1x_ser_p1(double eps, double* GKYL_RESTRICT fluid); 
GKYL_CU_DH double iso_euler_surfx_1x_ser_p1(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_limiterx_1x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 

GKYL_CU_DH double iso_euler_vol_1x_ser_p2(const double *w, const double *dxv, double vt, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_pressure_1x_ser_p2(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_iso_source_1x_ser_p2(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_iso_em_source_1x_ser_p2(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double iso_euler_surfx_1x_ser_p2(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_limiterx_1x_ser_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 

GKYL_CU_DH double iso_euler_vol_1x_ser_p3(const double *w, const double *dxv, double vt, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_pressure_1x_ser_p3(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_iso_source_1x_ser_p3(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_iso_em_source_1x_ser_p3(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double iso_euler_surfx_1x_ser_p3(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_limiterx_1x_ser_p3(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 

GKYL_CU_DH double iso_euler_vol_2x_ser_p1(const double *w, const double *dxv, double vt, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_pressure_2x_ser_p1(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_iso_source_2x_ser_p1(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_iso_em_source_2x_ser_p1(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH int fluid_vars_iso_positivity_2x_ser_p1(double eps, double* GKYL_RESTRICT fluid); 
GKYL_CU_DH double iso_euler_surfx_2x_ser_p1(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_limiterx_2x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double iso_euler_surfy_2x_ser_p1(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_limitery_2x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 

GKYL_CU_DH double iso_euler_vol_3x_ser_p1(const double *w, const double *dxv, double vt, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_pressure_3x_ser_p1(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_iso_source_3x_ser_p1(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_iso_em_source_3x_ser_p1(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH int fluid_vars_iso_positivity_3x_ser_p1(double eps, double* GKYL_RESTRICT fluid); 
GKYL_CU_DH double iso_euler_surfx_3x_ser_p1(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_limiterx_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double iso_euler_surfy_3x_ser_p1(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_limitery_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double iso_euler_surfz_3x_ser_p1(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_limiterz_3x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 

GKYL_CU_DH double iso_euler_vol_2x_tensor_p2(const double *w, const double *dxv, double vt, const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_pressure_2x_tensor_p2(double param, const double *fluid, const double *u, double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf); 
GKYL_CU_DH void fluid_vars_iso_source_2x_tensor_p2(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double fluid_vars_iso_em_source_2x_tensor_p2(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out); 
GKYL_CU_DH double iso_euler_surfx_2x_tensor_p2(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_limiterx_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 
GKYL_CU_DH double iso_euler_surfy_2x_tensor_p2(const double *w, const double *dxv, const struct gkyl_wv_eqn *wv_eqn, 
    const struct gkyl_wave_cell_geom *geom_l, const struct gkyl_wave_cell_geom *geom_r, 
    const double *fluid_l, const double *fluid_c, const double *fluid_r, double* GKYL_RESTRICT out); 
GKYL_CU_DH void fluid_vars_iso_limitery_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_cell_geom *geom, double *fluid_l, double *fluid_c, double *fluid_r); 

EXTERN_C_END 
