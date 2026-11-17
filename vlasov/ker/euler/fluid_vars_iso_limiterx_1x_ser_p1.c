#include <gkyl_iso_euler_kernels.h> 
#include <gkyl_wv_iso_euler_priv.h> 
GKYL_CU_DH void fluid_vars_iso_limiterx_1x_ser_p1(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
  const struct gkyl_wave_cell_geom *geom, 
  double *fluid_l, double *fluid_c, double *fluid_r)
{ 
  // limiter_fac:            Factor for relationship between cell slopes and cell average differences (by default: 1/sqrt(3)).
  // wv_eqn:                 Wave equation for computing waves for limiting characteristics.
  // geom:                   Geometry on the left (*only works with Cartesian components*).
  // fluid_l/c/r: [rho, rho ux, rho uy, rho uz], Fluid input state vector in left/center/right cells.

  const double *norm = geom->norm[0]; 
  const double *tau1 = geom->tau1[0]; 
  const double *tau2 = geom->tau2[0]; 

  double *rho_l = &fluid_l[0]; 
  double *rhoux_l = &fluid_l[2]; 
  double *rhouy_l = &fluid_l[4]; 
  double *rhouz_l = &fluid_l[6]; 

  double *rho_c = &fluid_c[0]; 
  double *rhoux_c = &fluid_c[2]; 
  double *rhouy_c = &fluid_c[4]; 
  double *rhouz_c = &fluid_c[6]; 

  double *rho_r = &fluid_r[0]; 
  double *rhoux_r = &fluid_r[2]; 
  double *rhouy_r = &fluid_r[4]; 
  double *rhouz_r = &fluid_r[6]; 

  double q_avg_l[4] = {0.0}; 
  q_avg_l[0] = rho_l[0]; 
  q_avg_l[1] = rhoux_l[0]; 
  q_avg_l[2] = rhouy_l[0]; 
  q_avg_l[3] = rhouz_l[0]; 
  double q_avg_c[4] = {0.0}; 
  q_avg_c[0] = rho_c[0]; 
  q_avg_c[1] = rhoux_c[0]; 
  q_avg_c[2] = rhouy_c[0]; 
  q_avg_c[3] = rhouz_c[0]; 
  double q_avg_r[4] = {0.0}; 
  q_avg_r[0] = rho_r[0]; 
  q_avg_r[1] = rhoux_r[0]; 
  q_avg_r[2] = rhouy_r[0]; 
  q_avg_r[3] = rhouz_r[0]; 
  double q_l_local[4] = {0.0}; 
  double q_c_local[4] = {0.0}; 
  double q_r_local[4] = {0.0}; 
  iso_euler_rot_to_local(wv_eqn, tau1, tau2, norm, q_avg_l, q_l_local); 
  iso_euler_rot_to_local(wv_eqn, tau1, tau2, norm, q_avg_c, q_c_local); 
  iso_euler_rot_to_local(wv_eqn, tau1, tau2, norm, q_avg_r, q_r_local); 

  double delta_l[4] = {0.0}; 
  double delta_c[4] = {0.0}; 
  double delta_r[4] = {0.0}; 
  delta_l[0] = limiter_fac*(q_c_local[0] - q_l_local[0]); 
  delta_l[1] = limiter_fac*(q_c_local[1] - q_l_local[1]); 
  delta_l[2] = limiter_fac*(q_c_local[2] - q_l_local[2]); 
  delta_l[3] = limiter_fac*(q_c_local[3] - q_l_local[3]); 
  delta_c[0] = rho_c[1]; 
  delta_c[1] = rhoux_c[1]; 
  delta_c[2] = rhouy_c[1]; 
  delta_c[3] = rhouz_c[1]; 
  delta_r[0] = limiter_fac*(q_r_local[0] - q_c_local[0]); 
  delta_r[1] = limiter_fac*(q_r_local[1] - q_c_local[1]); 
  delta_r[2] = limiter_fac*(q_r_local[2] - q_c_local[2]); 
  delta_r[3] = limiter_fac*(q_r_local[3] - q_c_local[3]); 

  double delta_c_local[4] = {0.0}; 
  iso_euler_rot_to_local(wv_eqn, tau1, tau2, norm, delta_c, delta_c_local); 

  double waves_slope_l[12] = {0.0}; 
  double waves_slope_c[12] = {0.0}; 
  double waves_slope_r[12] = {0.0}; 
  // Waves of the isothermal flux Jacobian at the cell average: the two sound waves 
  // (u -+ vt) and the shear wave (u). 
  const double vt = container_of(wv_eqn, struct wv_iso_euler, eqn)->vt; 
  const double u_n = q_c_local[1]/q_c_local[0]; 
  const double u_t1 = q_c_local[2]/q_c_local[0]; 
  const double u_t2 = q_c_local[3]/q_c_local[0]; 
  const double am_l = 0.5*(delta_l[0]*(vt + u_n) - delta_l[1])/vt; 
  const double ap_l = 0.5*(delta_l[0]*(vt - u_n) + delta_l[1])/vt; 
  waves_slope_l[0] = am_l; 
  waves_slope_l[1] = am_l*(u_n - vt); 
  waves_slope_l[2] = am_l*u_t1; 
  waves_slope_l[3] = am_l*u_t2; 
  waves_slope_l[6] = delta_l[2] - delta_l[0]*u_t1; 
  waves_slope_l[7] = delta_l[3] - delta_l[0]*u_t2; 
  waves_slope_l[8] = ap_l; 
  waves_slope_l[9] = ap_l*(u_n + vt); 
  waves_slope_l[10] = ap_l*u_t1; 
  waves_slope_l[11] = ap_l*u_t2; 
  const double am_c = 0.5*(delta_c_local[0]*(vt + u_n) - delta_c_local[1])/vt; 
  const double ap_c = 0.5*(delta_c_local[0]*(vt - u_n) + delta_c_local[1])/vt; 
  waves_slope_c[0] = am_c; 
  waves_slope_c[1] = am_c*(u_n - vt); 
  waves_slope_c[2] = am_c*u_t1; 
  waves_slope_c[3] = am_c*u_t2; 
  waves_slope_c[6] = delta_c_local[2] - delta_c_local[0]*u_t1; 
  waves_slope_c[7] = delta_c_local[3] - delta_c_local[0]*u_t2; 
  waves_slope_c[8] = ap_c; 
  waves_slope_c[9] = ap_c*(u_n + vt); 
  waves_slope_c[10] = ap_c*u_t1; 
  waves_slope_c[11] = ap_c*u_t2; 
  const double am_r = 0.5*(delta_r[0]*(vt + u_n) - delta_r[1])/vt; 
  const double ap_r = 0.5*(delta_r[0]*(vt - u_n) + delta_r[1])/vt; 
  waves_slope_r[0] = am_r; 
  waves_slope_r[1] = am_r*(u_n - vt); 
  waves_slope_r[2] = am_r*u_t1; 
  waves_slope_r[3] = am_r*u_t2; 
  waves_slope_r[6] = delta_r[2] - delta_r[0]*u_t1; 
  waves_slope_r[7] = delta_r[3] - delta_r[0]*u_t2; 
  waves_slope_r[8] = ap_r; 
  waves_slope_r[9] = ap_r*(u_n + vt); 
  waves_slope_r[10] = ap_r*u_t1; 
  waves_slope_r[11] = ap_r*u_t2; 

  double mm[12] = {0.0}; 
  double slope[4] = {0.0}; 
  for (int i = 0; i < 3; ++i) { 
    mm[4*i] = gkyl_minmod(waves_slope_c[4*i], waves_slope_l[4*i], waves_slope_r[4*i]); 
    mm[4*i+1] = gkyl_minmod(waves_slope_c[4*i+1], waves_slope_l[4*i+1], waves_slope_r[4*i+1]); 
    mm[4*i+2] = gkyl_minmod(waves_slope_c[4*i+2], waves_slope_l[4*i+2], waves_slope_r[4*i+2]); 
    mm[4*i+3] = gkyl_minmod(waves_slope_c[4*i+3], waves_slope_l[4*i+3], waves_slope_r[4*i+3]); 
    slope[0] += mm[4*i]; 
    slope[1] += mm[4*i+1]; 
    slope[2] += mm[4*i+2]; 
    slope[3] += mm[4*i+3]; 
  } 

  // Rotate limited slope back to global coordinates 
  rho_c[1] = slope[0]; 
  rhoux_c[1] = slope[1]*norm[0] + slope[2]*tau1[0] + slope[3]*tau2[0]; 
  rhouy_c[1] = slope[1]*norm[1] + slope[2]*tau1[1] + slope[3]*tau2[1]; 
  rhouz_c[1] = slope[1]*norm[2] + slope[2]*tau1[2] + slope[3]*tau2[2]; 
} 
