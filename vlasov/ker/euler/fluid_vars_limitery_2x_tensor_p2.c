#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_limitery_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
  const struct gkyl_wave_cell_geom *geom, 
  double *fluid_l, double *fluid_c, double *fluid_r)
{ 
  // limiter_fac:            Factor for relationship between cell slopes and cell average differences (by default: 1/sqrt(3)).
  // wv_eqn:                 Wave equation for computing waves for limiting characteristics.
  // geom:                   Geometry on the left (*only works with Cartesian components*).
  // fluid_l/c/r: [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector in left/center/right cells.

  const double *norm = geom->norm[1]; 
  const double *tau1 = geom->tau1[1]; 
  const double *tau2 = geom->tau2[1]; 

  double *rho_l = &fluid_l[0]; 
  double *rhoux_l = &fluid_l[9]; 
  double *rhouy_l = &fluid_l[18]; 
  double *rhouz_l = &fluid_l[27]; 
  double *energy_l = &fluid_l[36]; 

  double *rho_c = &fluid_c[0]; 
  double *rhoux_c = &fluid_c[9]; 
  double *rhouy_c = &fluid_c[18]; 
  double *rhouz_c = &fluid_c[27]; 
  double *energy_c = &fluid_c[36]; 

  double *rho_r = &fluid_r[0]; 
  double *rhoux_r = &fluid_r[9]; 
  double *rhouy_r = &fluid_r[18]; 
  double *rhouz_r = &fluid_r[27]; 
  double *energy_r = &fluid_r[36]; 

  double q_avg_l[5] = {0.0}; 
  q_avg_l[0] = rho_l[0]; 
  q_avg_l[1] = rhoux_l[0]; 
  q_avg_l[2] = rhouy_l[0]; 
  q_avg_l[3] = rhouz_l[0]; 
  q_avg_l[4] = energy_l[0]; 
  double q_avg_c[5] = {0.0}; 
  q_avg_c[0] = rho_c[0]; 
  q_avg_c[1] = rhoux_c[0]; 
  q_avg_c[2] = rhouy_c[0]; 
  q_avg_c[3] = rhouz_c[0]; 
  q_avg_c[4] = energy_c[0]; 
  double q_avg_r[5] = {0.0}; 
  q_avg_r[0] = rho_r[0]; 
  q_avg_r[1] = rhoux_r[0]; 
  q_avg_r[2] = rhouy_r[0]; 
  q_avg_r[3] = rhouz_r[0]; 
  q_avg_r[4] = energy_r[0]; 
  double q_l_local[5] = {0.0}; 
  double q_c_local[5] = {0.0}; 
  double q_r_local[5] = {0.0}; 
  rot_to_local(wv_eqn, tau1, tau2, norm, q_avg_l, q_l_local); 
  rot_to_local(wv_eqn, tau1, tau2, norm, q_avg_c, q_c_local); 
  rot_to_local(wv_eqn, tau1, tau2, norm, q_avg_r, q_r_local); 

  double delta_l[5] = {0.0}; 
  double delta_c[5] = {0.0}; 
  double delta_r[5] = {0.0}; 
  delta_l[0] = limiter_fac*(q_c_local[0] - q_l_local[0]); 
  delta_l[1] = limiter_fac*(q_c_local[1] - q_l_local[1]); 
  delta_l[2] = limiter_fac*(q_c_local[2] - q_l_local[2]); 
  delta_l[3] = limiter_fac*(q_c_local[3] - q_l_local[3]); 
  delta_l[4] = limiter_fac*(q_c_local[4] - q_l_local[4]); 
  delta_c[0] = rho_c[2]; 
  delta_c[1] = rhoux_c[2]; 
  delta_c[2] = rhouy_c[2]; 
  delta_c[3] = rhouz_c[2]; 
  delta_c[4] = energy_c[2]; 
  delta_r[0] = limiter_fac*(q_r_local[0] - q_c_local[0]); 
  delta_r[1] = limiter_fac*(q_r_local[1] - q_c_local[1]); 
  delta_r[2] = limiter_fac*(q_r_local[2] - q_c_local[2]); 
  delta_r[3] = limiter_fac*(q_r_local[3] - q_c_local[3]); 
  delta_r[4] = limiter_fac*(q_r_local[4] - q_c_local[4]); 

  double delta_c_local[5] = {0.0}; 
  rot_to_local(wv_eqn, tau1, tau2, norm, delta_c, delta_c_local); 

  double waves_slope_l[15] = {0.0}; 
  double waves_slope_c[15] = {0.0}; 
  double waves_slope_r[15] = {0.0}; 
  double speeds[3] = {0.0}; 
  double my_max_speed_l = wave_roe(wv_eqn, delta_l, q_c_local, q_c_local, waves_slope_l, speeds); 
  double my_max_speed_c = wave_roe(wv_eqn, delta_c_local, q_c_local, q_c_local, waves_slope_c, speeds); 
  double my_max_speed_r = wave_roe(wv_eqn, delta_r, q_c_local, q_c_local, waves_slope_r, speeds); 

  double mm[15] = {0.0}; 
  double slope[5] = {0.0}; 
  for (int i = 0; i < 3; ++i) { 
    mm[5*i] = gkyl_minmod(waves_slope_c[5*i], waves_slope_l[5*i], waves_slope_r[5*i]); 
    mm[5*i+1] = gkyl_minmod(waves_slope_c[5*i+1], waves_slope_l[5*i+1], waves_slope_r[5*i+1]); 
    mm[5*i+2] = gkyl_minmod(waves_slope_c[5*i+2], waves_slope_l[5*i+2], waves_slope_r[5*i+2]); 
    mm[5*i+3] = gkyl_minmod(waves_slope_c[5*i+3], waves_slope_l[5*i+3], waves_slope_r[5*i+3]); 
    mm[5*i+4] = gkyl_minmod(waves_slope_c[5*i+4], waves_slope_l[5*i+4], waves_slope_r[5*i+4]); 
    slope[0] += mm[5*i]; 
    slope[1] += mm[5*i+1]; 
    slope[2] += mm[5*i+2]; 
    slope[3] += mm[5*i+3]; 
    slope[4] += mm[5*i+4]; 
  } 

  // Rotate limited slope back to global coordinates 
  rho_c[2] = slope[0]; 
  rhoux_c[2] = slope[1]*norm[0] + slope[2]*tau1[0] + slope[3]*tau2[0]; 
  rhouy_c[2] = slope[1]*norm[1] + slope[2]*tau1[1] + slope[3]*tau2[1]; 
  rhouz_c[2] = slope[1]*norm[2] + slope[2]*tau1[2] + slope[3]*tau2[2]; 
  energy_c[2] = slope[4]; 
  // Components (in the local frame) whose slope the limiter changed by more than round-off: 
  // the change summed over the waves against the size of the component and of its slopes in 
  // the three cells. 
  int is_limited[5] = {0}; 
  for (int k = 0; k < 5; ++k) { 
    double slope_change = 0.0; 
    double slope_size = fabs(q_l_local[k]) + fabs(q_c_local[k]) + fabs(q_r_local[k]); 
    for (int i = 0; i < 3; ++i) { 
      slope_change += fabs(mm[5*i+k] - waves_slope_c[5*i+k]); 
      slope_size += fabs(waves_slope_c[5*i+k]) + fabs(waves_slope_l[5*i+k]) + fabs(waves_slope_r[5*i+k]); 
    } 
    is_limited[k] = slope_change > 1.0e-12*slope_size; 
  } 
  // A limited local momentum component drops the higher coefficients of the global 
  // momentum components it projects onto. 
  int is_limited_ux = (is_limited[1] && norm[0] != 0.0) || (is_limited[2] && tau1[0] != 0.0) || (is_limited[3] && tau2[0] != 0.0); 
  int is_limited_uy = (is_limited[1] && norm[1] != 0.0) || (is_limited[2] && tau1[1] != 0.0) || (is_limited[3] && tau2[1] != 0.0); 
  int is_limited_uz = (is_limited[1] && norm[2] != 0.0) || (is_limited[2] && tau1[2] != 0.0) || (is_limited[3] && tau2[2] != 0.0); 
  if (is_limited[0]) { 
    rho_c[3] = 0.0; 
    rho_c[4] = 0.0; 
    rho_c[5] = 0.0; 
    rho_c[6] = 0.0; 
    rho_c[7] = 0.0; 
    rho_c[8] = 0.0; 
  } 
  if (is_limited_ux) { 
    rhoux_c[3] = 0.0; 
    rhoux_c[4] = 0.0; 
    rhoux_c[5] = 0.0; 
    rhoux_c[6] = 0.0; 
    rhoux_c[7] = 0.0; 
    rhoux_c[8] = 0.0; 
  } 
  if (is_limited_uy) { 
    rhouy_c[3] = 0.0; 
    rhouy_c[4] = 0.0; 
    rhouy_c[5] = 0.0; 
    rhouy_c[6] = 0.0; 
    rhouy_c[7] = 0.0; 
    rhouy_c[8] = 0.0; 
  } 
  if (is_limited_uz) { 
    rhouz_c[3] = 0.0; 
    rhouz_c[4] = 0.0; 
    rhouz_c[5] = 0.0; 
    rhouz_c[6] = 0.0; 
    rhouz_c[7] = 0.0; 
    rhouz_c[8] = 0.0; 
  } 
  if (is_limited[4]) { 
    energy_c[3] = 0.0; 
    energy_c[4] = 0.0; 
    energy_c[5] = 0.0; 
    energy_c[6] = 0.0; 
    energy_c[7] = 0.0; 
    energy_c[8] = 0.0; 
  } 
} 
