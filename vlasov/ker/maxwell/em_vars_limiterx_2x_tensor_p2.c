#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH void em_vars_limiterx_2x_tensor_p2(double limiter_fac, const struct gkyl_wv_eqn *wv_eqn, 
  const struct gkyl_wave_cell_geom *geom, 
  double *em_l, double *em_c, double *em_r)
{ 
  // limiter_fac:            Factor for relationship between cell slopes and cell average differences (by default: 1/sqrt(3)).
  // wv_eqn:                 Wave equation for computing waves for limiting characteristics.
  // geom:                   Geometry on the left (*only works with Cartesian components*).
  // em_l/c/r: [Ex, Ey, Ez, Bx, By, Bz, phi, psi], em input state vector in left/center/right cells.

  const double *norm = geom->norm[0]; 
  const double *tau1 = geom->tau1[0]; 
  const double *tau2 = geom->tau2[0]; 

  double *ex_l = &em_l[0]; 
  double *ey_l = &em_l[9]; 
  double *ez_l = &em_l[18]; 
  double *bx_l = &em_l[27]; 
  double *by_l = &em_l[36]; 
  double *bz_l = &em_l[45]; 
  double *phi_l = &em_l[54]; 
  double *psi_l = &em_l[63]; 

  double *ex_c = &em_c[0]; 
  double *ey_c = &em_c[9]; 
  double *ez_c = &em_c[18]; 
  double *bx_c = &em_c[27]; 
  double *by_c = &em_c[36]; 
  double *bz_c = &em_c[45]; 
  double *phi_c = &em_c[54]; 
  double *psi_c = &em_c[63]; 

  double *ex_r = &em_r[0]; 
  double *ey_r = &em_r[9]; 
  double *ez_r = &em_r[18]; 
  double *bx_r = &em_r[27]; 
  double *by_r = &em_r[36]; 
  double *bz_r = &em_r[45]; 
  double *phi_r = &em_r[54]; 
  double *psi_r = &em_r[63]; 

  double q_avg_l[8] = {0.0}; 
  q_avg_l[0] = ex_l[0]; 
  q_avg_l[1] = ey_l[0]; 
  q_avg_l[2] = ez_l[0]; 
  q_avg_l[3] = bx_l[0]; 
  q_avg_l[4] = by_l[0]; 
  q_avg_l[5] = bz_l[0]; 
  q_avg_l[6] = phi_l[0]; 
  q_avg_l[7] = psi_l[0]; 
  double q_avg_c[8] = {0.0}; 
  q_avg_c[0] = ex_c[0]; 
  q_avg_c[1] = ey_c[0]; 
  q_avg_c[2] = ez_c[0]; 
  q_avg_c[3] = bx_c[0]; 
  q_avg_c[4] = by_c[0]; 
  q_avg_c[5] = bz_c[0]; 
  q_avg_c[6] = phi_c[0]; 
  q_avg_c[7] = psi_c[0]; 
  double q_avg_r[8] = {0.0}; 
  q_avg_r[0] = ex_r[0]; 
  q_avg_r[1] = ey_r[0]; 
  q_avg_r[2] = ez_r[0]; 
  q_avg_r[3] = bx_r[0]; 
  q_avg_r[4] = by_r[0]; 
  q_avg_r[5] = bz_r[0]; 
  q_avg_r[6] = phi_r[0]; 
  q_avg_r[7] = psi_r[0]; 
  double q_l_local[8] = {0.0}; 
  double q_c_local[8] = {0.0}; 
  double q_r_local[8] = {0.0}; 
  rot_to_local(wv_eqn, tau1, tau2, norm, q_avg_l, q_l_local); 
  rot_to_local(wv_eqn, tau1, tau2, norm, q_avg_c, q_c_local); 
  rot_to_local(wv_eqn, tau1, tau2, norm, q_avg_r, q_r_local); 

  double delta_l[8] = {0.0}; 
  double delta_c[8] = {0.0}; 
  double delta_r[8] = {0.0}; 
  delta_l[0] = limiter_fac*(q_c_local[0] - q_l_local[0]); 
  delta_l[1] = limiter_fac*(q_c_local[1] - q_l_local[1]); 
  delta_l[2] = limiter_fac*(q_c_local[2] - q_l_local[2]); 
  delta_l[3] = limiter_fac*(q_c_local[3] - q_l_local[3]); 
  delta_l[4] = limiter_fac*(q_c_local[4] - q_l_local[4]); 
  delta_l[5] = limiter_fac*(q_c_local[5] - q_l_local[5]); 
  delta_l[6] = limiter_fac*(q_c_local[6] - q_l_local[6]); 
  delta_l[7] = limiter_fac*(q_c_local[7] - q_l_local[7]); 
  delta_c[0] = ex_c[1]; 
  delta_c[1] = ey_c[1]; 
  delta_c[2] = ez_c[1]; 
  delta_c[3] = bx_c[1]; 
  delta_c[4] = by_c[1]; 
  delta_c[5] = bz_c[1]; 
  delta_c[6] = phi_c[1]; 
  delta_c[7] = psi_c[1]; 
  delta_r[0] = limiter_fac*(q_r_local[0] - q_c_local[0]); 
  delta_r[1] = limiter_fac*(q_r_local[1] - q_c_local[1]); 
  delta_r[2] = limiter_fac*(q_r_local[2] - q_c_local[2]); 
  delta_r[3] = limiter_fac*(q_r_local[3] - q_c_local[3]); 
  delta_r[4] = limiter_fac*(q_r_local[4] - q_c_local[4]); 
  delta_r[5] = limiter_fac*(q_r_local[5] - q_c_local[5]); 
  delta_r[6] = limiter_fac*(q_r_local[6] - q_c_local[6]); 
  delta_r[7] = limiter_fac*(q_r_local[7] - q_c_local[7]); 

  double delta_c_local[8] = {0.0}; 
  rot_to_local(wv_eqn, tau1, tau2, norm, delta_c, delta_c_local); 

  double waves_slope_l[48] = {0.0}; 
  double waves_slope_c[48] = {0.0}; 
  double waves_slope_r[48] = {0.0}; 
  double speeds[6] = {0.0}; 
  double my_max_speed_l = wave(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_l, q_c_local, q_c_local, 1.0, 1.0, waves_slope_l, speeds); 
  double my_max_speed_c = wave(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_c_local, q_c_local, q_c_local, 1.0, 1.0, waves_slope_c, speeds); 
  double my_max_speed_r = wave(wv_eqn, GKYL_WV_HIGH_ORDER_FLUX, delta_r, q_c_local, q_c_local, 1.0, 1.0, waves_slope_r, speeds); 

  double mm[48] = {0.0}; 
  double slope[8] = {0.0}; 
  for (int i = 0; i < 6; ++i) { 
    mm[8*i] = gkyl_minmod(waves_slope_c[8*i], waves_slope_l[8*i], waves_slope_r[8*i]); 
    mm[8*i+1] = gkyl_minmod(waves_slope_c[8*i+1], waves_slope_l[8*i+1], waves_slope_r[8*i+1]); 
    mm[8*i+2] = gkyl_minmod(waves_slope_c[8*i+2], waves_slope_l[8*i+2], waves_slope_r[8*i+2]); 
    mm[8*i+3] = gkyl_minmod(waves_slope_c[8*i+3], waves_slope_l[8*i+3], waves_slope_r[8*i+3]); 
    mm[8*i+4] = gkyl_minmod(waves_slope_c[8*i+4], waves_slope_l[8*i+4], waves_slope_r[8*i+4]); 
    mm[8*i+5] = gkyl_minmod(waves_slope_c[8*i+5], waves_slope_l[8*i+5], waves_slope_r[8*i+5]); 
    mm[8*i+6] = gkyl_minmod(waves_slope_c[8*i+6], waves_slope_l[8*i+6], waves_slope_r[8*i+6]); 
    mm[8*i+7] = gkyl_minmod(waves_slope_c[8*i+7], waves_slope_l[8*i+7], waves_slope_r[8*i+7]); 
    slope[0] += mm[8*i]; 
    slope[1] += mm[8*i+1]; 
    slope[2] += mm[8*i+2]; 
    slope[3] += mm[8*i+3]; 
    slope[4] += mm[8*i+4]; 
    slope[5] += mm[8*i+5]; 
    slope[6] += mm[8*i+6]; 
    slope[7] += mm[8*i+7]; 
  } 

  // Rotate limited slope back to global coordinates 
  ex_c[1] = slope[0]*norm[0] + slope[1]*tau1[0] + slope[2]*tau2[0]; 
  ey_c[1] = slope[0]*norm[1] + slope[1]*tau1[1] + slope[2]*tau2[1]; 
  ez_c[1] = slope[0]*norm[2] + slope[1]*tau1[2] + slope[2]*tau2[2]; 
  bx_c[1] = slope[3]*norm[0] + slope[4]*tau1[0] + slope[5]*tau2[0]; 
  by_c[1] = slope[3]*norm[1] + slope[4]*tau1[1] + slope[5]*tau2[1]; 
  bz_c[1] = slope[3]*norm[2] + slope[4]*tau1[2] + slope[5]*tau2[2]; 
  phi_c[1] = slope[6]; 
  psi_c[1] = slope[7]; 
  // Components (in the local frame) whose slope the limiter changed by more than round-off: 
  // the change summed over the waves against the size of the component and of its slopes in 
  // the three cells. 
  int is_limited[8] = {0}; 
  for (int k = 0; k < 8; ++k) { 
    double slope_change = 0.0; 
    double slope_size = fabs(q_l_local[k]) + fabs(q_c_local[k]) + fabs(q_r_local[k]); 
    for (int i = 0; i < 6; ++i) { 
      slope_change += fabs(mm[8*i+k] - waves_slope_c[8*i+k]); 
      slope_size += fabs(waves_slope_c[8*i+k]) + fabs(waves_slope_l[8*i+k]) + fabs(waves_slope_r[8*i+k]); 
    } 
    is_limited[k] = slope_change > 1.0e-12*slope_size; 
  } 
  // A limited local component of E or B drops the higher coefficients of the global 
  // components it projects onto. 
  int is_limited_ex = (is_limited[0] && norm[0] != 0.0) || (is_limited[1] && tau1[0] != 0.0) || (is_limited[2] && tau2[0] != 0.0); 
  int is_limited_ey = (is_limited[0] && norm[1] != 0.0) || (is_limited[1] && tau1[1] != 0.0) || (is_limited[2] && tau2[1] != 0.0); 
  int is_limited_ez = (is_limited[0] && norm[2] != 0.0) || (is_limited[1] && tau1[2] != 0.0) || (is_limited[2] && tau2[2] != 0.0); 
  int is_limited_bx = (is_limited[3] && norm[0] != 0.0) || (is_limited[4] && tau1[0] != 0.0) || (is_limited[5] && tau2[0] != 0.0); 
  int is_limited_by = (is_limited[3] && norm[1] != 0.0) || (is_limited[4] && tau1[1] != 0.0) || (is_limited[5] && tau2[1] != 0.0); 
  int is_limited_bz = (is_limited[3] && norm[2] != 0.0) || (is_limited[4] && tau1[2] != 0.0) || (is_limited[5] && tau2[2] != 0.0); 
  if (is_limited_ex) { 
    ex_c[3] = 0.0; 
    ex_c[4] = 0.0; 
    ex_c[5] = 0.0; 
    ex_c[6] = 0.0; 
    ex_c[7] = 0.0; 
    ex_c[8] = 0.0; 
  } 
  if (is_limited_ey) { 
    ey_c[3] = 0.0; 
    ey_c[4] = 0.0; 
    ey_c[5] = 0.0; 
    ey_c[6] = 0.0; 
    ey_c[7] = 0.0; 
    ey_c[8] = 0.0; 
  } 
  if (is_limited_ez) { 
    ez_c[3] = 0.0; 
    ez_c[4] = 0.0; 
    ez_c[5] = 0.0; 
    ez_c[6] = 0.0; 
    ez_c[7] = 0.0; 
    ez_c[8] = 0.0; 
  } 
  if (is_limited_bx) { 
    bx_c[3] = 0.0; 
    bx_c[4] = 0.0; 
    bx_c[5] = 0.0; 
    bx_c[6] = 0.0; 
    bx_c[7] = 0.0; 
    bx_c[8] = 0.0; 
  } 
  if (is_limited_by) { 
    by_c[3] = 0.0; 
    by_c[4] = 0.0; 
    by_c[5] = 0.0; 
    by_c[6] = 0.0; 
    by_c[7] = 0.0; 
    by_c[8] = 0.0; 
  } 
  if (is_limited_bz) { 
    bz_c[3] = 0.0; 
    bz_c[4] = 0.0; 
    bz_c[5] = 0.0; 
    bz_c[6] = 0.0; 
    bz_c[7] = 0.0; 
    bz_c[8] = 0.0; 
  } 
  if (is_limited[6]) { 
    phi_c[3] = 0.0; 
    phi_c[4] = 0.0; 
    phi_c[5] = 0.0; 
    phi_c[6] = 0.0; 
    phi_c[7] = 0.0; 
    phi_c[8] = 0.0; 
  } 
  if (is_limited[7]) { 
    psi_c[3] = 0.0; 
    psi_c[4] = 0.0; 
    psi_c[5] = 0.0; 
    psi_c[6] = 0.0; 
    psi_c[7] = 0.0; 
    psi_c[8] = 0.0; 
  } 
} 
