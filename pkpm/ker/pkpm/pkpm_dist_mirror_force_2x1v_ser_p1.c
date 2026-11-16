#include <gkyl_vlasov_pkpm_kernels.h> 
GKYL_CU_DH void pkpm_dist_mirror_force_2x1v_ser_p1(const double *w, const double *dxv, 
  const double *pkpm_prim, const double *nu_prim_moms_sum, 
  const double *div_b, const double *pkpm_accel_vars, 
  const double *f, const double *F_k_p_1, 
  double* GKYL_RESTRICT g_dist_source, double* GKYL_RESTRICT F_k_m_1) 
{ 
  // w[NDIM]:          Cell-center coordinates. 
  // dxv[NDIM]:        Cell spacing. 
  // pkpm_prim:        Input primitive variables [ux, uy, uz, 1/rho div(p_par b), T_perp/m, m/T_perp]. 
  // nu_prim_moms_sum: Input sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // div_b:            Input volume expansion of div(b). 
  // pkpm_accel_vars:  Input pkpm acceleration variables [T_perp/m*div(b), bb:grad(u), p_force, p_perp_source]. 
  // f:                Input distribution functions [F_0, T_perp/m G = T_perp/m (F_0 - F_1)].
  // F_k_p_1:          Input k+1 distribution function. F_2 expansion is the first NP coefficients. 
  // g_dist_source:    Output [2.0*T_perp/m*(2.0*T_perp/m G + T_perp/m (F_2 - F_0)),  
  //                   (-vpar div(b) + bb:grad(u) - div(u) - 2 nu) T_perp/m G + 2 nu vth^2 F_0 ].
  //                   First output is mirror force source, second output is vperp characteristics source.
  // F_k_m_1:          Output k-1 distribution function. F_1 expansion is the first NP coefficients. 

  const double dvpar = dxv[2], wvpar = w[2]; 
  const double *T_perp_over_m = &pkpm_prim[16]; 
  const double *T_perp_over_m_inv = &pkpm_prim[20]; 

  const double *p_perp_source = &pkpm_accel_vars[12]; 

  double alpha_G_1_source[12] = {0.0}; 

  alpha_G_1_source[0] = 1.4142135623730951*p_perp_source[0]-1.4142135623730951*div_b[0]*wvpar; 
  alpha_G_1_source[1] = 1.4142135623730951*p_perp_source[1]-1.4142135623730951*div_b[1]*wvpar; 
  alpha_G_1_source[2] = 1.4142135623730951*p_perp_source[2]-1.4142135623730951*div_b[2]*wvpar; 
  alpha_G_1_source[3] = -(0.408248290463863*div_b[0]*dvpar); 
  alpha_G_1_source[4] = 1.4142135623730951*p_perp_source[3]-1.4142135623730951*div_b[3]*wvpar; 
  alpha_G_1_source[5] = -(0.408248290463863*div_b[1]*dvpar); 
  alpha_G_1_source[6] = -(0.408248290463863*div_b[2]*dvpar); 
  alpha_G_1_source[8] = -(0.408248290463863*div_b[3]*dvpar); 

  const double *F_0 = &f[0]; 
  const double *G_1 = &f[12]; 
  const double *F_2 = &F_k_p_1[0]; 

  const double *nu_vtsq_sum = &nu_prim_moms_sum[4];

  double *out_G_1_mirror = &g_dist_source[0]; 
  double *out_G_1_vperp = &g_dist_source[12]; 
  double *out_F_1 = &F_k_m_1[0]; 

  double tmp_F_2_m_F_0[12] = {0.0}; 
  double tmp_T_perp_F_0[12] = {0.0}; 
  double tmp_T_perp_g_dist[12] = {0.0}; 
  double tmp_F_0_m_F_1[12] = {0.0}; 

  tmp_F_2_m_F_0[0] = F_2[0]-1.0*F_0[0]; 
  tmp_F_2_m_F_0[1] = F_2[1]-1.0*F_0[1]; 
  tmp_F_2_m_F_0[2] = F_2[2]-1.0*F_0[2]; 
  tmp_F_2_m_F_0[3] = F_2[3]-1.0*F_0[3]; 
  tmp_F_2_m_F_0[4] = F_2[4]-1.0*F_0[4]; 
  tmp_F_2_m_F_0[5] = F_2[5]-1.0*F_0[5]; 
  tmp_F_2_m_F_0[6] = F_2[6]-1.0*F_0[6]; 
  tmp_F_2_m_F_0[7] = F_2[7]-1.0*F_0[7]; 
  tmp_F_2_m_F_0[8] = F_2[8]-1.0*F_0[8]; 
  tmp_F_2_m_F_0[9] = F_2[9]-1.0*F_0[9]; 
  tmp_F_2_m_F_0[10] = F_2[10]-1.0*F_0[10]; 
  tmp_F_2_m_F_0[11] = F_2[11]-1.0*F_0[11]; 

  tmp_T_perp_F_0[0] = 0.5*T_perp_over_m[3]*tmp_F_2_m_F_0[4]+0.5*T_perp_over_m[2]*tmp_F_2_m_F_0[2]+0.5*T_perp_over_m[1]*tmp_F_2_m_F_0[1]+0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[0]; 
  tmp_T_perp_F_0[1] = 0.5*T_perp_over_m[2]*tmp_F_2_m_F_0[4]+0.5*tmp_F_2_m_F_0[2]*T_perp_over_m[3]+0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[1]+0.5*tmp_F_2_m_F_0[0]*T_perp_over_m[1]; 
  tmp_T_perp_F_0[2] = 0.5*T_perp_over_m[1]*tmp_F_2_m_F_0[4]+0.5*tmp_F_2_m_F_0[1]*T_perp_over_m[3]+0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[2]+0.5*tmp_F_2_m_F_0[0]*T_perp_over_m[2]; 
  tmp_T_perp_F_0[3] = 0.5*T_perp_over_m[3]*tmp_F_2_m_F_0[8]+0.5*T_perp_over_m[2]*tmp_F_2_m_F_0[6]+0.5*T_perp_over_m[1]*tmp_F_2_m_F_0[5]+0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[3]; 
  tmp_T_perp_F_0[4] = 0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[4]+0.5*tmp_F_2_m_F_0[0]*T_perp_over_m[3]+0.5*T_perp_over_m[1]*tmp_F_2_m_F_0[2]+0.5*tmp_F_2_m_F_0[1]*T_perp_over_m[2]; 
  tmp_T_perp_F_0[5] = 0.5*T_perp_over_m[2]*tmp_F_2_m_F_0[8]+0.5*T_perp_over_m[3]*tmp_F_2_m_F_0[6]+0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[5]+0.5*T_perp_over_m[1]*tmp_F_2_m_F_0[3]; 
  tmp_T_perp_F_0[6] = 0.5*T_perp_over_m[1]*tmp_F_2_m_F_0[8]+0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[6]+0.5*T_perp_over_m[3]*tmp_F_2_m_F_0[5]+0.5*T_perp_over_m[2]*tmp_F_2_m_F_0[3]; 
  tmp_T_perp_F_0[7] = 0.5*T_perp_over_m[3]*tmp_F_2_m_F_0[11]+0.5000000000000001*T_perp_over_m[2]*tmp_F_2_m_F_0[10]+0.5000000000000001*T_perp_over_m[1]*tmp_F_2_m_F_0[9]+0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[7]; 
  tmp_T_perp_F_0[8] = 0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[8]+0.5*T_perp_over_m[1]*tmp_F_2_m_F_0[6]+0.5*T_perp_over_m[2]*tmp_F_2_m_F_0[5]+0.5*T_perp_over_m[3]*tmp_F_2_m_F_0[3]; 
  tmp_T_perp_F_0[9] = 0.5000000000000001*T_perp_over_m[2]*tmp_F_2_m_F_0[11]+0.5*T_perp_over_m[3]*tmp_F_2_m_F_0[10]+0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[9]+0.5000000000000001*T_perp_over_m[1]*tmp_F_2_m_F_0[7]; 
  tmp_T_perp_F_0[10] = 0.5000000000000001*T_perp_over_m[1]*tmp_F_2_m_F_0[11]+0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[10]+0.5*T_perp_over_m[3]*tmp_F_2_m_F_0[9]+0.5000000000000001*T_perp_over_m[2]*tmp_F_2_m_F_0[7]; 
  tmp_T_perp_F_0[11] = 0.5*T_perp_over_m[0]*tmp_F_2_m_F_0[11]+0.5000000000000001*T_perp_over_m[1]*tmp_F_2_m_F_0[10]+0.5000000000000001*T_perp_over_m[2]*tmp_F_2_m_F_0[9]+0.5*T_perp_over_m[3]*tmp_F_2_m_F_0[7]; 

  tmp_T_perp_g_dist[0] = tmp_T_perp_F_0[0]+2.0*G_1[0]; 
  tmp_T_perp_g_dist[1] = tmp_T_perp_F_0[1]+2.0*G_1[1]; 
  tmp_T_perp_g_dist[2] = tmp_T_perp_F_0[2]+2.0*G_1[2]; 
  tmp_T_perp_g_dist[3] = tmp_T_perp_F_0[3]+2.0*G_1[3]; 
  tmp_T_perp_g_dist[4] = tmp_T_perp_F_0[4]+2.0*G_1[4]; 
  tmp_T_perp_g_dist[5] = tmp_T_perp_F_0[5]+2.0*G_1[5]; 
  tmp_T_perp_g_dist[6] = tmp_T_perp_F_0[6]+2.0*G_1[6]; 
  tmp_T_perp_g_dist[7] = tmp_T_perp_F_0[7]+2.0*G_1[7]; 
  tmp_T_perp_g_dist[8] = tmp_T_perp_F_0[8]+2.0*G_1[8]; 
  tmp_T_perp_g_dist[9] = tmp_T_perp_F_0[9]+2.0*G_1[9]; 
  tmp_T_perp_g_dist[10] = tmp_T_perp_F_0[10]+2.0*G_1[10]; 
  tmp_T_perp_g_dist[11] = tmp_T_perp_F_0[11]+2.0*G_1[11]; 

  out_G_1_mirror[0] = T_perp_over_m[3]*tmp_T_perp_g_dist[4]+T_perp_over_m[2]*tmp_T_perp_g_dist[2]+T_perp_over_m[1]*tmp_T_perp_g_dist[1]+T_perp_over_m[0]*tmp_T_perp_g_dist[0]; 
  out_G_1_mirror[1] = T_perp_over_m[2]*tmp_T_perp_g_dist[4]+tmp_T_perp_g_dist[2]*T_perp_over_m[3]+T_perp_over_m[0]*tmp_T_perp_g_dist[1]+tmp_T_perp_g_dist[0]*T_perp_over_m[1]; 
  out_G_1_mirror[2] = T_perp_over_m[1]*tmp_T_perp_g_dist[4]+tmp_T_perp_g_dist[1]*T_perp_over_m[3]+T_perp_over_m[0]*tmp_T_perp_g_dist[2]+tmp_T_perp_g_dist[0]*T_perp_over_m[2]; 
  out_G_1_mirror[3] = T_perp_over_m[3]*tmp_T_perp_g_dist[8]+T_perp_over_m[2]*tmp_T_perp_g_dist[6]+T_perp_over_m[1]*tmp_T_perp_g_dist[5]+T_perp_over_m[0]*tmp_T_perp_g_dist[3]; 
  out_G_1_mirror[4] = T_perp_over_m[0]*tmp_T_perp_g_dist[4]+tmp_T_perp_g_dist[0]*T_perp_over_m[3]+T_perp_over_m[1]*tmp_T_perp_g_dist[2]+tmp_T_perp_g_dist[1]*T_perp_over_m[2]; 
  out_G_1_mirror[5] = T_perp_over_m[2]*tmp_T_perp_g_dist[8]+T_perp_over_m[3]*tmp_T_perp_g_dist[6]+T_perp_over_m[0]*tmp_T_perp_g_dist[5]+T_perp_over_m[1]*tmp_T_perp_g_dist[3]; 
  out_G_1_mirror[6] = T_perp_over_m[1]*tmp_T_perp_g_dist[8]+T_perp_over_m[0]*tmp_T_perp_g_dist[6]+T_perp_over_m[3]*tmp_T_perp_g_dist[5]+T_perp_over_m[2]*tmp_T_perp_g_dist[3]; 
  out_G_1_mirror[7] = T_perp_over_m[3]*tmp_T_perp_g_dist[11]+1.0000000000000002*T_perp_over_m[2]*tmp_T_perp_g_dist[10]+1.0000000000000002*T_perp_over_m[1]*tmp_T_perp_g_dist[9]+T_perp_over_m[0]*tmp_T_perp_g_dist[7]; 
  out_G_1_mirror[8] = T_perp_over_m[0]*tmp_T_perp_g_dist[8]+T_perp_over_m[1]*tmp_T_perp_g_dist[6]+T_perp_over_m[2]*tmp_T_perp_g_dist[5]+T_perp_over_m[3]*tmp_T_perp_g_dist[3]; 
  out_G_1_mirror[9] = 1.0000000000000002*T_perp_over_m[2]*tmp_T_perp_g_dist[11]+T_perp_over_m[3]*tmp_T_perp_g_dist[10]+T_perp_over_m[0]*tmp_T_perp_g_dist[9]+1.0000000000000002*T_perp_over_m[1]*tmp_T_perp_g_dist[7]; 
  out_G_1_mirror[10] = 1.0000000000000002*T_perp_over_m[1]*tmp_T_perp_g_dist[11]+T_perp_over_m[0]*tmp_T_perp_g_dist[10]+T_perp_over_m[3]*tmp_T_perp_g_dist[9]+1.0000000000000002*T_perp_over_m[2]*tmp_T_perp_g_dist[7]; 
  out_G_1_mirror[11] = T_perp_over_m[0]*tmp_T_perp_g_dist[11]+1.0000000000000002*T_perp_over_m[1]*tmp_T_perp_g_dist[10]+1.0000000000000002*T_perp_over_m[2]*tmp_T_perp_g_dist[9]+T_perp_over_m[3]*tmp_T_perp_g_dist[7]; 

  out_G_1_vperp[0] = 0.3535533905932737*G_1[8]*alpha_G_1_source[8]+0.3535533905932737*G_1[6]*alpha_G_1_source[6]+0.3535533905932737*G_1[5]*alpha_G_1_source[5]+0.3535533905932737*G_1[4]*alpha_G_1_source[4]+nu_vtsq_sum[3]*F_0[4]+0.3535533905932737*G_1[3]*alpha_G_1_source[3]+F_0[2]*nu_vtsq_sum[2]+0.3535533905932737*G_1[2]*alpha_G_1_source[2]+F_0[1]*nu_vtsq_sum[1]+0.3535533905932737*G_1[1]*alpha_G_1_source[1]+F_0[0]*nu_vtsq_sum[0]+0.3535533905932737*G_1[0]*alpha_G_1_source[0]; 
  out_G_1_vperp[1] = 0.3535533905932737*G_1[6]*alpha_G_1_source[8]+0.3535533905932737*alpha_G_1_source[6]*G_1[8]+0.3535533905932737*G_1[3]*alpha_G_1_source[5]+0.3535533905932737*alpha_G_1_source[3]*G_1[5]+0.3535533905932737*G_1[2]*alpha_G_1_source[4]+0.3535533905932737*alpha_G_1_source[2]*G_1[4]+nu_vtsq_sum[2]*F_0[4]+F_0[2]*nu_vtsq_sum[3]+F_0[0]*nu_vtsq_sum[1]+0.3535533905932737*G_1[0]*alpha_G_1_source[1]+0.3535533905932737*alpha_G_1_source[0]*G_1[1]+nu_vtsq_sum[0]*F_0[1]; 
  out_G_1_vperp[2] = 0.3535533905932737*G_1[5]*alpha_G_1_source[8]+0.3535533905932737*alpha_G_1_source[5]*G_1[8]+0.3535533905932737*G_1[3]*alpha_G_1_source[6]+0.3535533905932737*alpha_G_1_source[3]*G_1[6]+0.3535533905932737*G_1[1]*alpha_G_1_source[4]+0.3535533905932737*alpha_G_1_source[1]*G_1[4]+nu_vtsq_sum[1]*F_0[4]+F_0[1]*nu_vtsq_sum[3]+F_0[0]*nu_vtsq_sum[2]+0.3535533905932737*G_1[0]*alpha_G_1_source[2]+0.3535533905932737*alpha_G_1_source[0]*G_1[2]+nu_vtsq_sum[0]*F_0[2]; 
  out_G_1_vperp[3] = 0.3162277660168379*alpha_G_1_source[8]*G_1[11]+0.31622776601683794*alpha_G_1_source[6]*G_1[10]+0.31622776601683794*alpha_G_1_source[5]*G_1[9]+0.3535533905932737*G_1[4]*alpha_G_1_source[8]+0.3535533905932737*alpha_G_1_source[4]*G_1[8]+nu_vtsq_sum[3]*F_0[8]+0.3162277660168379*alpha_G_1_source[3]*G_1[7]+0.3535533905932737*G_1[2]*alpha_G_1_source[6]+0.3535533905932737*alpha_G_1_source[2]*G_1[6]+nu_vtsq_sum[2]*F_0[6]+0.3535533905932737*G_1[1]*alpha_G_1_source[5]+0.3535533905932737*alpha_G_1_source[1]*G_1[5]+nu_vtsq_sum[1]*F_0[5]+0.3535533905932737*G_1[0]*alpha_G_1_source[3]+0.3535533905932737*alpha_G_1_source[0]*G_1[3]+nu_vtsq_sum[0]*F_0[3]; 
  out_G_1_vperp[4] = 0.3535533905932737*G_1[3]*alpha_G_1_source[8]+0.3535533905932737*alpha_G_1_source[3]*G_1[8]+0.3535533905932737*G_1[5]*alpha_G_1_source[6]+0.3535533905932737*alpha_G_1_source[5]*G_1[6]+0.3535533905932737*G_1[0]*alpha_G_1_source[4]+0.3535533905932737*alpha_G_1_source[0]*G_1[4]+nu_vtsq_sum[0]*F_0[4]+F_0[0]*nu_vtsq_sum[3]+F_0[1]*nu_vtsq_sum[2]+0.3535533905932737*G_1[1]*alpha_G_1_source[2]+0.3535533905932737*alpha_G_1_source[1]*G_1[2]+nu_vtsq_sum[1]*F_0[2]; 
  out_G_1_vperp[5] = 0.3162277660168379*alpha_G_1_source[6]*G_1[11]+0.31622776601683794*alpha_G_1_source[8]*G_1[10]+0.31622776601683794*alpha_G_1_source[3]*G_1[9]+0.3535533905932737*G_1[2]*alpha_G_1_source[8]+0.3535533905932737*alpha_G_1_source[2]*G_1[8]+nu_vtsq_sum[2]*F_0[8]+0.3162277660168379*alpha_G_1_source[5]*G_1[7]+0.3535533905932737*G_1[4]*alpha_G_1_source[6]+0.3535533905932737*alpha_G_1_source[4]*G_1[6]+nu_vtsq_sum[3]*F_0[6]+0.3535533905932737*G_1[0]*alpha_G_1_source[5]+0.3535533905932737*alpha_G_1_source[0]*G_1[5]+nu_vtsq_sum[0]*F_0[5]+0.3535533905932737*G_1[1]*alpha_G_1_source[3]+0.3535533905932737*alpha_G_1_source[1]*G_1[3]+nu_vtsq_sum[1]*F_0[3]; 
  out_G_1_vperp[6] = 0.3162277660168379*alpha_G_1_source[5]*G_1[11]+0.31622776601683794*alpha_G_1_source[3]*G_1[10]+0.31622776601683794*alpha_G_1_source[8]*G_1[9]+0.3535533905932737*G_1[1]*alpha_G_1_source[8]+0.3535533905932737*alpha_G_1_source[1]*G_1[8]+nu_vtsq_sum[1]*F_0[8]+0.3162277660168379*alpha_G_1_source[6]*G_1[7]+0.3535533905932737*G_1[0]*alpha_G_1_source[6]+0.3535533905932737*alpha_G_1_source[0]*G_1[6]+nu_vtsq_sum[0]*F_0[6]+0.3535533905932737*G_1[4]*alpha_G_1_source[5]+0.3535533905932737*alpha_G_1_source[4]*G_1[5]+nu_vtsq_sum[3]*F_0[5]+0.3535533905932737*G_1[2]*alpha_G_1_source[3]+0.3535533905932737*alpha_G_1_source[2]*G_1[3]+nu_vtsq_sum[2]*F_0[3]; 
  out_G_1_vperp[7] = 0.3535533905932737*alpha_G_1_source[4]*G_1[11]+nu_vtsq_sum[3]*F_0[11]+0.3535533905932737*alpha_G_1_source[2]*G_1[10]+1.0000000000000002*nu_vtsq_sum[2]*F_0[10]+0.3535533905932737*alpha_G_1_source[1]*G_1[9]+1.0000000000000002*nu_vtsq_sum[1]*F_0[9]+0.3162277660168379*G_1[8]*alpha_G_1_source[8]+0.3535533905932737*alpha_G_1_source[0]*G_1[7]+nu_vtsq_sum[0]*F_0[7]+0.3162277660168379*G_1[6]*alpha_G_1_source[6]+0.3162277660168379*G_1[5]*alpha_G_1_source[5]+0.3162277660168379*G_1[3]*alpha_G_1_source[3]; 
  out_G_1_vperp[8] = 0.3162277660168379*alpha_G_1_source[3]*G_1[11]+0.31622776601683794*alpha_G_1_source[5]*G_1[10]+0.31622776601683794*alpha_G_1_source[6]*G_1[9]+0.3162277660168379*G_1[7]*alpha_G_1_source[8]+0.3535533905932737*G_1[0]*alpha_G_1_source[8]+0.3535533905932737*alpha_G_1_source[0]*G_1[8]+nu_vtsq_sum[0]*F_0[8]+0.3535533905932737*G_1[1]*alpha_G_1_source[6]+0.3535533905932737*alpha_G_1_source[1]*G_1[6]+nu_vtsq_sum[1]*F_0[6]+0.3535533905932737*G_1[2]*alpha_G_1_source[5]+0.3535533905932737*alpha_G_1_source[2]*G_1[5]+nu_vtsq_sum[2]*F_0[5]+0.3535533905932737*G_1[3]*alpha_G_1_source[4]+0.3535533905932737*alpha_G_1_source[3]*G_1[4]+F_0[3]*nu_vtsq_sum[3]; 
  out_G_1_vperp[9] = 0.3535533905932737*alpha_G_1_source[2]*G_1[11]+1.0000000000000002*nu_vtsq_sum[2]*F_0[11]+0.3535533905932737*alpha_G_1_source[4]*G_1[10]+nu_vtsq_sum[3]*F_0[10]+0.3535533905932737*alpha_G_1_source[0]*G_1[9]+nu_vtsq_sum[0]*F_0[9]+0.31622776601683794*G_1[6]*alpha_G_1_source[8]+0.31622776601683794*alpha_G_1_source[6]*G_1[8]+0.3535533905932737*alpha_G_1_source[1]*G_1[7]+1.0000000000000002*nu_vtsq_sum[1]*F_0[7]+0.31622776601683794*G_1[3]*alpha_G_1_source[5]+0.31622776601683794*alpha_G_1_source[3]*G_1[5]; 
  out_G_1_vperp[10] = 0.3535533905932737*alpha_G_1_source[1]*G_1[11]+1.0000000000000002*nu_vtsq_sum[1]*F_0[11]+0.3535533905932737*alpha_G_1_source[0]*G_1[10]+nu_vtsq_sum[0]*F_0[10]+0.3535533905932737*alpha_G_1_source[4]*G_1[9]+nu_vtsq_sum[3]*F_0[9]+0.31622776601683794*G_1[5]*alpha_G_1_source[8]+0.31622776601683794*alpha_G_1_source[5]*G_1[8]+0.3535533905932737*alpha_G_1_source[2]*G_1[7]+1.0000000000000002*nu_vtsq_sum[2]*F_0[7]+0.31622776601683794*G_1[3]*alpha_G_1_source[6]+0.31622776601683794*alpha_G_1_source[3]*G_1[6]; 
  out_G_1_vperp[11] = 0.3535533905932737*alpha_G_1_source[0]*G_1[11]+nu_vtsq_sum[0]*F_0[11]+0.3535533905932737*alpha_G_1_source[1]*G_1[10]+1.0000000000000002*nu_vtsq_sum[1]*F_0[10]+0.3535533905932737*alpha_G_1_source[2]*G_1[9]+1.0000000000000002*nu_vtsq_sum[2]*F_0[9]+0.3162277660168379*G_1[3]*alpha_G_1_source[8]+0.3162277660168379*alpha_G_1_source[3]*G_1[8]+0.3535533905932737*alpha_G_1_source[4]*G_1[7]+nu_vtsq_sum[3]*F_0[7]+0.3162277660168379*G_1[5]*alpha_G_1_source[6]+0.3162277660168379*alpha_G_1_source[5]*G_1[6]; 

  tmp_F_0_m_F_1[0] = 0.5*T_perp_over_m_inv[3]*G_1[4]+0.5*G_1[2]*T_perp_over_m_inv[2]+0.5*G_1[1]*T_perp_over_m_inv[1]+0.5*G_1[0]*T_perp_over_m_inv[0]; 
  tmp_F_0_m_F_1[1] = 0.5*T_perp_over_m_inv[2]*G_1[4]+0.5*G_1[2]*T_perp_over_m_inv[3]+0.5*G_1[0]*T_perp_over_m_inv[1]+0.5*T_perp_over_m_inv[0]*G_1[1]; 
  tmp_F_0_m_F_1[2] = 0.5*T_perp_over_m_inv[1]*G_1[4]+0.5*G_1[1]*T_perp_over_m_inv[3]+0.5*G_1[0]*T_perp_over_m_inv[2]+0.5*T_perp_over_m_inv[0]*G_1[2]; 
  tmp_F_0_m_F_1[3] = 0.5*T_perp_over_m_inv[3]*G_1[8]+0.5*T_perp_over_m_inv[2]*G_1[6]+0.5*T_perp_over_m_inv[1]*G_1[5]+0.5*T_perp_over_m_inv[0]*G_1[3]; 
  tmp_F_0_m_F_1[4] = 0.5*T_perp_over_m_inv[0]*G_1[4]+0.5*G_1[0]*T_perp_over_m_inv[3]+0.5*G_1[1]*T_perp_over_m_inv[2]+0.5*T_perp_over_m_inv[1]*G_1[2]; 
  tmp_F_0_m_F_1[5] = 0.5*T_perp_over_m_inv[2]*G_1[8]+0.5*T_perp_over_m_inv[3]*G_1[6]+0.5*T_perp_over_m_inv[0]*G_1[5]+0.5*T_perp_over_m_inv[1]*G_1[3]; 
  tmp_F_0_m_F_1[6] = 0.5*T_perp_over_m_inv[1]*G_1[8]+0.5*T_perp_over_m_inv[0]*G_1[6]+0.5*T_perp_over_m_inv[3]*G_1[5]+0.5*T_perp_over_m_inv[2]*G_1[3]; 
  tmp_F_0_m_F_1[7] = 0.5*T_perp_over_m_inv[3]*G_1[11]+0.5000000000000001*T_perp_over_m_inv[2]*G_1[10]+0.5000000000000001*T_perp_over_m_inv[1]*G_1[9]+0.5*T_perp_over_m_inv[0]*G_1[7]; 
  tmp_F_0_m_F_1[8] = 0.5*T_perp_over_m_inv[0]*G_1[8]+0.5*T_perp_over_m_inv[1]*G_1[6]+0.5*T_perp_over_m_inv[2]*G_1[5]+0.5*G_1[3]*T_perp_over_m_inv[3]; 
  tmp_F_0_m_F_1[9] = 0.5000000000000001*T_perp_over_m_inv[2]*G_1[11]+0.5*T_perp_over_m_inv[3]*G_1[10]+0.5*T_perp_over_m_inv[0]*G_1[9]+0.5000000000000001*T_perp_over_m_inv[1]*G_1[7]; 
  tmp_F_0_m_F_1[10] = 0.5000000000000001*T_perp_over_m_inv[1]*G_1[11]+0.5*T_perp_over_m_inv[0]*G_1[10]+0.5*T_perp_over_m_inv[3]*G_1[9]+0.5000000000000001*T_perp_over_m_inv[2]*G_1[7]; 
  tmp_F_0_m_F_1[11] = 0.5*T_perp_over_m_inv[0]*G_1[11]+0.5000000000000001*T_perp_over_m_inv[1]*G_1[10]+0.5000000000000001*T_perp_over_m_inv[2]*G_1[9]+0.5*T_perp_over_m_inv[3]*G_1[7]; 

  out_F_1[0] = F_0[0]-1.0*tmp_F_0_m_F_1[0]; 
  out_F_1[1] = F_0[1]-1.0*tmp_F_0_m_F_1[1]; 
  out_F_1[2] = F_0[2]-1.0*tmp_F_0_m_F_1[2]; 
  out_F_1[3] = F_0[3]-1.0*tmp_F_0_m_F_1[3]; 
  out_F_1[4] = F_0[4]-1.0*tmp_F_0_m_F_1[4]; 
  out_F_1[5] = F_0[5]-1.0*tmp_F_0_m_F_1[5]; 
  out_F_1[6] = F_0[6]-1.0*tmp_F_0_m_F_1[6]; 
  out_F_1[7] = F_0[7]-1.0*tmp_F_0_m_F_1[7]; 
  out_F_1[8] = F_0[8]-1.0*tmp_F_0_m_F_1[8]; 
  out_F_1[9] = F_0[9]-1.0*tmp_F_0_m_F_1[9]; 
  out_F_1[10] = F_0[10]-1.0*tmp_F_0_m_F_1[10]; 
  out_F_1[11] = F_0[11]-1.0*tmp_F_0_m_F_1[11]; 
} 
