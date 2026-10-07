#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_tensor_p1.h> 
GKYL_CU_DH int nc_hamil_phase_alpha_quad_vy_1x3v_tensor_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 14; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double *poisson_tensor_conf_x1 = &poisson_tensor_conf[2]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv12 = 2.0/dxv[3]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[18]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[22]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[24]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[28]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[30]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[34]; 
  for (int i = tid; i < 2; i += nthreads) { 
    double px0_q = 0.0; 
    for (int a = 0; a < 2; ++a) { 
      px0_q += vst_1x3v_tensor_p1_conf_ev[i*2 + a]*poisson_tensor_conf_x1[a]; 
    } 
    for (int a = 0; a < 2; ++a) O[a*2 + i] = -(px0_q*vst_1x3v_tensor_p1_ph_v1_CmDx0[i*2 + a]*dx10*jacob_cx_inv); 
    double p0_q = 0.0; 
    double p2_q = 0.0; 
    double p3_q = 0.0; 
    double p5_q = 0.0; 
    double p6_q = 0.0; 
    double p8_q = 0.0; 
    for (int a = 0; a < 2; ++a) { 
      p0_q += vst_1x3v_tensor_p1_conf_ev[i*2 + a]*poisson_tensor_conf_0[a]; 
      p2_q += vst_1x3v_tensor_p1_conf_ev[i*2 + a]*poisson_tensor_conf_2[a]; 
      p3_q += vst_1x3v_tensor_p1_conf_ev[i*2 + a]*poisson_tensor_conf_3[a]; 
      p5_q += vst_1x3v_tensor_p1_conf_ev[i*2 + a]*poisson_tensor_conf_5[a]; 
      p6_q += vst_1x3v_tensor_p1_conf_ev[i*2 + a]*poisson_tensor_conf_6[a]; 
      p8_q += vst_1x3v_tensor_p1_conf_ev[i*2 + a]*poisson_tensor_conf_8[a]; 
    } 
    for (int a = 0; a < 2; ++a) O[(2 + a)*2 + i] = -p0_q*vst_1x3v_tensor_p1_ph_v1_Cm[i*2 + a]; 
    for (int a = 0; a < 2; ++a) O[(4 + a)*2 + i] = -p3_q*vst_1x3v_tensor_p1_ph_v1_Cm[i*2 + a]; 
    for (int a = 0; a < 2; ++a) O[(6 + a)*2 + i] = -p6_q*vst_1x3v_tensor_p1_ph_v1_Cm[i*2 + a]; 
    for (int a = 0; a < 2; ++a) O[(8 + a)*2 + i] = p2_q*vst_1x3v_tensor_p1_ph_v1_Cm[i*2 + a]; 
    for (int a = 0; a < 2; ++a) O[(10 + a)*2 + i] = p5_q*vst_1x3v_tensor_p1_ph_v1_Cm[i*2 + a]; 
    for (int a = 0; a < 2; ++a) O[(12 + a)*2 + i] = p8_q*vst_1x3v_tensor_p1_ph_v1_Cm[i*2 + a]; 
  } 
  for (int j = tid; j < 9; j += nthreads) { 
    double G[2]; 
    double Gd0[2]; 
    double Gd2[2]; 
    for (int a = 0; a < 2; ++a) { G[a] = 0.0; Gd0[a] = 0.0; Gd2[a] = 0.0; } 
    for (int k = 0; k < 54; ++k) { 
      const int a = vst_1x3v_tensor_p1_ph_v1_cmap[k]; 
      G[a] += vst_1x3v_tensor_p1_ph_v1_V[j*27 + vst_1x3v_tensor_p1_ph_v1_vrmap[k]]*(vst_1x3v_tensor_p1_ph_v1_coefr[k]*hamil[k]); 
      Gd0[a] += vst_1x3v_tensor_p1_ph_v1_Vd0[j*19 + vst_1x3v_tensor_p1_ph_v1_vrd0map[k]]*(vst_1x3v_tensor_p1_ph_v1_dcoefr0[k]*hamil[k]); 
      Gd2[a] += vst_1x3v_tensor_p1_ph_v1_Vd2[j*19 + vst_1x3v_tensor_p1_ph_v1_vrd2map[k]]*(vst_1x3v_tensor_p1_ph_v1_dcoefr2[k]*hamil[k]); 
    } 
    for (int a = 0; a < 2; ++a) I[a*9 + j] = G[a]; 
    const double jacob_vx_inv = 1.0/jacob_vel_surf[0 + j/3]; 
    const double jacob_vz_inv = 1.0/jacob_vel_surf[6 + j%3]; 
    const double xn0 = vst_1x3v_tensor_p1_vel_nodes_v1[j*2 + 0]; 
    const double xn0_sq = xn0*xn0; 
    const double vt1 = 4.677071733467426*vmap_v0[3]*xn0*xn0_sq+2.371708245126284*vmap_v0[2]*xn0_sq-2.806243040080455*vmap_v0[3]*xn0+1.224744871391589*vmap_v0[1]*xn0-0.7905694150420947*vmap_v0[2]+0.7071067811865475*vmap_v0[0]; 
    const double vt2 = -(1.8708286933869707*vmap_v1[3])+1.5811388300841895*vmap_v1[2]-1.224744871391589*vmap_v1[1]+0.7071067811865475*vmap_v1[0]; 
    const double xn2 = vst_1x3v_tensor_p1_vel_nodes_v1[j*2 + 1]; 
    const double xn2_sq = xn2*xn2; 
    const double vt3 = 4.677071733467426*vmap_v2[3]*xn2*xn2_sq+2.371708245126284*vmap_v2[2]*xn2_sq-2.806243040080455*vmap_v2[3]*xn2+1.224744871391589*vmap_v2[1]*xn2-0.7905694150420947*vmap_v2[2]+0.7071067811865475*vmap_v2[0]; 
    for (int a = 0; a < 2; ++a) I[(2 + a)*9 + j] = vt1*Gd0[a]*dv10*jacob_vx_inv; 
    for (int a = 0; a < 2; ++a) I[(4 + a)*9 + j] = vt2*Gd0[a]*dv10*jacob_vx_inv; 
    for (int a = 0; a < 2; ++a) I[(6 + a)*9 + j] = vt3*Gd0[a]*dv10*jacob_vx_inv; 
    for (int a = 0; a < 2; ++a) I[(8 + a)*9 + j] = vt1*Gd2[a]*dv12*jacob_vz_inv; 
    for (int a = 0; a < 2; ++a) I[(10 + a)*9 + j] = vt2*Gd2[a]*dv12*jacob_vz_inv; 
    for (int a = 0; a < 2; ++a) I[(12 + a)*9 + j] = vt3*Gd2[a]*dv12*jacob_vz_inv; 
  } 
  return 14; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[28]; 
  double I[126]; 
  nc_hamil_phase_alpha_quad_vy_1x3v_tensor_p1_shared(0, 1, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 14; ++t) alpha += O[t*2 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
