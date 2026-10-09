#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH int nc_hamil_phase_alpha_quad_vy_1x3v_ser_p2_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 21; 
  O += off*3; 
  I += off*9; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double *poisson_tensor_conf_x1 = &poisson_tensor_conf[3]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv12 = 2.0/dxv[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double jacob_vz_inv = 1.0/jacob_vel_surf[8]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[27]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[33]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[36]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[42]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[45]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[51]; 
  if (tid < 3) { 
    const int i = tid; 
    double px0_q = 0.0; 
    for (int a = 0; a < 3; ++a) { 
      px0_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_x1[a]; 
    } 
    for (int a = 0; a < 3; ++a) O[a*3 + i] = -(px0_q*vst_1x3v_ser_p2_ph_v1_CmDx0[i*3 + a]*dx10*jacob_cx_inv); 
    double p0_q = 0.0; 
    double p2_q = 0.0; 
    double p3_q = 0.0; 
    double p5_q = 0.0; 
    double p6_q = 0.0; 
    double p8_q = 0.0; 
    for (int a = 0; a < 3; ++a) { 
      p0_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_0[a]; 
      p2_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_2[a]; 
      p3_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_3[a]; 
      p5_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_5[a]; 
      p6_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_6[a]; 
      p8_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_8[a]; 
    } 
    for (int a = 0; a < 3; ++a) O[(3 + a)*3 + i] = -p0_q*vst_1x3v_ser_p2_ph_v1_Cm[i*3 + a]; 
    for (int a = 0; a < 3; ++a) O[(6 + a)*3 + i] = -p3_q*vst_1x3v_ser_p2_ph_v1_Cm[i*3 + a]; 
    for (int a = 0; a < 3; ++a) O[(9 + a)*3 + i] = -p6_q*vst_1x3v_ser_p2_ph_v1_Cm[i*3 + a]; 
    for (int a = 0; a < 3; ++a) O[(12 + a)*3 + i] = p2_q*vst_1x3v_ser_p2_ph_v1_Cm[i*3 + a]; 
    for (int a = 0; a < 3; ++a) O[(15 + a)*3 + i] = p5_q*vst_1x3v_ser_p2_ph_v1_Cm[i*3 + a]; 
    for (int a = 0; a < 3; ++a) O[(18 + a)*3 + i] = p8_q*vst_1x3v_ser_p2_ph_v1_Cm[i*3 + a]; 
  } 
  if (tid < 9) { 
    const int j = tid; 
    double G[3]; 
    double Gd0[3]; 
    double Gd2[3]; 
    for (int a = 0; a < 3; ++a) { G[a] = 0.0; Gd0[a] = 0.0; Gd2[a] = 0.0; } 
    for (int k = 0; k < 48; ++k) { 
      const int a = vst_1x3v_ser_p2_ph_v1_cmap[k]; 
      G[a] += vst_1x3v_ser_p2_ph_v1_V[j*20 + vst_1x3v_ser_p2_ph_v1_vrmap[k]]*(vst_1x3v_ser_p2_ph_v1_coefr[k]*hamil[k]); 
      Gd0[a] += vst_1x3v_ser_p2_ph_v1_Vd0[j*13 + vst_1x3v_ser_p2_ph_v1_vrd0map[k]]*(vst_1x3v_ser_p2_ph_v1_dcoefr0[k]*hamil[k]); 
      Gd2[a] += vst_1x3v_ser_p2_ph_v1_Vd2[j*13 + vst_1x3v_ser_p2_ph_v1_vrd2map[k]]*(vst_1x3v_ser_p2_ph_v1_dcoefr2[k]*hamil[k]); 
    } 
    for (int a = 0; a < 3; ++a) I[a*9 + j] = G[a]; 
    const double vt1 = 0.7071067811865475*vmap_v0[0] + 1.224744871391589*vmap_v0[1]*vst_1x3v_ser_p2_vel_nodes_v1[j*2 + 0]; 
    const double vt2 = 0.7071067811865475*vmap_v1[0] - 1.224744871391589*vmap_v1[1]; 
    const double vt3 = 0.7071067811865475*vmap_v2[0] + 1.224744871391589*vmap_v2[1]*vst_1x3v_ser_p2_vel_nodes_v1[j*2 + 1]; 
    for (int a = 0; a < 3; ++a) I[(3 + a)*9 + j] = vt1*Gd0[a]*dv10*jacob_vx_inv; 
    for (int a = 0; a < 3; ++a) I[(6 + a)*9 + j] = vt2*Gd0[a]*dv10*jacob_vx_inv; 
    for (int a = 0; a < 3; ++a) I[(9 + a)*9 + j] = vt3*Gd0[a]*dv10*jacob_vx_inv; 
    for (int a = 0; a < 3; ++a) I[(12 + a)*9 + j] = vt1*Gd2[a]*dv12*jacob_vz_inv; 
    for (int a = 0; a < 3; ++a) I[(15 + a)*9 + j] = vt2*Gd2[a]*dv12*jacob_vz_inv; 
    for (int a = 0; a < 3; ++a) I[(18 + a)*9 + j] = vt3*Gd2[a]*dv12*jacob_vz_inv; 
  } 
  return 21; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[63]; 
  double I[189]; 
  for (int tid = 0; tid < 9; ++tid) nc_hamil_phase_alpha_quad_vy_1x3v_ser_p2_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 21; ++t) alpha += O[t*3 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
