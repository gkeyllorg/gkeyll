#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH int nc_hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2_shared(int tid, int nthreads, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 21; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double *poisson_tensor_conf_x0 = &poisson_tensor_conf[0]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[4]; 
  const double jacob_vz_inv = 1.0/jacob_vel_surf[8]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[27]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[30]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[36]; 
  const double *poisson_tensor_conf_4 = &poisson_tensor_conf[39]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[45]; 
  const double *poisson_tensor_conf_7 = &poisson_tensor_conf[48]; 
  for (int i = tid; i < 4; i += nthreads) { 
    double px0_q = 0.0; 
    for (int a = 0; a < 3; ++a) { 
      px0_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_x0[a]; 
    } 
    for (int a = 0; a < 3; ++a) O[a*4 + i] = -(px0_q*vst_1x3v_ser_p2_ho_ph_v0_CmDx0[i*3 + a]*dx10*jacob_cx_inv); 
    double p0_q = 0.0; 
    double p1_q = 0.0; 
    double p3_q = 0.0; 
    double p4_q = 0.0; 
    double p6_q = 0.0; 
    double p7_q = 0.0; 
    for (int a = 0; a < 3; ++a) { 
      p0_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_0[a]; 
      p1_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_1[a]; 
      p3_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_3[a]; 
      p4_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_4[a]; 
      p6_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_6[a]; 
      p7_q += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*poisson_tensor_conf_7[a]; 
    } 
    for (int a = 0; a < 3; ++a) O[(3 + a)*4 + i] = p0_q*vst_1x3v_ser_p2_ho_ph_v0_Cm[i*3 + a]; 
    for (int a = 0; a < 3; ++a) O[(6 + a)*4 + i] = p3_q*vst_1x3v_ser_p2_ho_ph_v0_Cm[i*3 + a]; 
    for (int a = 0; a < 3; ++a) O[(9 + a)*4 + i] = p6_q*vst_1x3v_ser_p2_ho_ph_v0_Cm[i*3 + a]; 
    for (int a = 0; a < 3; ++a) O[(12 + a)*4 + i] = p1_q*vst_1x3v_ser_p2_ho_ph_v0_Cm[i*3 + a]; 
    for (int a = 0; a < 3; ++a) O[(15 + a)*4 + i] = p4_q*vst_1x3v_ser_p2_ho_ph_v0_Cm[i*3 + a]; 
    for (int a = 0; a < 3; ++a) O[(18 + a)*4 + i] = p7_q*vst_1x3v_ser_p2_ho_ph_v0_Cm[i*3 + a]; 
  } 
  for (int j = tid; j < 16; j += nthreads) { 
    double G[3]; 
    double Gd1[3]; 
    double Gd2[3]; 
    for (int a = 0; a < 3; ++a) { G[a] = 0.0; Gd1[a] = 0.0; Gd2[a] = 0.0; } 
    for (int k = 0; k < 48; ++k) { 
      const int a = vst_1x3v_ser_p2_ho_ph_v0_cmap[k]; 
      G[a] += vst_1x3v_ser_p2_ho_ph_v0_V[j*20 + vst_1x3v_ser_p2_ho_ph_v0_vrmap[k]]*(vst_1x3v_ser_p2_ho_ph_v0_coefr[k]*hamil[k]); 
      Gd1[a] += vst_1x3v_ser_p2_ho_ph_v0_Vd1[j*13 + vst_1x3v_ser_p2_ho_ph_v0_vrd1map[k]]*(vst_1x3v_ser_p2_ho_ph_v0_dcoefr1[k]*hamil[k]); 
      Gd2[a] += vst_1x3v_ser_p2_ho_ph_v0_Vd2[j*13 + vst_1x3v_ser_p2_ho_ph_v0_vrd2map[k]]*(vst_1x3v_ser_p2_ho_ph_v0_dcoefr2[k]*hamil[k]); 
    } 
    for (int a = 0; a < 3; ++a) I[a*16 + j] = G[a]; 
    const double vt1 = 0.7071067811865475*vmap_v0[0] - 1.224744871391589*vmap_v0[1]; 
    const double vt2 = 0.7071067811865475*vmap_v1[0] + 1.224744871391589*vmap_v1[1]*vst_1x3v_ser_p2_ho_vel_nodes_v0[j*2 + 0]; 
    const double vt3 = 0.7071067811865475*vmap_v2[0] + 1.224744871391589*vmap_v2[1]*vst_1x3v_ser_p2_ho_vel_nodes_v0[j*2 + 1]; 
    for (int a = 0; a < 3; ++a) I[(3 + a)*16 + j] = vt1*Gd1[a]*dv11*jacob_vy_inv; 
    for (int a = 0; a < 3; ++a) I[(6 + a)*16 + j] = vt2*Gd1[a]*dv11*jacob_vy_inv; 
    for (int a = 0; a < 3; ++a) I[(9 + a)*16 + j] = vt3*Gd1[a]*dv11*jacob_vy_inv; 
    for (int a = 0; a < 3; ++a) I[(12 + a)*16 + j] = vt1*Gd2[a]*dv12*jacob_vz_inv; 
    for (int a = 0; a < 3; ++a) I[(15 + a)*16 + j] = vt2*Gd2[a]*dv12*jacob_vz_inv; 
    for (int a = 0; a < 3; ++a) I[(18 + a)*16 + j] = vt3*Gd2[a]*dv12*jacob_vz_inv; 
  } 
  return 21; 
} 

GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[84]; 
  double I[336]; 
  nc_hamil_phase_ho_alpha_quad_vx_1x3v_ser_p2_shared(0, 1, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 21; ++t) alpha += O[t*4 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
