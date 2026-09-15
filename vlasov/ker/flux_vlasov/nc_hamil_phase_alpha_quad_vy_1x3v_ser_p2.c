#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double nc_hamil_phase_alpha_quad_vy_1x3v_ser_p2_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
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
  double dH_dv0 = 0.0; 
  for (int a = 0; a < 3; ++a) dH_dv0 += vst_1x3v_ser_p2_ph_v1_Cm[i*3 + a]*Gd0[a]; 
  double dH_dv2 = 0.0; 
  for (int a = 0; a < 3; ++a) dH_dv2 += vst_1x3v_ser_p2_ph_v1_Cm[i*3 + a]*Gd2[a]; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double *poisson_tensor_conf_x1 = &poisson_tensor_conf[3]; 
  double px0_q = 0.0; 
  for (int a = 0; a < 3; ++a) { 
    px0_q += vst_1x3v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_x1[a]; 
  } 
  double dH_dx0 = 0.0; 
  for (int a = 0; a < 3; ++a) dH_dx0 += vst_1x3v_ser_p2_ph_v1_CmDx0[i*3 + a]*G[a]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv12 = 2.0/dxv[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double jacob_vz_inv = 1.0/jacob_vel_surf[8]; 
  const double vt1 = 0.7071067811865475*vmap_v0[0] + 1.224744871391589*vmap_v0[1]*vst_1x3v_ser_p2_vel_nodes_v1[j*2 + 0]; 
  const double vt2 = 0.7071067811865475*vmap_v1[0] - 1.224744871391589*vmap_v1[1]; 
  const double vt3 = 0.7071067811865475*vmap_v2[0] + 1.224744871391589*vmap_v2[1]*vst_1x3v_ser_p2_vel_nodes_v1[j*2 + 1]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[27]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[33]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[36]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[42]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[45]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[51]; 
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
  const double om12 = vt1*p0_q + vt2*p3_q + vt3*p6_q; 
  const double om23 = vt1*p2_q + vt2*p5_q + vt3*p8_q; 
  return -(px0_q*dH_dx0*dx10*jacob_cx_inv) + -om12*dH_dv0*dv10*jacob_vx_inv + om23*dH_dv2*dv12*jacob_vz_inv; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x3v_ser_p2(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 9; ++j) alpha_quad[i*9 + j] += nc_hamil_phase_alpha_quad_vy_1x3v_ser_p2_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
