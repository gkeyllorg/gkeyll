#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH double nc_hamil_phase_ho_alpha_quad_vz_2x3v_tensor_p1_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double G[4]; 
  double Gd0[4]; 
  double Gd1[4]; 
  for (int a = 0; a < 4; ++a) { G[a] = 0.0; Gd0[a] = 0.0; Gd1[a] = 0.0; } 
  for (int k = 0; k < 108; ++k) { 
    const int a = vst_2x3v_tensor_p1_ho_ph_v2_cmap[k]; 
    G[a] += vst_2x3v_tensor_p1_ho_ph_v2_V[j*27 + vst_2x3v_tensor_p1_ho_ph_v2_vrmap[k]]*(vst_2x3v_tensor_p1_ho_ph_v2_coefr[k]*hamil[k]); 
    Gd0[a] += vst_2x3v_tensor_p1_ho_ph_v2_Vd0[j*19 + vst_2x3v_tensor_p1_ho_ph_v2_vrd0map[k]]*(vst_2x3v_tensor_p1_ho_ph_v2_dcoefr0[k]*hamil[k]); 
    Gd1[a] += vst_2x3v_tensor_p1_ho_ph_v2_Vd1[j*19 + vst_2x3v_tensor_p1_ho_ph_v2_vrd1map[k]]*(vst_2x3v_tensor_p1_ho_ph_v2_dcoefr1[k]*hamil[k]); 
  } 
  double dH_dv0 = 0.0; 
  for (int a = 0; a < 4; ++a) dH_dv0 += vst_2x3v_tensor_p1_ho_ph_v2_Cm[i*4 + a]*Gd0[a]; 
  double dH_dv1 = 0.0; 
  for (int a = 0; a < 4; ++a) dH_dv1 += vst_2x3v_tensor_p1_ho_ph_v2_Cm[i*4 + a]*Gd1[a]; 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  const double *poisson_tensor_conf_x2 = &poisson_tensor_conf[8]; 
  const double *poisson_tensor_conf_x5 = &poisson_tensor_conf[20]; 
  double px0_q = 0.0; 
  double px1_q = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    px0_q += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*poisson_tensor_conf_x2[a]; 
    px1_q += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*poisson_tensor_conf_x5[a]; 
  } 
  double dH_dx0 = 0.0; 
  for (int a = 0; a < 4; ++a) dH_dx0 += vst_2x3v_tensor_p1_ho_ph_v2_CmDx0[i*4 + a]*G[a]; 
  double dH_dx1 = 0.0; 
  for (int a = 0; a < 4; ++a) dH_dx1 += vst_2x3v_tensor_p1_ho_ph_v2_CmDx1[i*4 + a]*G[a]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0 + j/4]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[4 + j%4]; 
  const double xn0 = vst_2x3v_tensor_p1_ho_vel_nodes_v2[j*2 + 0]; 
  const double xn0_sq = xn0*xn0; 
  const double vt1 = 4.677071733467426*vmap_v0[3]*xn0*xn0_sq+2.371708245126284*vmap_v0[2]*xn0_sq-2.806243040080455*vmap_v0[3]*xn0+1.224744871391589*vmap_v0[1]*xn0-0.7905694150420947*vmap_v0[2]+0.7071067811865475*vmap_v0[0]; 
  const double xn1 = vst_2x3v_tensor_p1_ho_vel_nodes_v2[j*2 + 1]; 
  const double xn1_sq = xn1*xn1; 
  const double vt2 = 4.677071733467426*vmap_v1[3]*xn1*xn1_sq+2.371708245126284*vmap_v1[2]*xn1_sq-2.806243040080455*vmap_v1[3]*xn1+1.224744871391589*vmap_v1[1]*xn1-0.7905694150420947*vmap_v1[2]+0.7071067811865475*vmap_v1[0]; 
  const double vt3 = -(1.8708286933869707*vmap_v2[3])+1.5811388300841895*vmap_v2[2]-1.224744871391589*vmap_v2[1]+0.7071067811865475*vmap_v2[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[40]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[44]; 
  const double *poisson_tensor_conf_4 = &poisson_tensor_conf[52]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[56]; 
  const double *poisson_tensor_conf_7 = &poisson_tensor_conf[64]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[68]; 
  double p1_q = 0.0; 
  double p2_q = 0.0; 
  double p4_q = 0.0; 
  double p5_q = 0.0; 
  double p7_q = 0.0; 
  double p8_q = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    p1_q += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*poisson_tensor_conf_1[a]; 
    p2_q += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*poisson_tensor_conf_2[a]; 
    p4_q += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*poisson_tensor_conf_4[a]; 
    p5_q += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*poisson_tensor_conf_5[a]; 
    p7_q += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*poisson_tensor_conf_7[a]; 
    p8_q += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*poisson_tensor_conf_8[a]; 
  } 
  const double om13 = vt1*p1_q + vt2*p4_q + vt3*p7_q; 
  const double om23 = vt1*p2_q + vt2*p5_q + vt3*p8_q; 
  return -(px0_q*dH_dx0*dx10*jacob_cx_inv + px1_q*dH_dx1*dx11*jacob_cy_inv) + -om13*dH_dv0*dv10*jacob_vx_inv + -om23*dH_dv1*dv11*jacob_vy_inv; 
} 

GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vz_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += nc_hamil_phase_ho_alpha_quad_vz_2x3v_tensor_p1_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
