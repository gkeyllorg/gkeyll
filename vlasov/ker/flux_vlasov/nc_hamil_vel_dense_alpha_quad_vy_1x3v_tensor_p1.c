#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_tensor_p1.h> 
GKYL_CU_DH double nc_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dH_dv0 = 0.0; 
  for (int b = 0; b < 27; ++b) dH_dv0 += vst_1x3v_tensor_p1_vel_dv0_v1[j*27 + b]*hamil[b]; 
  double dH_dv2 = 0.0; 
  for (int b = 0; b < 27; ++b) dH_dv2 += vst_1x3v_tensor_p1_vel_dv2_v1[j*27 + b]*hamil[b]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv12 = 2.0/dxv[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0 + j/3]; 
  const double jacob_vz_inv = 1.0/jacob_vel_surf[6 + j%3]; 
  const double xn0 = vst_1x3v_tensor_p1_vel_nodes_v1[j*2 + 0]; 
  const double xn0_sq = xn0*xn0; 
  const double vt1 = 4.677071733467426*vmap_v0[3]*xn0*xn0_sq+2.371708245126284*vmap_v0[2]*xn0_sq-2.806243040080455*vmap_v0[3]*xn0+1.224744871391589*vmap_v0[1]*xn0-0.7905694150420947*vmap_v0[2]+0.7071067811865475*vmap_v0[0]; 
  const double vt2 = -(1.8708286933869707*vmap_v1[3])+1.5811388300841895*vmap_v1[2]-1.224744871391589*vmap_v1[1]+0.7071067811865475*vmap_v1[0]; 
  const double xn2 = vst_1x3v_tensor_p1_vel_nodes_v1[j*2 + 1]; 
  const double xn2_sq = xn2*xn2; 
  const double vt3 = 4.677071733467426*vmap_v2[3]*xn2*xn2_sq+2.371708245126284*vmap_v2[2]*xn2_sq-2.806243040080455*vmap_v2[3]*xn2+1.224744871391589*vmap_v2[1]*xn2-0.7905694150420947*vmap_v2[2]+0.7071067811865475*vmap_v2[0]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[18]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[22]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[24]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[28]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[30]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[34]; 
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
  const double om12 = vt1*p0_q + vt2*p3_q + vt3*p6_q; 
  const double om23 = vt1*p2_q + vt2*p5_q + vt3*p8_q; 
  return -om12*dH_dv0*dv10*jacob_vx_inv + om23*dH_dv2*dv12*jacob_vz_inv; 
} 

GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 9; ++j) alpha_quad[i*9 + j] += nc_hamil_vel_dense_alpha_quad_vy_1x3v_tensor_p1_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
