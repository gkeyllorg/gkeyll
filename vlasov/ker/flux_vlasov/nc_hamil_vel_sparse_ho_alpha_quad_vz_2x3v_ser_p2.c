#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_ser_p2.h> 
GKYL_CU_DH double nc_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_ser_p2_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dH_dv0 = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_2x3v_ser_p2_ho_vel_sparse_idx[s]; 
    dH_dv0 += vst_2x3v_ser_p2_ho_vel_dv0_v2[j*20 + b]*hamil[b]; 
  } 
  double dH_dv1 = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_2x3v_ser_p2_ho_vel_sparse_idx[s]; 
    dH_dv1 += vst_2x3v_ser_p2_ho_vel_dv1_v2[j*20 + b]*hamil[b]; 
  } 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[4]; 
  const double vt1 = 0.7071067811865475*vmap_v0[0] + 1.224744871391589*vmap_v0[1]*vst_2x3v_ser_p2_ho_vel_nodes_v2[j*2 + 0]; 
  const double vt2 = 0.7071067811865475*vmap_v1[0] + 1.224744871391589*vmap_v1[1]*vst_2x3v_ser_p2_ho_vel_nodes_v2[j*2 + 1]; 
  const double vt3 = 0.7071067811865475*vmap_v2[0] - 1.224744871391589*vmap_v2[1]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[80]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[88]; 
  const double *poisson_tensor_conf_4 = &poisson_tensor_conf[104]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[112]; 
  const double *poisson_tensor_conf_7 = &poisson_tensor_conf[128]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[136]; 
  double p1_q = 0.0; 
  double p2_q = 0.0; 
  double p4_q = 0.0; 
  double p5_q = 0.0; 
  double p7_q = 0.0; 
  double p8_q = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    p1_q += vst_2x3v_ser_p2_ho_conf_ev[i*8 + a]*poisson_tensor_conf_1[a]; 
    p2_q += vst_2x3v_ser_p2_ho_conf_ev[i*8 + a]*poisson_tensor_conf_2[a]; 
    p4_q += vst_2x3v_ser_p2_ho_conf_ev[i*8 + a]*poisson_tensor_conf_4[a]; 
    p5_q += vst_2x3v_ser_p2_ho_conf_ev[i*8 + a]*poisson_tensor_conf_5[a]; 
    p7_q += vst_2x3v_ser_p2_ho_conf_ev[i*8 + a]*poisson_tensor_conf_7[a]; 
    p8_q += vst_2x3v_ser_p2_ho_conf_ev[i*8 + a]*poisson_tensor_conf_8[a]; 
  } 
  const double om13 = vt1*p1_q + vt2*p4_q + vt3*p7_q; 
  const double om23 = vt1*p2_q + vt2*p5_q + vt3*p8_q; 
  return -om13*dH_dv0*dv10*jacob_vx_inv + -om23*dH_dv1*dv11*jacob_vy_inv; 
} 

GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_ser_p2(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += nc_hamil_vel_sparse_ho_alpha_quad_vz_2x3v_ser_p2_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
