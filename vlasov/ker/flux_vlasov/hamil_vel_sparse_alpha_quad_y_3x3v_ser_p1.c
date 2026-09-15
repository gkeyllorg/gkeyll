#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_ser_p1.h> 
GKYL_CU_DH double hamil_vel_sparse_alpha_quad_y_3x3v_ser_p1_node(int i, int m, int hamil_pt_edge, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  const double dv10 = 2.0/dxv[3]; 
  const double dv11 = 2.0/dxv[4]; 
  const double dv12 = 2.0/dxv[5]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[3]; 
  const double jacob_vz_inv = 1.0/jacob_vel_surf[6]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[24]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[32]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[40]; 
  double dH_dv0 = 0.0; 
  for (int s = 0; s < 4; ++s) { 
    const int b = vst_3x3v_ser_p1_vel_sparse_idx[s]; 
    dH_dv0 += vst_3x3v_ser_p1_vel_vol_dv0[m*8 + b]*hamil[b]; 
  } 
  double dH_dv1 = 0.0; 
  for (int s = 0; s < 4; ++s) { 
    const int b = vst_3x3v_ser_p1_vel_sparse_idx[s]; 
    dH_dv1 += vst_3x3v_ser_p1_vel_vol_dv1[m*8 + b]*hamil[b]; 
  } 
  double dH_dv2 = 0.0; 
  for (int s = 0; s < 4; ++s) { 
    const int b = vst_3x3v_ser_p1_vel_sparse_idx[s]; 
    dH_dv2 += vst_3x3v_ser_p1_vel_vol_dv2[m*8 + b]*hamil[b]; 
  } 
  if (hamil_pt_edge == -1) { 
    double P0 = 0.0; 
    for (int a = 0; a < 8; ++a) P0 += vst_3x3v_ser_p1_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 8; ++a) P1 += vst_3x3v_ser_p1_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_1[a]; 
    double P2 = 0.0; 
    for (int a = 0; a < 8; ++a) P2 += vst_3x3v_ser_p1_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_2[a]; 
    return P0*dH_dv0*dv10*jacob_vx_inv + P1*dH_dv1*dv11*jacob_vy_inv + P2*dH_dv2*dv12*jacob_vz_inv; 
  } 
  else if (hamil_pt_edge == 1) { 
    double P0 = 0.0; 
    for (int a = 0; a < 8; ++a) P0 += vst_3x3v_ser_p1_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 8; ++a) P1 += vst_3x3v_ser_p1_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_1[a]; 
    double P2 = 0.0; 
    for (int a = 0; a < 8; ++a) P2 += vst_3x3v_ser_p1_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_2[a]; 
    return P0*dH_dv0*dv10*jacob_vx_inv + P1*dH_dv1*dv11*jacob_vy_inv + P2*dH_dv2*dv12*jacob_vz_inv; 
  } 
  return 0.0; 
} 

GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_3x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int m = 0; m < 8; ++m) alpha_quad[i*8 + m] += hamil_vel_sparse_alpha_quad_y_3x3v_ser_p1_node(i, m, hamil_pt_edge, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
