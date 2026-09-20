#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_tensor_p1.h> 
GKYL_CU_DH double hamil_vel_sparse_ho_alpha_quad_x_1x3v_tensor_p1_node(int i, int m, int hamil_pt_edge, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[2]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[4]; 
  double dH_dv0 = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_1x3v_tensor_p1_ho_vel_sparse_idx[s]; 
    dH_dv0 += vst_1x3v_tensor_p1_ho_vel_vol_dv0[m*27 + b]*hamil[b]; 
  } 
  double dH_dv1 = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_1x3v_tensor_p1_ho_vel_sparse_idx[s]; 
    dH_dv1 += vst_1x3v_tensor_p1_ho_vel_vol_dv1[m*27 + b]*hamil[b]; 
  } 
  double dH_dv2 = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_1x3v_tensor_p1_ho_vel_sparse_idx[s]; 
    dH_dv2 += vst_1x3v_tensor_p1_ho_vel_vol_dv2[m*27 + b]*hamil[b]; 
  } 
  if (hamil_pt_edge == -1) { 
    double P0 = 0.0; 
    for (int a = 0; a < 2; ++a) P0 += vst_1x3v_tensor_p1_ho_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 2; ++a) P1 += vst_1x3v_tensor_p1_ho_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_1[a]; 
    double P2 = 0.0; 
    for (int a = 0; a < 2; ++a) P2 += vst_1x3v_tensor_p1_ho_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_2[a]; 
    return P0*dH_dv0*dv10*(1.0/jacob_vel_surf[0 + m/16]) + P1*dH_dv1*dv11*(1.0/jacob_vel_surf[4 + m/4%4]) + P2*dH_dv2*dv12*(1.0/jacob_vel_surf[8 + m%4]); 
  } 
  else if (hamil_pt_edge == 1) { 
    double P0 = 0.0; 
    for (int a = 0; a < 2; ++a) P0 += vst_1x3v_tensor_p1_ho_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 2; ++a) P1 += vst_1x3v_tensor_p1_ho_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_1[a]; 
    double P2 = 0.0; 
    for (int a = 0; a < 2; ++a) P2 += vst_1x3v_tensor_p1_ho_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_2[a]; 
    return P0*dH_dv0*dv10*(1.0/jacob_vel_surf[0 + m/16]) + P1*dH_dv1*dv11*(1.0/jacob_vel_surf[4 + m/4%4]) + P2*dH_dv2*dv12*(1.0/jacob_vel_surf[8 + m%4]); 
  } 
  return 0.0; 
} 

GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_x_1x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 1; ++i) { 
    for (int m = 0; m < 64; ++m) alpha_quad[i*64 + m] += hamil_vel_sparse_ho_alpha_quad_x_1x3v_tensor_p1_node(i, m, hamil_pt_edge, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
