#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH double hamil_vel_sparse_alpha_quad_x_2x2v_tensor_p1_node(int i, int m, int hamil_pt_edge, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[4]; 
  double dH_dv0 = 0.0; 
  for (int s = 0; s < 5; ++s) { 
    const int b = vst_2x2v_tensor_p1_vel_sparse_idx[s]; 
    dH_dv0 += vst_2x2v_tensor_p1_vel_vol_dv0[m*9 + b]*hamil[b]; 
  } 
  double dH_dv1 = 0.0; 
  for (int s = 0; s < 5; ++s) { 
    const int b = vst_2x2v_tensor_p1_vel_sparse_idx[s]; 
    dH_dv1 += vst_2x2v_tensor_p1_vel_vol_dv1[m*9 + b]*hamil[b]; 
  } 
  if (hamil_pt_edge == -1) { 
    double P0 = 0.0; 
    for (int a = 0; a < 4; ++a) P0 += vst_2x2v_tensor_p1_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 4; ++a) P1 += vst_2x2v_tensor_p1_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_1[a]; 
    return P0*dH_dv0*dv10*(1.0/jacob_vel_surf[0 + m/3]) + P1*dH_dv1*dv11*(1.0/jacob_vel_surf[3 + m%3]); 
  } 
  else if (hamil_pt_edge == 1) { 
    double P0 = 0.0; 
    for (int a = 0; a < 4; ++a) P0 += vst_2x2v_tensor_p1_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 4; ++a) P1 += vst_2x2v_tensor_p1_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_1[a]; 
    return P0*dH_dv0*dv10*(1.0/jacob_vel_surf[0 + m/3]) + P1*dH_dv1*dv11*(1.0/jacob_vel_surf[3 + m%3]); 
  } 
  return 0.0; 
} 

GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 2; ++i) { 
    for (int m = 0; m < 9; ++m) alpha_quad[i*9 + m] += hamil_vel_sparse_alpha_quad_x_2x2v_tensor_p1_node(i, m, hamil_pt_edge, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
