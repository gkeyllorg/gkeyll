#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH int hamil_vel_sparse_alpha_quad_y_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 3; 
  O += off*4; 
  I += off*27; 
  const double dv10 = 2.0/dxv[3]; 
  const double dv11 = 2.0/dxv[4]; 
  const double dv12 = 2.0/dxv[5]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[24]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[32]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[40]; 
  if (tid < 4) { 
    const int i = tid; 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_3x3v_tensor_p1_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_0[a]; 
      O[0*4 + i] = P0; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_3x3v_tensor_p1_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_1[a]; 
      O[1*4 + i] = P1; 
      double P2 = 0.0; 
      for (int a = 0; a < 8; ++a) P2 += vst_3x3v_tensor_p1_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_2[a]; 
      O[2*4 + i] = P2; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_3x3v_tensor_p1_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_0[a]; 
      O[0*4 + i] = P0; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_3x3v_tensor_p1_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_1[a]; 
      O[1*4 + i] = P1; 
      double P2 = 0.0; 
      for (int a = 0; a < 8; ++a) P2 += vst_3x3v_tensor_p1_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_2[a]; 
      O[2*4 + i] = P2; 
    } 
  } 
  if (tid < 27) { 
    const int j = tid; 
    double dH_dv0 = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_3x3v_tensor_p1_vel_sparse_idx[s]; 
      dH_dv0 += vst_3x3v_tensor_p1_vel_vol_dv0[j*27 + b]*hamil[b]; 
    } 
    I[0*27 + j] = dH_dv0*dv10*(1.0/jacob_vel_surf[0 + j/9]); 
    double dH_dv1 = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_3x3v_tensor_p1_vel_sparse_idx[s]; 
      dH_dv1 += vst_3x3v_tensor_p1_vel_vol_dv1[j*27 + b]*hamil[b]; 
    } 
    I[1*27 + j] = dH_dv1*dv11*(1.0/jacob_vel_surf[3 + j/3%3]); 
    double dH_dv2 = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_3x3v_tensor_p1_vel_sparse_idx[s]; 
      dH_dv2 += vst_3x3v_tensor_p1_vel_vol_dv2[j*27 + b]*hamil[b]; 
    } 
    I[2*27 + j] = dH_dv2*dv12*(1.0/jacob_vel_surf[6 + j%3]); 
  } 
  return 3; 
} 

GKYL_CU_DH void hamil_vel_sparse_alpha_quad_y_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[12]; 
  double I[81]; 
  for (int tid = 0; tid < 27; ++tid) hamil_vel_sparse_alpha_quad_y_3x3v_tensor_p1_shared(tid, 0, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 27; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 3; ++t) alpha += O[t*4 + i]*I[t*27 + j]; 
      alpha_quad[i*27 + j] += alpha; 
    } 
  } 
} 
