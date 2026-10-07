#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH int hamil_vel_sparse_alpha_quad_x_2x3v_tensor_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 3; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[4]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[8]; 
  for (int i = tid; i < 2; i += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 4; ++a) P0 += vst_2x3v_tensor_p1_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_0[a]; 
      O[0*2 + i] = P0; 
      double P1 = 0.0; 
      for (int a = 0; a < 4; ++a) P1 += vst_2x3v_tensor_p1_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_1[a]; 
      O[1*2 + i] = P1; 
      double P2 = 0.0; 
      for (int a = 0; a < 4; ++a) P2 += vst_2x3v_tensor_p1_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_2[a]; 
      O[2*2 + i] = P2; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 4; ++a) P0 += vst_2x3v_tensor_p1_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_0[a]; 
      O[0*2 + i] = P0; 
      double P1 = 0.0; 
      for (int a = 0; a < 4; ++a) P1 += vst_2x3v_tensor_p1_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_1[a]; 
      O[1*2 + i] = P1; 
      double P2 = 0.0; 
      for (int a = 0; a < 4; ++a) P2 += vst_2x3v_tensor_p1_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_2[a]; 
      O[2*2 + i] = P2; 
    } 
  } 
  for (int j = tid; j < 27; j += nthreads) { 
    double dH_dv0 = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_2x3v_tensor_p1_vel_sparse_idx[s]; 
      dH_dv0 += vst_2x3v_tensor_p1_vel_vol_dv0[j*27 + b]*hamil[b]; 
    } 
    I[0*27 + j] = dH_dv0*dv10*(1.0/jacob_vel_surf[0 + j/9]); 
    double dH_dv1 = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_2x3v_tensor_p1_vel_sparse_idx[s]; 
      dH_dv1 += vst_2x3v_tensor_p1_vel_vol_dv1[j*27 + b]*hamil[b]; 
    } 
    I[1*27 + j] = dH_dv1*dv11*(1.0/jacob_vel_surf[3 + j/3%3]); 
    double dH_dv2 = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_2x3v_tensor_p1_vel_sparse_idx[s]; 
      dH_dv2 += vst_2x3v_tensor_p1_vel_vol_dv2[j*27 + b]*hamil[b]; 
    } 
    I[2*27 + j] = dH_dv2*dv12*(1.0/jacob_vel_surf[6 + j%3]); 
  } 
  return 3; 
} 

GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[6]; 
  double I[81]; 
  hamil_vel_sparse_alpha_quad_x_2x3v_tensor_p1_shared(0, 1, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 27; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 3; ++t) alpha += O[t*2 + i]*I[t*27 + j]; 
      alpha_quad[i*27 + j] += alpha; 
    } 
  } 
} 
