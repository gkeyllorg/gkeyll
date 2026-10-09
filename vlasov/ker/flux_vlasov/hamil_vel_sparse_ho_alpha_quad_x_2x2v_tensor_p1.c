#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH int hamil_vel_sparse_ho_alpha_quad_x_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  O += off*2; 
  I += off*16; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[4]; 
  if (tid < 2) { 
    const int i = tid; 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 4; ++a) P0 += vst_2x2v_tensor_p1_ho_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_0[a]; 
      O[0*2 + i] = P0; 
      double P1 = 0.0; 
      for (int a = 0; a < 4; ++a) P1 += vst_2x2v_tensor_p1_ho_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_1[a]; 
      O[1*2 + i] = P1; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 4; ++a) P0 += vst_2x2v_tensor_p1_ho_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_0[a]; 
      O[0*2 + i] = P0; 
      double P1 = 0.0; 
      for (int a = 0; a < 4; ++a) P1 += vst_2x2v_tensor_p1_ho_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_1[a]; 
      O[1*2 + i] = P1; 
    } 
  } 
  if (tid < 16) { 
    const int j = tid; 
    double dH_dv0 = 0.0; 
    for (int s = 0; s < 5; ++s) { 
      const int b = vst_2x2v_tensor_p1_ho_vel_sparse_idx[s]; 
      dH_dv0 += vst_2x2v_tensor_p1_ho_vel_vol_dv0[j*9 + b]*hamil[b]; 
    } 
    I[0*16 + j] = dH_dv0*dv10*(1.0/jacob_vel_surf[0 + j/4]); 
    double dH_dv1 = 0.0; 
    for (int s = 0; s < 5; ++s) { 
      const int b = vst_2x2v_tensor_p1_ho_vel_sparse_idx[s]; 
      dH_dv1 += vst_2x2v_tensor_p1_ho_vel_vol_dv1[j*9 + b]*hamil[b]; 
    } 
    I[1*16 + j] = dH_dv1*dv11*(1.0/jacob_vel_surf[4 + j%4]); 
  } 
  return 2; 
} 

GKYL_CU_DH void hamil_vel_sparse_ho_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[4]; 
  double I[32]; 
  for (int tid = 0; tid < 16; ++tid) hamil_vel_sparse_ho_alpha_quad_x_2x2v_tensor_p1_shared(tid, 0, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*2 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
