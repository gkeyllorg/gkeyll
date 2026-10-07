#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_ser_p1.h> 
GKYL_CU_DH int hamil_vel_sparse_alpha_quad_x_1x2v_ser_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[3]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[2]; 
  for (int i = tid; i < 1; i += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 2; ++a) P0 += vst_1x2v_ser_p1_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_0[a]; 
      O[0*1 + i] = P0; 
      double P1 = 0.0; 
      for (int a = 0; a < 2; ++a) P1 += vst_1x2v_ser_p1_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_1[a]; 
      O[1*1 + i] = P1; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 2; ++a) P0 += vst_1x2v_ser_p1_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_0[a]; 
      O[0*1 + i] = P0; 
      double P1 = 0.0; 
      for (int a = 0; a < 2; ++a) P1 += vst_1x2v_ser_p1_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_1[a]; 
      O[1*1 + i] = P1; 
    } 
  } 
  for (int j = tid; j < 4; j += nthreads) { 
    double dH_dv0 = 0.0; 
    for (int s = 0; s < 3; ++s) { 
      const int b = vst_1x2v_ser_p1_vel_sparse_idx[s]; 
      dH_dv0 += vst_1x2v_ser_p1_vel_vol_dv0[j*4 + b]*hamil[b]; 
    } 
    I[0*4 + j] = dH_dv0*dv10*jacob_vx_inv; 
    double dH_dv1 = 0.0; 
    for (int s = 0; s < 3; ++s) { 
      const int b = vst_1x2v_ser_p1_vel_sparse_idx[s]; 
      dH_dv1 += vst_1x2v_ser_p1_vel_vol_dv1[j*4 + b]*hamil[b]; 
    } 
    I[1*4 + j] = dH_dv1*dv11*jacob_vy_inv; 
  } 
  return 2; 
} 

GKYL_CU_DH void hamil_vel_sparse_alpha_quad_x_1x2v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[2]; 
  double I[8]; 
  hamil_vel_sparse_alpha_quad_x_1x2v_ser_p1_shared(0, 1, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 1; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*1 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
