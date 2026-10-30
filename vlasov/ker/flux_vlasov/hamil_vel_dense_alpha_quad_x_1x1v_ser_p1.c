#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x1v_ser_p1.h> 
GKYL_CU_DH int hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*1; 
  I += off*2; 
  const double dv10 = 2.0/dxv[1]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  if (tid < 1) { 
    const int i = tid; 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 2; ++a) P0 += vst_1x1v_ser_p1_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_0[a]; 
      O[0*1 + i] = P0; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 2; ++a) P0 += vst_1x1v_ser_p1_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_0[a]; 
      O[0*1 + i] = P0; 
    } 
  } 
  if (tid < 2) { 
    const int j = tid; 
    double dH_dv0 = 0.0; 
    for (int b = 0; b < 2; ++b) dH_dv0 += vst_1x1v_ser_p1_vel_vol_dv0[j*2 + b]*hamil[b]; 
    I[0*2 + j] = dH_dv0*dv10*jacob_vx_inv; 
  } 
  return 1; 
} 

GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x1v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[1]; 
  double I[2]; 
  for (int tid = 0; tid < 2; ++tid) hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_shared(tid, 0, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 1; ++i) { 
    for (int j = 0; j < 2; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*1 + i]*I[t*2 + j]; 
      alpha_quad[i*2 + j] += alpha; 
    } 
  } 
} 
