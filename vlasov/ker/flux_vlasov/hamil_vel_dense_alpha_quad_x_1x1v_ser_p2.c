#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x1v_ser_p2.h> 
GKYL_CU_DH int hamil_vel_dense_alpha_quad_x_1x1v_ser_p2_shared(int tid, int nthreads, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double dv10 = 2.0/dxv[1]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  for (int i = tid; i < 1; i += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 3; ++a) P0 += vst_1x1v_ser_p2_confsurf_x0_ev_r[i*3 + a]*poisson_tensor_conf_0[a]; 
      O[0*1 + i] = P0; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 3; ++a) P0 += vst_1x1v_ser_p2_confsurf_x0_ev_l[i*3 + a]*poisson_tensor_conf_0[a]; 
      O[0*1 + i] = P0; 
    } 
  } 
  for (int j = tid; j < 3; j += nthreads) { 
    double dH_dv0 = 0.0; 
    for (int b = 0; b < 3; ++b) dH_dv0 += vst_1x1v_ser_p2_vel_vol_dv0[j*3 + b]*hamil[b]; 
    I[0*3 + j] = dH_dv0*dv10*jacob_vx_inv; 
  } 
  return 1; 
} 

GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x1v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[1]; 
  double I[3]; 
  hamil_vel_dense_alpha_quad_x_1x1v_ser_p2_shared(0, 1, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 1; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*1 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
