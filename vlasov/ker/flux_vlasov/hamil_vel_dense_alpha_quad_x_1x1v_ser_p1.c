#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x1v_ser_p1.h> 
GKYL_CU_DH double hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_node(int i, int m, int hamil_pt_edge, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  const double dv10 = 2.0/dxv[1]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  double dH_dv0 = 0.0; 
  for (int b = 0; b < 2; ++b) dH_dv0 += vst_1x1v_ser_p1_vel_vol_dv0[m*2 + b]*hamil[b]; 
  if (hamil_pt_edge == -1) { 
    double P0 = 0.0; 
    for (int a = 0; a < 2; ++a) P0 += vst_1x1v_ser_p1_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_0[a]; 
    return P0*dH_dv0*dv10*jacob_vx_inv; 
  } 
  else if (hamil_pt_edge == 1) { 
    double P0 = 0.0; 
    for (int a = 0; a < 2; ++a) P0 += vst_1x1v_ser_p1_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_0[a]; 
    return P0*dH_dv0*dv10*jacob_vx_inv; 
  } 
  return 0.0; 
} 

GKYL_CU_DH void hamil_vel_dense_alpha_quad_x_1x1v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 1; ++i) { 
    for (int m = 0; m < 2; ++m) alpha_quad[i*2 + m] += hamil_vel_dense_alpha_quad_x_1x1v_ser_p1_node(i, m, hamil_pt_edge, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
