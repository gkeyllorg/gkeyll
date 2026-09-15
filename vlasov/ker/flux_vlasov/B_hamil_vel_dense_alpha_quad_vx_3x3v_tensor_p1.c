#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double B_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvy = 0.0; 
  for (int b = 0; b < 27; ++b) dH_dvy += vst_3x3v_tensor_p1_vel_dv1_v0[j*27 + b]*hamil[b]; 
  double dH_dvz = 0.0; 
  for (int b = 0; b < 27; ++b) dH_dvz += vst_3x3v_tensor_p1_vel_dv2_v0[j*27 + b]*hamil[b]; 
  double dv11 = 2.0/dxv[4]; 
  double dv12 = 2.0/dxv[5]; 
  const double *By = &qmem[32]; 
  const double *Bz = &qmem[40]; 
  double By_quad = 0.0; 
  double Bz_quad = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    By_quad += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*By[a]; 
    Bz_quad += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*Bz[a]; 
  } 
  return dv11*dH_dvy*Bz_quad/jacob_vel_surf[3 + j/3] - dv12*dH_dvz*By_quad/jacob_vel_surf[6 + j%3]; 
} 

GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 9; ++j) alpha_quad[i*9 + j] += B_hamil_vel_dense_alpha_quad_vx_3x3v_tensor_p1_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
