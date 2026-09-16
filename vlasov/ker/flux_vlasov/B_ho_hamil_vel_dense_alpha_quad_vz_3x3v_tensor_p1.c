#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double B_ho_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvx = 0.0; 
  for (int b = 0; b < 27; ++b) dH_dvx += vst_3x3v_tensor_p1_ho_vel_dv0_v2[j*27 + b]*hamil[b]; 
  double dH_dvy = 0.0; 
  for (int b = 0; b < 27; ++b) dH_dvy += vst_3x3v_tensor_p1_ho_vel_dv1_v2[j*27 + b]*hamil[b]; 
  double dv10 = 2.0/dxv[3]; 
  double dv11 = 2.0/dxv[4]; 
  const double *Bx = &qmem[24]; 
  const double *By = &qmem[32]; 
  double Bx_quad = 0.0; 
  double By_quad = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    Bx_quad += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*Bx[a]; 
    By_quad += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*By[a]; 
  } 
  return dv10*dH_dvx*By_quad/jacob_vel_surf[0 + j/4] - dv11*dH_dvy*Bx_quad/jacob_vel_surf[4 + j%4]; 
} 

GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += B_ho_hamil_vel_dense_alpha_quad_vz_3x3v_tensor_p1_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
