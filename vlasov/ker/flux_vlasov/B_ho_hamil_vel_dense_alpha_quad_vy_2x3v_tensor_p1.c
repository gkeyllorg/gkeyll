#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH double B_ho_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvx = 0.0; 
  for (int b = 0; b < 27; ++b) dH_dvx += vst_2x3v_tensor_p1_ho_vel_dv0_v1[j*27 + b]*hamil[b]; 
  double dH_dvz = 0.0; 
  for (int b = 0; b < 27; ++b) dH_dvz += vst_2x3v_tensor_p1_ho_vel_dv2_v1[j*27 + b]*hamil[b]; 
  double dv10 = 2.0/dxv[2]; 
  double dv12 = 2.0/dxv[4]; 
  const double *Bx = &qmem[12]; 
  const double *Bz = &qmem[20]; 
  double Bx_quad = 0.0; 
  double Bz_quad = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    Bx_quad += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*Bx[a]; 
    Bz_quad += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*Bz[a]; 
  } 
  return dv12*dH_dvz*Bx_quad/jacob_vel_surf[8 + j%4] - dv10*dH_dvx*Bz_quad/jacob_vel_surf[0 + j/4]; 
} 

GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += B_ho_hamil_vel_dense_alpha_quad_vy_2x3v_tensor_p1_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
