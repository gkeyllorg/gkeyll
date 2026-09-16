#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double B_ho_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvy = 0.0; 
  for (int b = 0; b < 20; ++b) dH_dvy += vst_1x3v_ser_p2_ho_vel_dv1_v0[j*20 + b]*hamil[b]; 
  double dH_dvz = 0.0; 
  for (int b = 0; b < 20; ++b) dH_dvz += vst_1x3v_ser_p2_ho_vel_dv2_v0[j*20 + b]*hamil[b]; 
  double dv11 = 2.0/dxv[2]; 
  double dv12 = 2.0/dxv[3]; 
  const double *By = &qmem[12]; 
  const double *Bz = &qmem[15]; 
  double By_quad = 0.0; 
  double Bz_quad = 0.0; 
  for (int a = 0; a < 3; ++a) { 
    By_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*By[a]; 
    Bz_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*Bz[a]; 
  } 
  return dv11*dH_dvy*Bz_quad/jacob_vel_surf[4] - dv12*dH_dvz*By_quad/jacob_vel_surf[8]; 
} 

GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += B_ho_hamil_vel_dense_alpha_quad_vx_1x3v_ser_p2_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
