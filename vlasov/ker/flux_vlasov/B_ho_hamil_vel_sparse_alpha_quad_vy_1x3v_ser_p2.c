#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvx = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_1x3v_ser_p2_ho_vel_sparse_idx[s]; 
    dH_dvx += vst_1x3v_ser_p2_ho_vel_dv0_v1[j*20 + b]*hamil[b]; 
  } 
  double dH_dvz = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_1x3v_ser_p2_ho_vel_sparse_idx[s]; 
    dH_dvz += vst_1x3v_ser_p2_ho_vel_dv2_v1[j*20 + b]*hamil[b]; 
  } 
  double dv10 = 2.0/dxv[1]; 
  double dv12 = 2.0/dxv[3]; 
  const double *Bx = &qmem[9]; 
  const double *Bz = &qmem[15]; 
  double Bx_quad = 0.0; 
  double Bz_quad = 0.0; 
  for (int a = 0; a < 3; ++a) { 
    Bx_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*Bx[a]; 
    Bz_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*Bz[a]; 
  } 
  return dv12*dH_dvz*Bx_quad/jacob_vel_surf[8] - dv10*dH_dvx*Bz_quad/jacob_vel_surf[0]; 
} 

GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
