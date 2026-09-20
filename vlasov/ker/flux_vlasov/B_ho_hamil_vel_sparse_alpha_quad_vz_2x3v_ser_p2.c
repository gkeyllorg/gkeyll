#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_ser_p2.h> 
GKYL_CU_DH double B_ho_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvx = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_2x3v_ser_p2_ho_vel_sparse_idx[s]; 
    dH_dvx += vst_2x3v_ser_p2_ho_vel_dv0_v2[j*20 + b]*hamil[b]; 
  } 
  double dH_dvy = 0.0; 
  for (int s = 0; s < 7; ++s) { 
    const int b = vst_2x3v_ser_p2_ho_vel_sparse_idx[s]; 
    dH_dvy += vst_2x3v_ser_p2_ho_vel_dv1_v2[j*20 + b]*hamil[b]; 
  } 
  double dv10 = 2.0/dxv[2]; 
  double dv11 = 2.0/dxv[3]; 
  const double *Bx = &qmem[24]; 
  const double *By = &qmem[32]; 
  double Bx_quad = 0.0; 
  double By_quad = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    Bx_quad += vst_2x3v_ser_p2_ho_conf_ev[i*8 + a]*Bx[a]; 
    By_quad += vst_2x3v_ser_p2_ho_conf_ev[i*8 + a]*By[a]; 
  } 
  return dv10*dH_dvx*By_quad/jacob_vel_surf[0] - dv11*dH_dvy*Bx_quad/jacob_vel_surf[4]; 
} 

GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += B_ho_hamil_vel_sparse_alpha_quad_vz_2x3v_ser_p2_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
