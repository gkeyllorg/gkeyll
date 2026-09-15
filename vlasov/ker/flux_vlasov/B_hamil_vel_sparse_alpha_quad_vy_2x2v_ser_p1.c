#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p1.h> 
GKYL_CU_DH double B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvx = 0.0; 
  for (int s = 0; s < 3; ++s) { 
    const int b = vst_2x2v_ser_p1_vel_sparse_idx[s]; 
    dH_dvx += vst_2x2v_ser_p1_vel_dv0_v1[j*4 + b]*hamil[b]; 
  } 
  const double *Bz = &qmem[20]; 
  double Bz_quad = 0.0; 
  for (int a = 0; a < 4; ++a) Bz_quad += vst_2x2v_ser_p1_conf_ev[i*4 + a]*Bz[a]; 
  return -2.0/(dxv[2]*jacob_vel_surf[0])*dH_dvx*Bz_quad; 
} 

GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 2; ++j) alpha_quad[i*2 + j] += B_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p1_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
