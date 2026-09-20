#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_tensor_p2.h> 
GKYL_CU_DH double B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p2_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvx = 0.0; 
  for (int s = 0; s < 5; ++s) { 
    const int b = vst_1x2v_tensor_p2_vel_sparse_idx[s]; 
    dH_dvx += vst_1x2v_tensor_p2_vel_dv0_v1[j*9 + b]*hamil[b]; 
  } 
  const double *Bz = &qmem[15]; 
  double Bz_quad = 0.0; 
  for (int a = 0; a < 3; ++a) Bz_quad += vst_1x2v_tensor_p2_conf_ev[i*3 + a]*Bz[a]; 
  return -2.0/(dxv[1]*jacob_vel_surf[0 + j])*dH_dvx*Bz_quad; 
} 

GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += B_hamil_vel_sparse_alpha_quad_vy_1x2v_tensor_p2_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
