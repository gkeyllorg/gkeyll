#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p2.h> 
GKYL_CU_DH double B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p2_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvy = 0.0; 
  for (int b = 0; b < 9; ++b) dH_dvy += vst_2x2v_tensor_p2_vel_dv1_v0[j*9 + b]*hamil[b]; 
  const double *Bz = &qmem[45]; 
  double Bz_quad = 0.0; 
  for (int a = 0; a < 9; ++a) Bz_quad += vst_2x2v_tensor_p2_conf_ev[i*9 + a]*Bz[a]; 
  return 2.0/(dxv[3]*jacob_vel_surf[4 + j])*dH_dvy*Bz_quad; 
} 

GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += B_hamil_vel_dense_alpha_quad_vx_2x2v_tensor_p2_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
