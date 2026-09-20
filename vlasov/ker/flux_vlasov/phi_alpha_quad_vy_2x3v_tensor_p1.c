#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH double phi_alpha_quad_vy_2x3v_tensor_p1_node(int i, int j, const double *dxv, const double *jacob_pos, const double *phi) 
{ 
  double dx11 = 2.0/dxv[1]; 
  const double *jacob_cy = &jacob_pos[2]; 
  const double jacob_cy_inv = 1.0/jacob_cy[0]; 
  double force_quad = 0.0; 
  for (int a = 0; a < 4; ++a) force_quad += vst_2x3v_tensor_p1_conf_dx1[i*4 + a]*phi[a]; 
  force_quad *= jacob_cy_inv; 
  return -dx11*force_quad; 
} 

GKYL_CU_DH void phi_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 9; ++j) alpha_quad[i*9 + j] += phi_alpha_quad_vy_2x3v_tensor_p1_node(i, j, dxv, jacob_pos, phi); 
  } 
} 
