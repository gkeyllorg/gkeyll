#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double phi_ho_alpha_quad_vz_3x3v_tensor_p1_node(int i, int j, const double *dxv, const double *jacob_pos, const double *phi) 
{ 
  double dx12 = 2.0/dxv[2]; 
  const double *jacob_cz = &jacob_pos[4]; 
  const double jacob_cz_inv = 1.0/jacob_cz[0]; 
  double force_quad = 0.0; 
  for (int a = 0; a < 8; ++a) force_quad += vst_3x3v_tensor_p1_ho_conf_dx2[i*8 + a]*phi[a]; 
  force_quad *= jacob_cz_inv; 
  return -dx12*force_quad; 
} 

GKYL_CU_DH void phi_ho_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *jacob_pos, const double *phi, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += phi_ho_alpha_quad_vz_3x3v_tensor_p1_node(i, j, dxv, jacob_pos, phi); 
  } 
} 
