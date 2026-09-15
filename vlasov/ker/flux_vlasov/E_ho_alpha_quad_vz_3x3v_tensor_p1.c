#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double E_ho_alpha_quad_vz_3x3v_tensor_p1_node(int i, int j, const double *dxv, const double *qmem) 
{ 
  const double *Ez = &qmem[16]; 
  double force_quad = 0.0; 
  for (int a = 0; a < 8; ++a) force_quad += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*Ez[a]; 
  return force_quad; 
} 

GKYL_CU_DH void E_ho_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += E_ho_alpha_quad_vz_3x3v_tensor_p1_node(i, j, dxv, qmem); 
  } 
} 
