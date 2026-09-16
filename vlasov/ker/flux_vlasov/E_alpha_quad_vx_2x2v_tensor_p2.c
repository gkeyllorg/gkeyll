#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p2.h> 
GKYL_CU_DH double E_alpha_quad_vx_2x2v_tensor_p2_node(int i, int j, const double *dxv, const double *qmem) 
{ 
  const double *Ex = &qmem[0]; 
  double force_quad = 0.0; 
  for (int a = 0; a < 9; ++a) force_quad += vst_2x2v_tensor_p2_conf_ev[i*9 + a]*Ex[a]; 
  return force_quad; 
} 

GKYL_CU_DH void E_alpha_quad_vx_2x2v_tensor_p2(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += E_alpha_quad_vx_2x2v_tensor_p2_node(i, j, dxv, qmem); 
  } 
} 
