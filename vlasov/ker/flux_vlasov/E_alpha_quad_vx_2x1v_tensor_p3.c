#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x1v_tensor_p3.h> 
GKYL_CU_DH double E_alpha_quad_vx_2x1v_tensor_p3_node(int i, int j, const double *dxv, const double *qmem) 
{ 
  const double *Ex = &qmem[0]; 
  double force_quad = 0.0; 
  for (int a = 0; a < 16; ++a) force_quad += vst_2x1v_tensor_p3_conf_ev[i*16 + a]*Ex[a]; 
  return force_quad; 
} 

GKYL_CU_DH void E_alpha_quad_vx_2x1v_tensor_p3(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 25; ++i) { 
    for (int j = 0; j < 1; ++j) alpha_quad[i*1 + j] += E_alpha_quad_vx_2x1v_tensor_p3_node(i, j, dxv, qmem); 
  } 
} 
