#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p1.h> 
GKYL_CU_DH double E_alpha_quad_vz_1x3v_ser_p1_node(int i, int j, const double *dxv, const double *qmem) 
{ 
  const double *Ez = &qmem[4]; 
  double force_quad = 0.0; 
  for (int a = 0; a < 2; ++a) force_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*Ez[a]; 
  return force_quad; 
} 

GKYL_CU_DH void E_alpha_quad_vz_1x3v_ser_p1(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += E_alpha_quad_vz_1x3v_ser_p1_node(i, j, dxv, qmem); 
  } 
} 
