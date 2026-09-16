#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double E_ho_alpha_quad_vy_1x3v_ser_p2_node(int i, int j, const double *dxv, const double *qmem) 
{ 
  const double *Ey = &qmem[3]; 
  double force_quad = 0.0; 
  for (int a = 0; a < 3; ++a) force_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*Ey[a]; 
  return force_quad; 
} 

GKYL_CU_DH void E_ho_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += E_ho_alpha_quad_vy_1x3v_ser_p2_node(i, j, dxv, qmem); 
  } 
} 
