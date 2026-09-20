#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x1v_ser_p3.h> 
GKYL_CU_DH double rad_alpha_quad_vx_2x1v_ser_p3_node(int i, int j, const double *dxv, const double *rad) 
{ 
  const double *rad_vx = &rad[0]; 
  double rad_quad = 0.0; 
  for (int b = 0; b < 4; ++b) rad_quad += vst_2x1v_ser_p3_vel_ev_v0[j*4 + b]*rad_vx[b]; 
  return rad_quad; 
} 

GKYL_CU_DH void rad_alpha_quad_vx_2x1v_ser_p3(const double *dxv,
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 1; ++j) alpha_quad[i*1 + j] += rad_alpha_quad_vx_2x1v_ser_p3_node(i, j, dxv, rad); 
  } 
} 
