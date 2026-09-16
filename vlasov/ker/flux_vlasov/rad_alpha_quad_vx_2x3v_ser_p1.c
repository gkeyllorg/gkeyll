#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_ser_p1.h> 
GKYL_CU_DH double rad_alpha_quad_vx_2x3v_ser_p1_node(int i, int j, const double *dxv, const double *rad) 
{ 
  const double *rad_vx = &rad[0]; 
  double rad_quad = 0.0; 
  for (int b = 0; b < 8; ++b) rad_quad += vst_2x3v_ser_p1_vel_ev_v0[j*8 + b]*rad_vx[b]; 
  return rad_quad; 
} 

GKYL_CU_DH void rad_alpha_quad_vx_2x3v_ser_p1(const double *dxv,
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += rad_alpha_quad_vx_2x3v_ser_p1_node(i, j, dxv, rad); 
  } 
} 
