#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double rad_ho_alpha_quad_vx_1x3v_ser_p2_node(int i, int j, const double *dxv, const double *rad) 
{ 
  const double *rad_vx = &rad[0]; 
  double rad_quad = 0.0; 
  for (int b = 0; b < 20; ++b) rad_quad += vst_1x3v_ser_p2_ho_vel_ev_v0[j*20 + b]*rad_vx[b]; 
  return rad_quad; 
} 

GKYL_CU_DH void rad_ho_alpha_quad_vx_1x3v_ser_p2(const double *dxv,
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += rad_ho_alpha_quad_vx_1x3v_ser_p2_node(i, j, dxv, rad); 
  } 
} 
