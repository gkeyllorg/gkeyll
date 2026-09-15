#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_tensor_p1.h> 
GKYL_CU_DH double rad_ho_alpha_quad_vy_1x2v_tensor_p1_node(int i, int j, const double *dxv, const double *rad) 
{ 
  const double *rad_vy = &rad[9]; 
  double rad_quad = 0.0; 
  for (int b = 0; b < 9; ++b) rad_quad += vst_1x2v_tensor_p1_ho_vel_ev_v1[j*9 + b]*rad_vy[b]; 
  return rad_quad; 
} 

GKYL_CU_DH void rad_ho_alpha_quad_vy_1x2v_tensor_p1(const double *dxv,
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += rad_ho_alpha_quad_vy_1x2v_tensor_p1_node(i, j, dxv, rad); 
  } 
} 
