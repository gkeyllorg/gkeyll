#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double rad_alpha_quad_vz_3x3v_tensor_p1_node(int i, int j, const double *dxv, const double *rad) 
{ 
  const double *rad_vz = &rad[54]; 
  double rad_quad = 0.0; 
  for (int b = 0; b < 27; ++b) rad_quad += vst_3x3v_tensor_p1_vel_ev_v2[j*27 + b]*rad_vz[b]; 
  return rad_quad; 
} 

GKYL_CU_DH void rad_alpha_quad_vz_3x3v_tensor_p1(const double *dxv,
  const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 9; ++j) alpha_quad[i*9 + j] += rad_alpha_quad_vz_3x3v_tensor_p1_node(i, j, dxv, rad); 
  } 
} 
