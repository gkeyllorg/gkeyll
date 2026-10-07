#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH int rad_alpha_quad_vz_3x3v_tensor_p1_shared(int tid, int nthreads, const double *dxv, const double *rad,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double *rad_vz = &rad[54]; 
  for (int i = tid; i < 8; i += nthreads) { 
    O[i] = 1.0; 
  } 
  for (int j = tid; j < 9; j += nthreads) { 
    double rad_quad = 0.0; 
    for (int b = 0; b < 27; ++b) rad_quad += vst_3x3v_tensor_p1_vel_ev_v2[j*27 + b]*rad_vz[b]; 
    I[j] = rad_quad; 
  } 
  return 1; 
} 

GKYL_CU_DH void rad_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[9]; 
  rad_alpha_quad_vz_3x3v_tensor_p1_shared(0, 1, dxv, rad, O, I); 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*8 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
