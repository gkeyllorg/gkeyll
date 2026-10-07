#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH int rad_ho_alpha_quad_vx_2x3v_tensor_p1_shared(int tid, int nthreads, const double *dxv, const double *rad,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double *rad_vx = &rad[0]; 
  for (int i = tid; i < 4; i += nthreads) { 
    O[i] = 1.0; 
  } 
  for (int j = tid; j < 16; j += nthreads) { 
    double rad_quad = 0.0; 
    for (int b = 0; b < 27; ++b) rad_quad += vst_2x3v_tensor_p1_ho_vel_ev_v0[j*27 + b]*rad_vx[b]; 
    I[j] = rad_quad; 
  } 
  return 1; 
} 

GKYL_CU_DH void rad_ho_alpha_quad_vx_2x3v_tensor_p1(const double *dxv, const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[4]; 
  double I[16]; 
  rad_ho_alpha_quad_vx_2x3v_tensor_p1_shared(0, 1, dxv, rad, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*4 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
