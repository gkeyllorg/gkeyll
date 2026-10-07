#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH int rad_alpha_quad_vz_1x3v_ser_p2_shared(int tid, int nthreads, const double *dxv, const double *rad,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double *rad_vz = &rad[40]; 
  for (int i = tid; i < 3; i += nthreads) { 
    O[i] = 1.0; 
  } 
  for (int j = tid; j < 9; j += nthreads) { 
    double rad_quad = 0.0; 
    for (int b = 0; b < 20; ++b) rad_quad += vst_1x3v_ser_p2_vel_ev_v2[j*20 + b]*rad_vz[b]; 
    I[j] = rad_quad; 
  } 
  return 1; 
} 

GKYL_CU_DH void rad_alpha_quad_vz_1x3v_ser_p2(const double *dxv, const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[3]; 
  double I[9]; 
  rad_alpha_quad_vz_1x3v_ser_p2_shared(0, 1, dxv, rad, O, I); 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*3 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
