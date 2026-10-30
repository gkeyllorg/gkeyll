#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH int rad_alpha_quad_vy_2x3v_tensor_p1_shared(int tid, int off, const double *dxv, const double *rad,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*4; 
  I += off*9; 
  const double *rad_vy = &rad[27]; 
  if (tid < 4) { 
    const int i = tid; 
    O[i] = 1.0; 
  } 
  if (tid < 9) { 
    const int j = tid; 
    double rad_quad = 0.0; 
    for (int b = 0; b < 27; ++b) rad_quad += vst_2x3v_tensor_p1_vel_ev_v1[j*27 + b]*rad_vy[b]; 
    I[j] = rad_quad; 
  } 
  return 1; 
} 

GKYL_CU_DH void rad_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[4]; 
  double I[9]; 
  for (int tid = 0; tid < 9; ++tid) rad_alpha_quad_vy_2x3v_tensor_p1_shared(tid, 0, dxv, rad, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*4 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
