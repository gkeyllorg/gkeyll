#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_ser_p2.h> 
GKYL_CU_DH int rad_alpha_quad_vx_1x2v_ser_p2_shared(int tid, int off, const double *dxv, const double *rad,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*3; 
  I += off*3; 
  const double *rad_vx = &rad[0]; 
  if (tid < 3) { 
    const int i = tid; 
    O[i] = 1.0; 
  } 
  if (tid < 3) { 
    const int j = tid; 
    double rad_quad = 0.0; 
    for (int b = 0; b < 8; ++b) rad_quad += vst_1x2v_ser_p2_vel_ev_v0[j*8 + b]*rad_vx[b]; 
    I[j] = rad_quad; 
  } 
  return 1; 
} 

GKYL_CU_DH void rad_alpha_quad_vx_1x2v_ser_p2(const double *dxv, const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[3]; 
  double I[3]; 
  for (int tid = 0; tid < 3; ++tid) rad_alpha_quad_vx_1x2v_ser_p2_shared(tid, 0, dxv, rad, O, I); 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*3 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
