#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x1v_tensor_p1.h> 
GKYL_CU_DH int rad_ho_alpha_quad_vx_1x1v_tensor_p1_shared(int tid, int off, const double *dxv, const double *rad,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*2; 
  I += off*1; 
  const double *rad_vx = &rad[0]; 
  if (tid < 2) { 
    const int i = tid; 
    O[i] = 1.0; 
  } 
  if (tid < 1) { 
    const int j = tid; 
    double rad_quad = 0.0; 
    for (int b = 0; b < 3; ++b) rad_quad += vst_1x1v_tensor_p1_ho_vel_ev_v0[j*3 + b]*rad_vx[b]; 
    I[j] = rad_quad; 
  } 
  return 1; 
} 

GKYL_CU_DH void rad_ho_alpha_quad_vx_1x1v_tensor_p1(const double *dxv, const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[2]; 
  double I[1]; 
  for (int tid = 0; tid < 2; ++tid) rad_ho_alpha_quad_vx_1x1v_tensor_p1_shared(tid, 0, dxv, rad, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*2 + i]*I[t*1 + j]; 
      alpha_quad[i*1 + j] += alpha; 
    } 
  } 
} 
