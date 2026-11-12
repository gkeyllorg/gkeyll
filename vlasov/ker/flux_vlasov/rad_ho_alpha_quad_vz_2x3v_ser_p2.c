#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_ser_p2.h> 
GKYL_CU_DH int rad_ho_alpha_quad_vz_2x3v_ser_p2_shared(int tid, int off, const double *dxv, const double *rad,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*16; 
  I += off*16; 
  const double *rad_vz = &rad[40]; 
  if (tid < 16) { 
    const int i = tid; 
    O[i] = 1.0; 
  } 
  if (tid < 16) { 
    const int j = tid; 
    double rad_quad = 0.0; 
    for (int b = 0; b < 20; ++b) rad_quad += vst_2x3v_ser_p2_ho_vel_ev_v2[j*20 + b]*rad_vz[b]; 
    I[j] = rad_quad; 
  } 
  return 1; 
} 

GKYL_CU_DH void rad_ho_alpha_quad_vz_2x3v_ser_p2(const double *dxv, const double *rad, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[16]; 
  double I[16]; 
  for (int tid = 0; tid < 16; ++tid) rad_ho_alpha_quad_vz_2x3v_ser_p2_shared(tid, 0, dxv, rad, O, I); 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*16 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
