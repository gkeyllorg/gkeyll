#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH int E_alpha_quad_vy_1x3v_ser_p2_shared(int tid, int off, const double *dxv, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*3; 
  I += off*9; 
  const double *Ey = &qmem[3]; 
  if (tid < 3) { 
    const int i = tid; 
    double force_quad = 0.0; 
    for (int a = 0; a < 3; ++a) force_quad += vst_1x3v_ser_p2_conf_ev[i*3 + a]*Ey[a]; 
    O[i] = force_quad; 
  } 
  if (tid < 9) { 
    const int j = tid; 
    I[j] = 1.0; 
  } 
  return 1; 
} 

GKYL_CU_DH void E_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[3]; 
  double I[9]; 
  for (int tid = 0; tid < 9; ++tid) E_alpha_quad_vy_1x3v_ser_p2_shared(tid, 0, dxv, qmem, O, I); 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*3 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
