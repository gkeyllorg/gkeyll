#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x1v_ser_p3.h> 
GKYL_CU_DH int E_alpha_quad_vx_1x1v_ser_p3_shared(int tid, int off, const double *dxv, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*4; 
  I += off*1; 
  const double *Ex = &qmem[0]; 
  if (tid < 4) { 
    const int i = tid; 
    double force_quad = 0.0; 
    for (int a = 0; a < 4; ++a) force_quad += vst_1x1v_ser_p3_conf_ev[i*4 + a]*Ex[a]; 
    O[i] = force_quad; 
  } 
  if (tid < 1) { 
    const int j = tid; 
    I[j] = 1.0; 
  } 
  return 1; 
} 

GKYL_CU_DH void E_alpha_quad_vx_1x1v_ser_p3(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[4]; 
  double I[1]; 
  for (int tid = 0; tid < 4; ++tid) E_alpha_quad_vx_1x1v_ser_p3_shared(tid, 0, dxv, qmem, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*4 + i]*I[t*1 + j]; 
      alpha_quad[i*1 + j] += alpha; 
    } 
  } 
} 
