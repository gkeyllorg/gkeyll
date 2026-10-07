#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p1.h> 
GKYL_CU_DH int E_alpha_quad_vx_2x2v_ser_p1_shared(int tid, int nthreads, const double *dxv, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double *Ex = &qmem[0]; 
  for (int i = tid; i < 4; i += nthreads) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 4; ++a) force_quad += vst_2x2v_ser_p1_conf_ev[i*4 + a]*Ex[a]; 
    O[i] = force_quad; 
  } 
  for (int j = tid; j < 2; j += nthreads) { 
    I[j] = 1.0; 
  } 
  return 1; 
} 

GKYL_CU_DH void E_alpha_quad_vx_2x2v_ser_p1(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[4]; 
  double I[2]; 
  E_alpha_quad_vx_2x2v_ser_p1_shared(0, 1, dxv, qmem, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 2; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*4 + i]*I[t*2 + j]; 
      alpha_quad[i*2 + j] += alpha; 
    } 
  } 
} 
