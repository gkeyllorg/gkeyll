#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x1v_ser_p3.h> 
GKYL_CU_DH int E_ho_alpha_quad_vx_2x1v_ser_p3_shared(int tid, int nthreads, const double *dxv, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double *Ex = &qmem[0]; 
  for (int i = tid; i < 25; i += nthreads) { 
    double force_quad = 0.0; 
    for (int a = 0; a < 12; ++a) force_quad += vst_2x1v_ser_p3_ho_conf_ev[i*12 + a]*Ex[a]; 
    O[i] = force_quad; 
  } 
  for (int j = tid; j < 1; j += nthreads) { 
    I[j] = 1.0; 
  } 
  return 1; 
} 

GKYL_CU_DH void E_ho_alpha_quad_vx_2x1v_ser_p3(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[25]; 
  double I[1]; 
  E_ho_alpha_quad_vx_2x1v_ser_p3_shared(0, 1, dxv, qmem, O, I); 
  for (int i = 0; i < 25; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*25 + i]*I[t*1 + j]; 
      alpha_quad[i*1 + j] += alpha; 
    } 
  } 
} 
