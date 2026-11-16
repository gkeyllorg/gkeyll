#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_tensor_p1.h> 
GKYL_CU_DH int E_alpha_quad_vy_1x2v_tensor_p1_shared(int tid, int off, const double *dxv, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*2; 
  I += off*3; 
  const double *Ey = &qmem[2]; 
  if (tid < 2) { 
    const int i = tid; 
    double force_quad = 0.0; 
    for (int a = 0; a < 2; ++a) force_quad += vst_1x2v_tensor_p1_conf_ev[i*2 + a]*Ey[a]; 
    O[i] = force_quad; 
  } 
  if (tid < 3) { 
    const int j = tid; 
    I[j] = 1.0; 
  } 
  return 1; 
} 

GKYL_CU_DH void E_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[2]; 
  double I[3]; 
  for (int tid = 0; tid < 3; ++tid) E_alpha_quad_vy_1x2v_tensor_p1_shared(tid, 0, dxv, qmem, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*2 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
