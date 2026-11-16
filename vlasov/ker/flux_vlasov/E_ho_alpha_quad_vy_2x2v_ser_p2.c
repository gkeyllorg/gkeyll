#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p2.h> 
GKYL_CU_DH int E_ho_alpha_quad_vy_2x2v_ser_p2_shared(int tid, int off, const double *dxv, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  O += off*16; 
  I += off*4; 
  const double *Ey = &qmem[8]; 
  if (tid < 16) { 
    const int i = tid; 
    double force_quad = 0.0; 
    for (int a = 0; a < 8; ++a) force_quad += vst_2x2v_ser_p2_ho_conf_ev[i*8 + a]*Ey[a]; 
    O[i] = force_quad; 
  } 
  if (tid < 4) { 
    const int j = tid; 
    I[j] = 1.0; 
  } 
  return 1; 
} 

GKYL_CU_DH void E_ho_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[16]; 
  double I[4]; 
  for (int tid = 0; tid < 16; ++tid) E_ho_alpha_quad_vy_2x2v_ser_p2_shared(tid, 0, dxv, qmem, O, I); 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*16 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
