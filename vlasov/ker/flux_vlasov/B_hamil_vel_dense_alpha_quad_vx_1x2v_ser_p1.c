#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_ser_p1.h> 
GKYL_CU_DH int B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double *Bz = &qmem[10]; 
  for (int i = tid; i < 2; i += nthreads) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 2; ++a) Bz_quad += vst_1x2v_ser_p1_conf_ev[i*2 + a]*Bz[a]; 
    O[0*2 + i] = Bz_quad; 
  } 
  for (int j = tid; j < 2; j += nthreads) { 
    double dH_dvy = 0.0; 
    for (int b = 0; b < 4; ++b) dH_dvy += vst_1x2v_ser_p1_vel_dv1_v0[j*4 + b]*hamil[b]; 
    I[0*2 + j] = 2.0/(dxv[2]*jacob_vel_surf[3])*dH_dvy; 
  } 
  return 1; 
} 

GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[2]; 
  double I[2]; 
  B_hamil_vel_dense_alpha_quad_vx_1x2v_ser_p1_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 2; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*2 + i]*I[t*2 + j]; 
      alpha_quad[i*2 + j] += alpha; 
    } 
  } 
} 
