#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_ser_p2.h> 
GKYL_CU_DH int B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double *Bz = &qmem[15]; 
  for (int i = tid; i < 3; i += nthreads) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 3; ++a) Bz_quad += vst_1x2v_ser_p2_conf_ev[i*3 + a]*Bz[a]; 
    O[0*3 + i] = -Bz_quad; 
  } 
  for (int j = tid; j < 3; j += nthreads) { 
    double dH_dvx = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dvx += vst_1x2v_ser_p2_vel_dv0_v1[j*8 + b]*hamil[b]; 
    I[0*3 + j] = 2.0/(dxv[1]*jacob_vel_surf[0])*dH_dvx; 
  } 
  return 1; 
} 

GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[3]; 
  double I[3]; 
  B_hamil_vel_dense_alpha_quad_vy_1x2v_ser_p2_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*3 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
