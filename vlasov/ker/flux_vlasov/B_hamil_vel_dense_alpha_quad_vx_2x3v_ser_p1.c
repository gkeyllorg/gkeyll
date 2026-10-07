#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_ser_p1.h> 
GKYL_CU_DH int B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  const double *Bz = &qmem[20]; 
  const double *By = &qmem[16]; 
  for (int i = tid; i < 4; i += nthreads) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 4; ++a) Bz_quad += vst_2x3v_ser_p1_conf_ev[i*4 + a]*Bz[a]; 
    O[0*4 + i] = Bz_quad; 
    double By_quad = 0.0; 
    for (int a = 0; a < 4; ++a) By_quad += vst_2x3v_ser_p1_conf_ev[i*4 + a]*By[a]; 
    O[1*4 + i] = -By_quad; 
  } 
  for (int j = tid; j < 4; j += nthreads) { 
    double dH_dvy = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dvy += vst_2x3v_ser_p1_vel_dv1_v0[j*8 + b]*hamil[b]; 
    I[0*4 + j] = 2.0/(dxv[3]*jacob_vel_surf[3])*dH_dvy; 
    double dH_dvz = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dvz += vst_2x3v_ser_p1_vel_dv2_v0[j*8 + b]*hamil[b]; 
    I[1*4 + j] = 2.0/(dxv[4]*jacob_vel_surf[6])*dH_dvz; 
  } 
  return 2; 
} 

GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[8]; 
  B_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p1_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*4 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
