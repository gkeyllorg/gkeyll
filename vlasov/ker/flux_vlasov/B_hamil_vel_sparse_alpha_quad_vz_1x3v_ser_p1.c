#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p1.h> 
GKYL_CU_DH int B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  const double *By = &qmem[8]; 
  const double *Bx = &qmem[6]; 
  for (int i = tid; i < 2; i += nthreads) { 
    double By_quad = 0.0; 
    for (int a = 0; a < 2; ++a) By_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*By[a]; 
    O[0*2 + i] = By_quad; 
    double Bx_quad = 0.0; 
    for (int a = 0; a < 2; ++a) Bx_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*Bx[a]; 
    O[1*2 + i] = -Bx_quad; 
  } 
  for (int j = tid; j < 4; j += nthreads) { 
    double dH_dvx = 0.0; 
    for (int s = 0; s < 4; ++s) { 
      const int b = vst_1x3v_ser_p1_vel_sparse_idx[s]; 
      dH_dvx += vst_1x3v_ser_p1_vel_dv0_v2[j*8 + b]*hamil[b]; 
    } 
    I[0*4 + j] = 2.0/(dxv[1]*jacob_vel_surf[0])*dH_dvx; 
    double dH_dvy = 0.0; 
    for (int s = 0; s < 4; ++s) { 
      const int b = vst_1x3v_ser_p1_vel_sparse_idx[s]; 
      dH_dvy += vst_1x3v_ser_p1_vel_dv1_v2[j*8 + b]*hamil[b]; 
    } 
    I[1*4 + j] = 2.0/(dxv[2]*jacob_vel_surf[3])*dH_dvy; 
  } 
  return 2; 
} 

GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[4]; 
  double I[8]; 
  B_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p1_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*2 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
