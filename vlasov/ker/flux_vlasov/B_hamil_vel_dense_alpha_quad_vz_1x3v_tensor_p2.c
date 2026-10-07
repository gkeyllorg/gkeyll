#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_tensor_p2.h> 
GKYL_CU_DH int B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p2_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  const double *By = &qmem[12]; 
  const double *Bx = &qmem[9]; 
  for (int i = tid; i < 4; i += nthreads) { 
    double By_quad = 0.0; 
    for (int a = 0; a < 3; ++a) By_quad += vst_1x3v_tensor_p2_conf_ev[i*3 + a]*By[a]; 
    O[0*4 + i] = By_quad; 
    double Bx_quad = 0.0; 
    for (int a = 0; a < 3; ++a) Bx_quad += vst_1x3v_tensor_p2_conf_ev[i*3 + a]*Bx[a]; 
    O[1*4 + i] = -Bx_quad; 
  } 
  for (int j = tid; j < 16; j += nthreads) { 
    double dH_dvx = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dvx += vst_1x3v_tensor_p2_vel_dv0_v2[j*27 + b]*hamil[b]; 
    I[0*16 + j] = 2.0/(dxv[1]*jacob_vel_surf[0 + j/4])*dH_dvx; 
    double dH_dvy = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dvy += vst_1x3v_tensor_p2_vel_dv1_v2[j*27 + b]*hamil[b]; 
    I[1*16 + j] = 2.0/(dxv[2]*jacob_vel_surf[4 + j%4])*dH_dvy; 
  } 
  return 2; 
} 

GKYL_CU_DH void B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p2(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[32]; 
  B_hamil_vel_dense_alpha_quad_vz_1x3v_tensor_p2_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*4 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
