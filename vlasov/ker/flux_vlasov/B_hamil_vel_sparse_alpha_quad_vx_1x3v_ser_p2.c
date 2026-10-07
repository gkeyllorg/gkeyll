#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH int B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  const double *Bz = &qmem[15]; 
  const double *By = &qmem[12]; 
  for (int i = tid; i < 3; i += nthreads) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 3; ++a) Bz_quad += vst_1x3v_ser_p2_conf_ev[i*3 + a]*Bz[a]; 
    O[0*3 + i] = Bz_quad; 
    double By_quad = 0.0; 
    for (int a = 0; a < 3; ++a) By_quad += vst_1x3v_ser_p2_conf_ev[i*3 + a]*By[a]; 
    O[1*3 + i] = -By_quad; 
  } 
  for (int j = tid; j < 9; j += nthreads) { 
    double dH_dvy = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_1x3v_ser_p2_vel_sparse_idx[s]; 
      dH_dvy += vst_1x3v_ser_p2_vel_dv1_v0[j*20 + b]*hamil[b]; 
    } 
    I[0*9 + j] = 2.0/(dxv[2]*jacob_vel_surf[4])*dH_dvy; 
    double dH_dvz = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_1x3v_ser_p2_vel_sparse_idx[s]; 
      dH_dvz += vst_1x3v_ser_p2_vel_dv2_v0[j*20 + b]*hamil[b]; 
    } 
    I[1*9 + j] = 2.0/(dxv[3]*jacob_vel_surf[8])*dH_dvz; 
  } 
  return 2; 
} 

GKYL_CU_DH void B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[6]; 
  double I[18]; 
  B_hamil_vel_sparse_alpha_quad_vx_1x3v_ser_p2_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*3 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
