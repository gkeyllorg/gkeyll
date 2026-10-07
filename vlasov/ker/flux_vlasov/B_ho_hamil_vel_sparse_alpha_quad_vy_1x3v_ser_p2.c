#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH int B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  const double *Bx = &qmem[9]; 
  const double *Bz = &qmem[15]; 
  for (int i = tid; i < 4; i += nthreads) { 
    double Bx_quad = 0.0; 
    for (int a = 0; a < 3; ++a) Bx_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*Bx[a]; 
    O[0*4 + i] = Bx_quad; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 3; ++a) Bz_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*Bz[a]; 
    O[1*4 + i] = -Bz_quad; 
  } 
  for (int j = tid; j < 16; j += nthreads) { 
    double dH_dvz = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_1x3v_ser_p2_ho_vel_sparse_idx[s]; 
      dH_dvz += vst_1x3v_ser_p2_ho_vel_dv2_v1[j*20 + b]*hamil[b]; 
    } 
    I[0*16 + j] = 2.0/(dxv[3]*jacob_vel_surf[8])*dH_dvz; 
    double dH_dvx = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_1x3v_ser_p2_ho_vel_sparse_idx[s]; 
      dH_dvx += vst_1x3v_ser_p2_ho_vel_dv0_v1[j*20 + b]*hamil[b]; 
    } 
    I[1*16 + j] = 2.0/(dxv[1]*jacob_vel_surf[0])*dH_dvx; 
  } 
  return 2; 
} 

GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[32]; 
  B_ho_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p2_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*4 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
