#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_ser_p2.h> 
GKYL_CU_DH int B_ho_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2_shared(int tid, int off, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  O += off*16; 
  I += off*16; 
  const double *Bz = &qmem[40]; 
  const double *By = &qmem[32]; 
  if (tid < 16) { 
    const int i = tid; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 8; ++a) Bz_quad += vst_2x3v_ser_p2_ho_conf_ev[i*8 + a]*Bz[a]; 
    O[0*16 + i] = Bz_quad; 
    double By_quad = 0.0; 
    for (int a = 0; a < 8; ++a) By_quad += vst_2x3v_ser_p2_ho_conf_ev[i*8 + a]*By[a]; 
    O[1*16 + i] = -By_quad; 
  } 
  if (tid < 16) { 
    const int j = tid; 
    double dH_dvy = 0.0; 
    for (int b = 0; b < 20; ++b) dH_dvy += vst_2x3v_ser_p2_ho_vel_dv1_v0[j*20 + b]*hamil[b]; 
    I[0*16 + j] = 2.0/(dxv[3]*jacob_vel_surf[4])*dH_dvy; 
    double dH_dvz = 0.0; 
    for (int b = 0; b < 20; ++b) dH_dvz += vst_2x3v_ser_p2_ho_vel_dv2_v0[j*20 + b]*hamil[b]; 
    I[1*16 + j] = 2.0/(dxv[4]*jacob_vel_surf[8])*dH_dvz; 
  } 
  return 2; 
} 

GKYL_CU_DH void B_ho_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[32]; 
  double I[32]; 
  for (int tid = 0; tid < 16; ++tid) B_ho_hamil_vel_dense_alpha_quad_vx_2x3v_ser_p2_shared(tid, 0, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*16 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
