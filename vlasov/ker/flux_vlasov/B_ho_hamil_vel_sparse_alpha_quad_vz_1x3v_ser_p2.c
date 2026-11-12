#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH int B_ho_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2_shared(int tid, int off, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  O += off*4; 
  I += off*16; 
  const double *By = &qmem[12]; 
  const double *Bx = &qmem[9]; 
  if (tid < 4) { 
    const int i = tid; 
    double By_quad = 0.0; 
    for (int a = 0; a < 3; ++a) By_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*By[a]; 
    O[0*4 + i] = By_quad; 
    double Bx_quad = 0.0; 
    for (int a = 0; a < 3; ++a) Bx_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*Bx[a]; 
    O[1*4 + i] = -Bx_quad; 
  } 
  if (tid < 16) { 
    const int j = tid; 
    double dH_dvx = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_1x3v_ser_p2_ho_vel_sparse_idx[s]; 
      dH_dvx += vst_1x3v_ser_p2_ho_vel_dv0_v2[j*20 + b]*hamil[b]; 
    } 
    I[0*16 + j] = 2.0/(dxv[1]*jacob_vel_surf[0])*dH_dvx; 
    double dH_dvy = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_1x3v_ser_p2_ho_vel_sparse_idx[s]; 
      dH_dvy += vst_1x3v_ser_p2_ho_vel_dv1_v2[j*20 + b]*hamil[b]; 
    } 
    I[1*16 + j] = 2.0/(dxv[2]*jacob_vel_surf[4])*dH_dvy; 
  } 
  return 2; 
} 

GKYL_CU_DH void B_ho_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[32]; 
  for (int tid = 0; tid < 16; ++tid) B_ho_hamil_vel_sparse_alpha_quad_vz_1x3v_ser_p2_shared(tid, 0, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*4 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
