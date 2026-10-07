#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p2.h> 
GKYL_CU_DH int B_ho_hamil_phase_alpha_quad_vy_2x2v_ser_p2_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 8; 
  const double *Bz = &qmem[40]; 
  for (int i = tid; i < 16; i += nthreads) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 8; ++a) Bz_quad += vst_2x2v_ser_p2_ho_conf_ev[i*8 + a]*Bz[a]; 
    for (int a = 0; a < 8; ++a) O[(0 + a)*16 + i] = -Bz_quad*vst_2x2v_ser_p2_ho_ph_v1_Cm[i*8 + a]; 
  } 
  for (int j = tid; j < 4; j += nthreads) { 
    double Gd0[8]; 
    for (int a = 0; a < 8; ++a) Gd0[a] = 0.0; 
    for (int b = 0; b < 48; ++b) { 
      Gd0[vst_2x2v_ser_p2_ho_ph_v1_cmap[b]] += vst_2x2v_ser_p2_ho_ph_v1_Vd0[j*6 + vst_2x2v_ser_p2_ho_ph_v1_vrd0map[b]]*(vst_2x2v_ser_p2_ho_ph_v1_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 8; ++a) I[(0 + a)*4 + j] = 2.0/(dxv[2]*jacob_vel_surf[0])*Gd0[a]; 
  } 
  return 8; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_2x2v_ser_p2(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[128]; 
  double I[32]; 
  B_ho_hamil_phase_alpha_quad_vy_2x2v_ser_p2_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 8; ++t) alpha += O[t*16 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
