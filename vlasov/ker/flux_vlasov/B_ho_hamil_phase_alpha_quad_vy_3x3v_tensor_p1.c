#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH int B_ho_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 16; 
  const double *Bx = &qmem[24]; 
  const double *Bz = &qmem[40]; 
  for (int i = tid; i < 8; i += nthreads) { 
    double Bx_quad = 0.0; 
    for (int a = 0; a < 8; ++a) Bx_quad += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*Bx[a]; 
    for (int a = 0; a < 8; ++a) O[(0 + a)*8 + i] = Bx_quad*vst_3x3v_tensor_p1_ho_ph_v1_Cm[i*8 + a]; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 8; ++a) Bz_quad += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*Bz[a]; 
    for (int a = 0; a < 8; ++a) O[(8 + a)*8 + i] = -Bz_quad*vst_3x3v_tensor_p1_ho_ph_v1_Cm[i*8 + a]; 
  } 
  for (int j = tid; j < 16; j += nthreads) { 
    double Gd2[8]; 
    for (int a = 0; a < 8; ++a) Gd2[a] = 0.0; 
    for (int b = 0; b < 216; ++b) { 
      Gd2[vst_3x3v_tensor_p1_ho_ph_v1_cmap[b]] += vst_3x3v_tensor_p1_ho_ph_v1_Vd2[j*19 + vst_3x3v_tensor_p1_ho_ph_v1_vrd2map[b]]*(vst_3x3v_tensor_p1_ho_ph_v1_dcoefr2[b]*hamil[b]); 
    } 
    for (int a = 0; a < 8; ++a) I[(0 + a)*16 + j] = 2.0/(dxv[5]*jacob_vel_surf[8 + j%4])*Gd2[a]; 
    double Gd0[8]; 
    for (int a = 0; a < 8; ++a) Gd0[a] = 0.0; 
    for (int b = 0; b < 216; ++b) { 
      Gd0[vst_3x3v_tensor_p1_ho_ph_v1_cmap[b]] += vst_3x3v_tensor_p1_ho_ph_v1_Vd0[j*19 + vst_3x3v_tensor_p1_ho_ph_v1_vrd0map[b]]*(vst_3x3v_tensor_p1_ho_ph_v1_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 8; ++a) I[(8 + a)*16 + j] = 2.0/(dxv[3]*jacob_vel_surf[0 + j/4])*Gd0[a]; 
  } 
  return 16; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[128]; 
  double I[256]; 
  B_ho_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 16; ++t) alpha += O[t*8 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
