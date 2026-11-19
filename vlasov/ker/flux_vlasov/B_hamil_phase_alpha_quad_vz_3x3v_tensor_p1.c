#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH int B_hamil_phase_alpha_quad_vz_3x3v_tensor_p1_shared(int tid, int off, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 16; 
  O += off*8; 
  I += off*9; 
  const double *By = &qmem[32]; 
  const double *Bx = &qmem[24]; 
  if (tid < 8) { 
    const int i = tid; 
    double By_quad = 0.0; 
    for (int a = 0; a < 8; ++a) By_quad += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*By[a]; 
    for (int a = 0; a < 8; ++a) O[(0 + a)*8 + i] = By_quad*vst_3x3v_tensor_p1_ph_v2_Cm[i*8 + a]; 
    double Bx_quad = 0.0; 
    for (int a = 0; a < 8; ++a) Bx_quad += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*Bx[a]; 
    for (int a = 0; a < 8; ++a) O[(8 + a)*8 + i] = -Bx_quad*vst_3x3v_tensor_p1_ph_v2_Cm[i*8 + a]; 
  } 
  if (tid < 9) { 
    const int j = tid; 
    double Gd0[8]; 
    for (int a = 0; a < 8; ++a) Gd0[a] = 0.0; 
    for (int b = 0; b < 216; ++b) { 
      Gd0[vst_3x3v_tensor_p1_ph_v2_cmap[b]] += vst_3x3v_tensor_p1_ph_v2_Vd0[j*19 + vst_3x3v_tensor_p1_ph_v2_vrd0map[b]]*(vst_3x3v_tensor_p1_ph_v2_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 8; ++a) I[(0 + a)*9 + j] = 2.0/(dxv[3]*jacob_vel_surf[0 + j/3])*Gd0[a]; 
    double Gd1[8]; 
    for (int a = 0; a < 8; ++a) Gd1[a] = 0.0; 
    for (int b = 0; b < 216; ++b) { 
      Gd1[vst_3x3v_tensor_p1_ph_v2_cmap[b]] += vst_3x3v_tensor_p1_ph_v2_Vd1[j*19 + vst_3x3v_tensor_p1_ph_v2_vrd1map[b]]*(vst_3x3v_tensor_p1_ph_v2_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 8; ++a) I[(8 + a)*9 + j] = 2.0/(dxv[4]*jacob_vel_surf[3 + j%3])*Gd1[a]; 
  } 
  return 16; 
} 

GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_3x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[128]; 
  double I[144]; 
  for (int tid = 0; tid < 9; ++tid) B_hamil_phase_alpha_quad_vz_3x3v_tensor_p1_shared(tid, 0, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 16; ++t) alpha += O[t*8 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
