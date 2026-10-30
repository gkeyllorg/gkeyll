#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_tensor_p1.h> 
GKYL_CU_DH int B_ho_hamil_phase_alpha_quad_vz_1x3v_tensor_p1_shared(int tid, int off, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 4; 
  O += off*2; 
  I += off*16; 
  const double *By = &qmem[8]; 
  const double *Bx = &qmem[6]; 
  if (tid < 2) { 
    const int i = tid; 
    double By_quad = 0.0; 
    for (int a = 0; a < 2; ++a) By_quad += vst_1x3v_tensor_p1_ho_conf_ev[i*2 + a]*By[a]; 
    for (int a = 0; a < 2; ++a) O[(0 + a)*2 + i] = By_quad*vst_1x3v_tensor_p1_ho_ph_v2_Cm[i*2 + a]; 
    double Bx_quad = 0.0; 
    for (int a = 0; a < 2; ++a) Bx_quad += vst_1x3v_tensor_p1_ho_conf_ev[i*2 + a]*Bx[a]; 
    for (int a = 0; a < 2; ++a) O[(2 + a)*2 + i] = -Bx_quad*vst_1x3v_tensor_p1_ho_ph_v2_Cm[i*2 + a]; 
  } 
  if (tid < 16) { 
    const int j = tid; 
    double Gd0[2]; 
    for (int a = 0; a < 2; ++a) Gd0[a] = 0.0; 
    for (int b = 0; b < 54; ++b) { 
      Gd0[vst_1x3v_tensor_p1_ho_ph_v2_cmap[b]] += vst_1x3v_tensor_p1_ho_ph_v2_Vd0[j*19 + vst_1x3v_tensor_p1_ho_ph_v2_vrd0map[b]]*(vst_1x3v_tensor_p1_ho_ph_v2_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 2; ++a) I[(0 + a)*16 + j] = 2.0/(dxv[1]*jacob_vel_surf[0 + j/4])*Gd0[a]; 
    double Gd1[2]; 
    for (int a = 0; a < 2; ++a) Gd1[a] = 0.0; 
    for (int b = 0; b < 54; ++b) { 
      Gd1[vst_1x3v_tensor_p1_ho_ph_v2_cmap[b]] += vst_1x3v_tensor_p1_ho_ph_v2_Vd1[j*19 + vst_1x3v_tensor_p1_ho_ph_v2_vrd1map[b]]*(vst_1x3v_tensor_p1_ho_ph_v2_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 2; ++a) I[(2 + a)*16 + j] = 2.0/(dxv[2]*jacob_vel_surf[4 + j%4])*Gd1[a]; 
  } 
  return 4; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vz_1x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[64]; 
  for (int tid = 0; tid < 16; ++tid) B_ho_hamil_phase_alpha_quad_vz_1x3v_tensor_p1_shared(tid, 0, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 4; ++t) alpha += O[t*2 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
