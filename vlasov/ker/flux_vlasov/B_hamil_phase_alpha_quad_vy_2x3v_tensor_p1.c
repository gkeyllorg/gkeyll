#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH int B_hamil_phase_alpha_quad_vy_2x3v_tensor_p1_shared(int tid, int off, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 8; 
  O += off*4; 
  I += off*9; 
  const double *Bx = &qmem[12]; 
  const double *Bz = &qmem[20]; 
  if (tid < 4) { 
    const int i = tid; 
    double Bx_quad = 0.0; 
    for (int a = 0; a < 4; ++a) Bx_quad += vst_2x3v_tensor_p1_conf_ev[i*4 + a]*Bx[a]; 
    for (int a = 0; a < 4; ++a) O[(0 + a)*4 + i] = Bx_quad*vst_2x3v_tensor_p1_ph_v1_Cm[i*4 + a]; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 4; ++a) Bz_quad += vst_2x3v_tensor_p1_conf_ev[i*4 + a]*Bz[a]; 
    for (int a = 0; a < 4; ++a) O[(4 + a)*4 + i] = -Bz_quad*vst_2x3v_tensor_p1_ph_v1_Cm[i*4 + a]; 
  } 
  if (tid < 9) { 
    const int j = tid; 
    double Gd2[4]; 
    for (int a = 0; a < 4; ++a) Gd2[a] = 0.0; 
    for (int b = 0; b < 108; ++b) { 
      Gd2[vst_2x3v_tensor_p1_ph_v1_cmap[b]] += vst_2x3v_tensor_p1_ph_v1_Vd2[j*19 + vst_2x3v_tensor_p1_ph_v1_vrd2map[b]]*(vst_2x3v_tensor_p1_ph_v1_dcoefr2[b]*hamil[b]); 
    } 
    for (int a = 0; a < 4; ++a) I[(0 + a)*9 + j] = 2.0/(dxv[4]*jacob_vel_surf[6 + j%3])*Gd2[a]; 
    double Gd0[4]; 
    for (int a = 0; a < 4; ++a) Gd0[a] = 0.0; 
    for (int b = 0; b < 108; ++b) { 
      Gd0[vst_2x3v_tensor_p1_ph_v1_cmap[b]] += vst_2x3v_tensor_p1_ph_v1_Vd0[j*19 + vst_2x3v_tensor_p1_ph_v1_vrd0map[b]]*(vst_2x3v_tensor_p1_ph_v1_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 4; ++a) I[(4 + a)*9 + j] = 2.0/(dxv[2]*jacob_vel_surf[0 + j/3])*Gd0[a]; 
  } 
  return 8; 
} 

GKYL_CU_DH void B_hamil_phase_alpha_quad_vy_2x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[32]; 
  double I[72]; 
  for (int tid = 0; tid < 9; ++tid) B_hamil_phase_alpha_quad_vy_2x3v_tensor_p1_shared(tid, 0, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 8; ++t) alpha += O[t*4 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
