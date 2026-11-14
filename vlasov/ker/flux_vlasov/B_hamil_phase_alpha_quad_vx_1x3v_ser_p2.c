#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH int B_hamil_phase_alpha_quad_vx_1x3v_ser_p2_shared(int tid, int off, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 6; 
  O += off*3; 
  I += off*9; 
  const double *Bz = &qmem[15]; 
  const double *By = &qmem[12]; 
  if (tid < 3) { 
    const int i = tid; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 3; ++a) Bz_quad += vst_1x3v_ser_p2_conf_ev[i*3 + a]*Bz[a]; 
    for (int a = 0; a < 3; ++a) O[(0 + a)*3 + i] = Bz_quad*vst_1x3v_ser_p2_ph_v0_Cm[i*3 + a]; 
    double By_quad = 0.0; 
    for (int a = 0; a < 3; ++a) By_quad += vst_1x3v_ser_p2_conf_ev[i*3 + a]*By[a]; 
    for (int a = 0; a < 3; ++a) O[(3 + a)*3 + i] = -By_quad*vst_1x3v_ser_p2_ph_v0_Cm[i*3 + a]; 
  } 
  if (tid < 9) { 
    const int j = tid; 
    double Gd1[3]; 
    for (int a = 0; a < 3; ++a) Gd1[a] = 0.0; 
    for (int b = 0; b < 48; ++b) { 
      Gd1[vst_1x3v_ser_p2_ph_v0_cmap[b]] += vst_1x3v_ser_p2_ph_v0_Vd1[j*13 + vst_1x3v_ser_p2_ph_v0_vrd1map[b]]*(vst_1x3v_ser_p2_ph_v0_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 3; ++a) I[(0 + a)*9 + j] = 2.0/(dxv[2]*jacob_vel_surf[4])*Gd1[a]; 
    double Gd2[3]; 
    for (int a = 0; a < 3; ++a) Gd2[a] = 0.0; 
    for (int b = 0; b < 48; ++b) { 
      Gd2[vst_1x3v_ser_p2_ph_v0_cmap[b]] += vst_1x3v_ser_p2_ph_v0_Vd2[j*13 + vst_1x3v_ser_p2_ph_v0_vrd2map[b]]*(vst_1x3v_ser_p2_ph_v0_dcoefr2[b]*hamil[b]); 
    } 
    for (int a = 0; a < 3; ++a) I[(3 + a)*9 + j] = 2.0/(dxv[3]*jacob_vel_surf[8])*Gd2[a]; 
  } 
  return 6; 
} 

GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[18]; 
  double I[54]; 
  for (int tid = 0; tid < 9; ++tid) B_hamil_phase_alpha_quad_vx_1x3v_ser_p2_shared(tid, 0, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 6; ++t) alpha += O[t*3 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
