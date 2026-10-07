#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH int B_ho_hamil_phase_alpha_quad_vz_1x3v_ser_p2_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 6; 
  const double *By = &qmem[12]; 
  const double *Bx = &qmem[9]; 
  for (int i = tid; i < 4; i += nthreads) { 
    double By_quad = 0.0; 
    for (int a = 0; a < 3; ++a) By_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*By[a]; 
    for (int a = 0; a < 3; ++a) O[(0 + a)*4 + i] = By_quad*vst_1x3v_ser_p2_ho_ph_v2_Cm[i*3 + a]; 
    double Bx_quad = 0.0; 
    for (int a = 0; a < 3; ++a) Bx_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*Bx[a]; 
    for (int a = 0; a < 3; ++a) O[(3 + a)*4 + i] = -Bx_quad*vst_1x3v_ser_p2_ho_ph_v2_Cm[i*3 + a]; 
  } 
  for (int j = tid; j < 16; j += nthreads) { 
    double Gd0[3]; 
    for (int a = 0; a < 3; ++a) Gd0[a] = 0.0; 
    for (int b = 0; b < 48; ++b) { 
      Gd0[vst_1x3v_ser_p2_ho_ph_v2_cmap[b]] += vst_1x3v_ser_p2_ho_ph_v2_Vd0[j*13 + vst_1x3v_ser_p2_ho_ph_v2_vrd0map[b]]*(vst_1x3v_ser_p2_ho_ph_v2_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 3; ++a) I[(0 + a)*16 + j] = 2.0/(dxv[1]*jacob_vel_surf[0])*Gd0[a]; 
    double Gd1[3]; 
    for (int a = 0; a < 3; ++a) Gd1[a] = 0.0; 
    for (int b = 0; b < 48; ++b) { 
      Gd1[vst_1x3v_ser_p2_ho_ph_v2_cmap[b]] += vst_1x3v_ser_p2_ho_ph_v2_Vd1[j*13 + vst_1x3v_ser_p2_ho_ph_v2_vrd1map[b]]*(vst_1x3v_ser_p2_ho_ph_v2_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 3; ++a) I[(3 + a)*16 + j] = 2.0/(dxv[2]*jacob_vel_surf[4])*Gd1[a]; 
  } 
  return 6; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[24]; 
  double I[96]; 
  B_ho_hamil_phase_alpha_quad_vz_1x3v_ser_p2_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 6; ++t) alpha += O[t*4 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
