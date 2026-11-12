#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_ser_p2.h> 
GKYL_CU_DH int B_ho_hamil_phase_alpha_quad_vx_1x2v_ser_p2_shared(int tid, int off, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 3; 
  O += off*4; 
  I += off*4; 
  const double *Bz = &qmem[15]; 
  if (tid < 4) { 
    const int i = tid; 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 3; ++a) Bz_quad += vst_1x2v_ser_p2_ho_conf_ev[i*3 + a]*Bz[a]; 
    for (int a = 0; a < 3; ++a) O[(0 + a)*4 + i] = Bz_quad*vst_1x2v_ser_p2_ho_ph_v0_Cm[i*3 + a]; 
  } 
  if (tid < 4) { 
    const int j = tid; 
    double Gd1[3]; 
    for (int a = 0; a < 3; ++a) Gd1[a] = 0.0; 
    for (int b = 0; b < 20; ++b) { 
      Gd1[vst_1x2v_ser_p2_ho_ph_v0_cmap[b]] += vst_1x2v_ser_p2_ho_ph_v0_Vd1[j*6 + vst_1x2v_ser_p2_ho_ph_v0_vrd1map[b]]*(vst_1x2v_ser_p2_ho_ph_v0_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 3; ++a) I[(0 + a)*4 + j] = 2.0/(dxv[2]*jacob_vel_surf[4])*Gd1[a]; 
  } 
  return 3; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_1x2v_ser_p2(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[12]; 
  double I[12]; 
  for (int tid = 0; tid < 4; ++tid) B_ho_hamil_phase_alpha_quad_vx_1x2v_ser_p2_shared(tid, 0, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 3; ++t) alpha += O[t*4 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
