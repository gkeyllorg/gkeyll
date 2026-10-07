#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_tensor_p1.h> 
GKYL_CU_DH int B_ho_hamil_phase_alpha_quad_vy_1x2v_tensor_p1_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  const double *Bz = &qmem[10]; 
  for (int i = tid; i < 2; i += nthreads) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 2; ++a) Bz_quad += vst_1x2v_tensor_p1_ho_conf_ev[i*2 + a]*Bz[a]; 
    for (int a = 0; a < 2; ++a) O[(0 + a)*2 + i] = -Bz_quad*vst_1x2v_tensor_p1_ho_ph_v1_Cm[i*2 + a]; 
  } 
  for (int j = tid; j < 4; j += nthreads) { 
    double Gd0[2]; 
    for (int a = 0; a < 2; ++a) Gd0[a] = 0.0; 
    for (int b = 0; b < 18; ++b) { 
      Gd0[vst_1x2v_tensor_p1_ho_ph_v1_cmap[b]] += vst_1x2v_tensor_p1_ho_ph_v1_Vd0[j*7 + vst_1x2v_tensor_p1_ho_ph_v1_vrd0map[b]]*(vst_1x2v_tensor_p1_ho_ph_v1_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 2; ++a) I[(0 + a)*4 + j] = 2.0/(dxv[1]*jacob_vel_surf[0 + j])*Gd0[a]; 
  } 
  return 2; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_1x2v_tensor_p1(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[4]; 
  double I[8]; 
  B_ho_hamil_phase_alpha_quad_vy_1x2v_tensor_p1_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*2 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
