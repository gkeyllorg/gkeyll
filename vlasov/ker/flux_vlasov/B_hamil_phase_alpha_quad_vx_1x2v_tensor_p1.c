#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_tensor_p1.h> 
GKYL_CU_DH int B_hamil_phase_alpha_quad_vx_1x2v_tensor_p1_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  const double *Bz = &qmem[10]; 
  for (int i = tid; i < 2; i += nthreads) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 2; ++a) Bz_quad += vst_1x2v_tensor_p1_conf_ev[i*2 + a]*Bz[a]; 
    for (int a = 0; a < 2; ++a) O[(0 + a)*2 + i] = Bz_quad*vst_1x2v_tensor_p1_ph_v0_Cm[i*2 + a]; 
  } 
  for (int j = tid; j < 3; j += nthreads) { 
    double Gd1[2]; 
    for (int a = 0; a < 2; ++a) Gd1[a] = 0.0; 
    for (int b = 0; b < 18; ++b) { 
      Gd1[vst_1x2v_tensor_p1_ph_v0_cmap[b]] += vst_1x2v_tensor_p1_ph_v0_Vd1[j*7 + vst_1x2v_tensor_p1_ph_v0_vrd1map[b]]*(vst_1x2v_tensor_p1_ph_v0_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 2; ++a) I[(0 + a)*3 + j] = 2.0/(dxv[2]*jacob_vel_surf[3 + j])*Gd1[a]; 
  } 
  return 2; 
} 

GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[4]; 
  double I[6]; 
  B_hamil_phase_alpha_quad_vx_1x2v_tensor_p1_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*2 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
