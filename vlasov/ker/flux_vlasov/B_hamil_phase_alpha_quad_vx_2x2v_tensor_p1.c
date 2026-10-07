#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH int B_hamil_phase_alpha_quad_vx_2x2v_tensor_p1_shared(int tid, int nthreads, const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 4; 
  const double *Bz = &qmem[20]; 
  for (int i = tid; i < 4; i += nthreads) { 
    double Bz_quad = 0.0; 
    for (int a = 0; a < 4; ++a) Bz_quad += vst_2x2v_tensor_p1_conf_ev[i*4 + a]*Bz[a]; 
    for (int a = 0; a < 4; ++a) O[(0 + a)*4 + i] = Bz_quad*vst_2x2v_tensor_p1_ph_v0_Cm[i*4 + a]; 
  } 
  for (int j = tid; j < 3; j += nthreads) { 
    double Gd1[4]; 
    for (int a = 0; a < 4; ++a) Gd1[a] = 0.0; 
    for (int b = 0; b < 36; ++b) { 
      Gd1[vst_2x2v_tensor_p1_ph_v0_cmap[b]] += vst_2x2v_tensor_p1_ph_v0_Vd1[j*7 + vst_2x2v_tensor_p1_ph_v0_vrd1map[b]]*(vst_2x2v_tensor_p1_ph_v0_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 4; ++a) I[(0 + a)*3 + j] = 2.0/(dxv[3]*jacob_vel_surf[3 + j])*Gd1[a]; 
  } 
  return 4; 
} 

GKYL_CU_DH void B_hamil_phase_alpha_quad_vx_2x2v_tensor_p1(const double *dxv, const double *jacob_vel_surf, const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[16]; 
  double I[12]; 
  B_hamil_phase_alpha_quad_vx_2x2v_tensor_p1_shared(0, 1, dxv, jacob_vel_surf, hamil, qmem, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 4; ++t) alpha += O[t*4 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
