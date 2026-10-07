#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_tensor_p1.h> 
GKYL_CU_DH int hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  for (int i = tid; i < 2; i += nthreads) { 
    for (int a = 0; a < 2; ++a) O[a*2 + i] = -dx10*jacob_cx_inv*vst_1x2v_tensor_p1_ho_ph_v0_CmD[i*2 + a]; 
  } 
  for (int j = tid; j < 4; j += nthreads) { 
    double G[2]; 
    for (int a = 0; a < 2; ++a) G[a] = 0.0; 
    for (int k = 0; k < 18; ++k) { 
      G[vst_1x2v_tensor_p1_ho_ph_v0_cmap[k]] += vst_1x2v_tensor_p1_ho_ph_v0_V[j*9 + vst_1x2v_tensor_p1_ho_ph_v0_vrmap[k]]*(vst_1x2v_tensor_p1_ho_ph_v0_coefr[k]*hamil[k]); 
    } 
    for (int a = 0; a < 2; ++a) I[a*4 + j] = G[a]; 
  } 
  return 2; 
} 

GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[4]; 
  double I[8]; 
  hamil_phase_ho_alpha_quad_vx_1x2v_tensor_p1_shared(0, 1, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*2 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
