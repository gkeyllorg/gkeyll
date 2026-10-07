#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p1.h> 
GKYL_CU_DH int hamil_phase_alpha_quad_vy_2x2v_ser_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 4; 
  const double dx11 = 2.0/dxv[1]; 
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  for (int i = tid; i < 4; i += nthreads) { 
    for (int a = 0; a < 4; ++a) O[a*4 + i] = -dx11*jacob_cy_inv*vst_2x2v_ser_p1_ph_v1_CmD[i*4 + a]; 
  } 
  for (int j = tid; j < 2; j += nthreads) { 
    double G[4]; 
    for (int a = 0; a < 4; ++a) G[a] = 0.0; 
    for (int k = 0; k < 16; ++k) { 
      G[vst_2x2v_ser_p1_ph_v1_cmap[k]] += vst_2x2v_ser_p1_ph_v1_V[j*4 + vst_2x2v_ser_p1_ph_v1_vrmap[k]]*(vst_2x2v_ser_p1_ph_v1_coefr[k]*hamil[k]); 
    } 
    for (int a = 0; a < 4; ++a) I[a*2 + j] = G[a]; 
  } 
  return 4; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_vy_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[16]; 
  double I[8]; 
  hamil_phase_alpha_quad_vy_2x2v_ser_p1_shared(0, 1, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 2; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 4; ++t) alpha += O[t*4 + i]*I[t*2 + j]; 
      alpha_quad[i*2 + j] += alpha; 
    } 
  } 
} 
