#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_ser_p2.h> 
GKYL_CU_DH int hamil_phase_alpha_quad_vx_1x2v_ser_p2_shared(int tid, int nthreads, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 3; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  for (int i = tid; i < 3; i += nthreads) { 
    for (int a = 0; a < 3; ++a) O[a*3 + i] = -dx10*jacob_cx_inv*vst_1x2v_ser_p2_ph_v0_CmD[i*3 + a]; 
  } 
  for (int j = tid; j < 3; j += nthreads) { 
    double G[3]; 
    for (int a = 0; a < 3; ++a) G[a] = 0.0; 
    for (int k = 0; k < 20; ++k) { 
      G[vst_1x2v_ser_p2_ph_v0_cmap[k]] += vst_1x2v_ser_p2_ph_v0_V[j*8 + vst_1x2v_ser_p2_ph_v0_vrmap[k]]*(vst_1x2v_ser_p2_ph_v0_coefr[k]*hamil[k]); 
    } 
    for (int a = 0; a < 3; ++a) I[a*3 + j] = G[a]; 
  } 
  return 3; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[9]; 
  double I[9]; 
  hamil_phase_alpha_quad_vx_1x2v_ser_p2_shared(0, 1, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 3; ++t) alpha += O[t*3 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
