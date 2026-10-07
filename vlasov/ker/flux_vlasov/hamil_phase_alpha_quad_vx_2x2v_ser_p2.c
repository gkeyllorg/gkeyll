#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p2.h> 
GKYL_CU_DH int hamil_phase_alpha_quad_vx_2x2v_ser_p2_shared(int tid, int nthreads, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 8; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  for (int i = tid; i < 9; i += nthreads) { 
    for (int a = 0; a < 8; ++a) O[a*9 + i] = -dx10*jacob_cx_inv*vst_2x2v_ser_p2_ph_v0_CmD[i*8 + a]; 
  } 
  for (int j = tid; j < 3; j += nthreads) { 
    double G[8]; 
    for (int a = 0; a < 8; ++a) G[a] = 0.0; 
    for (int k = 0; k < 48; ++k) { 
      G[vst_2x2v_ser_p2_ph_v0_cmap[k]] += vst_2x2v_ser_p2_ph_v0_V[j*8 + vst_2x2v_ser_p2_ph_v0_vrmap[k]]*(vst_2x2v_ser_p2_ph_v0_coefr[k]*hamil[k]); 
    } 
    for (int a = 0; a < 8; ++a) I[a*3 + j] = G[a]; 
  } 
  return 8; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_vx_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[72]; 
  double I[24]; 
  hamil_phase_alpha_quad_vx_2x2v_ser_p2_shared(0, 1, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 9; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 8; ++t) alpha += O[t*9 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
