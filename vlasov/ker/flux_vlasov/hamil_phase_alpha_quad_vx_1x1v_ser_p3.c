#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x1v_ser_p3.h> 
GKYL_CU_DH int hamil_phase_alpha_quad_vx_1x1v_ser_p3_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 4; 
  O += off*4; 
  I += off*1; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  if (tid < 4) { 
    const int i = tid; 
    for (int a = 0; a < 4; ++a) O[a*4 + i] = -dx10*jacob_cx_inv*vst_1x1v_ser_p3_ph_v0_CmD[i*4 + a]; 
  } 
  if (tid < 1) { 
    const int j = tid; 
    double G[4]; 
    for (int a = 0; a < 4; ++a) G[a] = 0.0; 
    for (int k = 0; k < 12; ++k) { 
      G[vst_1x1v_ser_p3_ph_v0_cmap[k]] += vst_1x1v_ser_p3_ph_v0_V[j*4 + vst_1x1v_ser_p3_ph_v0_vrmap[k]]*(vst_1x1v_ser_p3_ph_v0_coefr[k]*hamil[k]); 
    } 
    for (int a = 0; a < 4; ++a) I[a*1 + j] = G[a]; 
  } 
  return 4; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_vx_1x1v_ser_p3(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[16]; 
  double I[4]; 
  for (int tid = 0; tid < 4; ++tid) hamil_phase_alpha_quad_vx_1x1v_ser_p3_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 4; ++t) alpha += O[t*4 + i]*I[t*1 + j]; 
      alpha_quad[i*1 + j] += alpha; 
    } 
  } 
} 
