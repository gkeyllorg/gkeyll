#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x1v_ser_p2.h> 
GKYL_CU_DH int hamil_phase_ho_alpha_quad_vx_2x1v_ser_p2_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 8; 
  O += off*16; 
  I += off*1; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  if (tid < 16) { 
    const int i = tid; 
    for (int a = 0; a < 8; ++a) O[a*16 + i] = -dx10*jacob_cx_inv*vst_2x1v_ser_p2_ho_ph_v0_CmD[i*8 + a]; 
  } 
  if (tid < 1) { 
    const int j = tid; 
    double G[8]; 
    for (int a = 0; a < 8; ++a) G[a] = 0.0; 
    for (int k = 0; k < 20; ++k) { 
      G[vst_2x1v_ser_p2_ho_ph_v0_cmap[k]] += vst_2x1v_ser_p2_ho_ph_v0_V[j*3 + vst_2x1v_ser_p2_ho_ph_v0_vrmap[k]]*(vst_2x1v_ser_p2_ho_ph_v0_coefr[k]*hamil[k]); 
    } 
    for (int a = 0; a < 8; ++a) I[a*1 + j] = G[a]; 
  } 
  return 8; 
} 

GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_2x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[128]; 
  double I[8]; 
  for (int tid = 0; tid < 16; ++tid) hamil_phase_ho_alpha_quad_vx_2x1v_ser_p2_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 8; ++t) alpha += O[t*16 + i]*I[t*1 + j]; 
      alpha_quad[i*1 + j] += alpha; 
    } 
  } 
} 
