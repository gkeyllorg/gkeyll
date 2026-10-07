#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_ser_p1.h> 
GKYL_CU_DH int hamil_phase_alpha_quad_vz_3x3v_ser_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 8; 
  const double dx12 = 2.0/dxv[2]; 
  const double jacob_cz_inv = 1.0/jacob_pos[4]; 
  for (int i = tid; i < 8; i += nthreads) { 
    for (int a = 0; a < 8; ++a) O[a*8 + i] = -dx12*jacob_cz_inv*vst_3x3v_ser_p1_ph_v2_CmD[i*8 + a]; 
  } 
  for (int j = tid; j < 4; j += nthreads) { 
    double G[8]; 
    for (int a = 0; a < 8; ++a) G[a] = 0.0; 
    for (int k = 0; k < 64; ++k) { 
      G[vst_3x3v_ser_p1_ph_v2_cmap[k]] += vst_3x3v_ser_p1_ph_v2_V[j*8 + vst_3x3v_ser_p1_ph_v2_vrmap[k]]*(vst_3x3v_ser_p1_ph_v2_coefr[k]*hamil[k]); 
    } 
    for (int a = 0; a < 8; ++a) I[a*4 + j] = G[a]; 
  } 
  return 8; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_vz_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[64]; 
  double I[32]; 
  hamil_phase_alpha_quad_vz_3x3v_ser_p1_shared(0, 1, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 8; ++t) alpha += O[t*8 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
