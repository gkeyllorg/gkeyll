#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_ser_p1.h> 
GKYL_CU_DH double hamil_phase_alpha_quad_vz_3x3v_ser_p1_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dx12 = 2.0/dxv[2]; 
  const double *jacob_cz = &jacob_pos[4]; 
  const double jacob_cz_inv = 1.0/jacob_cz[0]; 
  double G[8]; 
  for (int a = 0; a < 8; ++a) G[a] = 0.0; 
  for (int k = 0; k < 64; ++k) { 
    G[vst_3x3v_ser_p1_ph_v2_cmap[k]] += vst_3x3v_ser_p1_ph_v2_V[j*8 + vst_3x3v_ser_p1_ph_v2_vrmap[k]]*(vst_3x3v_ser_p1_ph_v2_coefr[k]*hamil[k]); 
  } 
  double dH_dx = 0.0; 
  for (int a = 0; a < 8; ++a) dH_dx += vst_3x3v_ser_p1_ph_v2_CmD[i*8 + a]*G[a]; 
  return -dx12*jacob_cz_inv*dH_dx; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_vz_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += hamil_phase_alpha_quad_vz_3x3v_ser_p1_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
