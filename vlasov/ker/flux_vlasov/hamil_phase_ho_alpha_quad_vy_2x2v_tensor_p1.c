#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH double hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dx11 = 2.0/dxv[1]; 
  const double *jacob_cy = &jacob_pos[2]; 
  const double jacob_cy_inv = 1.0/jacob_cy[0]; 
  double G[4]; 
  for (int a = 0; a < 4; ++a) G[a] = 0.0; 
  for (int k = 0; k < 36; ++k) { 
    G[vst_2x2v_tensor_p1_ho_ph_v1_cmap[k]] += vst_2x2v_tensor_p1_ho_ph_v1_V[j*9 + vst_2x2v_tensor_p1_ho_ph_v1_vrmap[k]]*(vst_2x2v_tensor_p1_ho_ph_v1_coefr[k]*hamil[k]); 
  } 
  double dH_dx = 0.0; 
  for (int a = 0; a < 4; ++a) dH_dx += vst_2x2v_tensor_p1_ho_ph_v1_CmD[i*4 + a]*G[a]; 
  return -dx11*jacob_cy_inv*dH_dx; 
} 

GKYL_CU_DH void hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += hamil_phase_ho_alpha_quad_vy_2x2v_tensor_p1_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
