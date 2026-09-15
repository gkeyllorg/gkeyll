#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x1v_ser_p3.h> 
GKYL_CU_DH double hamil_phase_ho_alpha_quad_vx_1x1v_ser_p3_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  double G[4]; 
  for (int a = 0; a < 4; ++a) G[a] = 0.0; 
  for (int k = 0; k < 12; ++k) { 
    G[vst_1x1v_ser_p3_ho_ph_v0_cmap[k]] += vst_1x1v_ser_p3_ho_ph_v0_V[j*4 + vst_1x1v_ser_p3_ho_ph_v0_vrmap[k]]*(vst_1x1v_ser_p3_ho_ph_v0_coefr[k]*hamil[k]); 
  } 
  double dH_dx = 0.0; 
  for (int a = 0; a < 4; ++a) dH_dx += vst_1x1v_ser_p3_ho_ph_v0_CmD[i*4 + a]*G[a]; 
  return -dx10*jacob_cx_inv*dH_dx; 
} 

GKYL_CU_DH void hamil_phase_ho_alpha_quad_vx_1x1v_ser_p3(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 5; ++i) { 
    for (int j = 0; j < 1; ++j) alpha_quad[i*1 + j] += hamil_phase_ho_alpha_quad_vx_1x1v_ser_p3_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
