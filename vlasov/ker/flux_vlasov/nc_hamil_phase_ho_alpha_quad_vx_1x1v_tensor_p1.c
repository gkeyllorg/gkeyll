#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x1v_tensor_p1.h> 
GKYL_CU_DH double nc_hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double G[2]; 
  for (int a = 0; a < 2; ++a) { G[a] = 0.0; } 
  for (int k = 0; k < 6; ++k) { 
    const int a = vst_1x1v_tensor_p1_ho_ph_v0_cmap[k]; 
    G[a] += vst_1x1v_tensor_p1_ho_ph_v0_V[j*3 + vst_1x1v_tensor_p1_ho_ph_v0_vrmap[k]]*(vst_1x1v_tensor_p1_ho_ph_v0_coefr[k]*hamil[k]); 
  } 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double *poisson_tensor_conf_x0 = &poisson_tensor_conf[0]; 
  double px0_q = 0.0; 
  for (int a = 0; a < 2; ++a) { 
    px0_q += vst_1x1v_tensor_p1_ho_conf_ev[i*2 + a]*poisson_tensor_conf_x0[a]; 
  } 
  double dH_dx0 = 0.0; 
  for (int a = 0; a < 2; ++a) dH_dx0 += vst_1x1v_tensor_p1_ho_ph_v0_CmDx0[i*2 + a]*G[a]; 
  return -(px0_q*dH_dx0*dx10*jacob_cx_inv); 
} 

GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 1; ++j) alpha_quad[i*1 + j] += nc_hamil_phase_ho_alpha_quad_vx_1x1v_tensor_p1_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
