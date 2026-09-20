#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_tensor_p1.h> 
GKYL_CU_DH double nc_hamil_phase_alpha_quad_vy_1x2v_tensor_p1_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double G[2]; 
  double Gd0[2]; 
  for (int a = 0; a < 2; ++a) { G[a] = 0.0; Gd0[a] = 0.0; } 
  for (int k = 0; k < 18; ++k) { 
    const int a = vst_1x2v_tensor_p1_ph_v1_cmap[k]; 
    G[a] += vst_1x2v_tensor_p1_ph_v1_V[j*9 + vst_1x2v_tensor_p1_ph_v1_vrmap[k]]*(vst_1x2v_tensor_p1_ph_v1_coefr[k]*hamil[k]); 
    Gd0[a] += vst_1x2v_tensor_p1_ph_v1_Vd0[j*7 + vst_1x2v_tensor_p1_ph_v1_vrd0map[k]]*(vst_1x2v_tensor_p1_ph_v1_dcoefr0[k]*hamil[k]); 
  } 
  double dH_dv0 = 0.0; 
  for (int a = 0; a < 2; ++a) dH_dv0 += vst_1x2v_tensor_p1_ph_v1_Cm[i*2 + a]*Gd0[a]; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double *poisson_tensor_conf_x1 = &poisson_tensor_conf[2]; 
  double px0_q = 0.0; 
  for (int a = 0; a < 2; ++a) { 
    px0_q += vst_1x2v_tensor_p1_conf_ev[i*2 + a]*poisson_tensor_conf_x1[a]; 
  } 
  double dH_dx0 = 0.0; 
  for (int a = 0; a < 2; ++a) dH_dx0 += vst_1x2v_tensor_p1_ph_v1_CmDx0[i*2 + a]*G[a]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double dv10 = 2.0/dxv[1]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0 + j]; 
  const double xn0 = vst_1x2v_tensor_p1_vel_nodes_v1[j*1 + 0]; 
  const double xn0_sq = xn0*xn0; 
  const double vt1 = 4.677071733467426*vmap_v0[3]*xn0*xn0_sq+2.371708245126284*vmap_v0[2]*xn0_sq-2.806243040080455*vmap_v0[3]*xn0+1.224744871391589*vmap_v0[1]*xn0-0.7905694150420947*vmap_v0[2]+0.7071067811865475*vmap_v0[0]; 
  const double vt2 = -(1.8708286933869707*vmap_v1[3])+1.5811388300841895*vmap_v1[2]-1.224744871391589*vmap_v1[1]+0.7071067811865475*vmap_v1[0]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[8]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[10]; 
  double p0_q = 0.0; 
  double p1_q = 0.0; 
  for (int a = 0; a < 2; ++a) { 
    p0_q += vst_1x2v_tensor_p1_conf_ev[i*2 + a]*poisson_tensor_conf_0[a]; 
    p1_q += vst_1x2v_tensor_p1_conf_ev[i*2 + a]*poisson_tensor_conf_1[a]; 
  } 
  const double om12 = vt1*p0_q + vt2*p1_q; 
  return -(px0_q*dH_dx0*dx10*jacob_cx_inv) + -om12*dH_dv0*dv10*jacob_vx_inv; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 3; ++j) alpha_quad[i*3 + j] += nc_hamil_phase_alpha_quad_vy_1x2v_tensor_p1_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
