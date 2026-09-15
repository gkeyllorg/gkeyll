#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p1.h> 
GKYL_CU_DH double nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double G[4]; 
  double Gd1[4]; 
  for (int a = 0; a < 4; ++a) { G[a] = 0.0; Gd1[a] = 0.0; } 
  for (int k = 0; k < 16; ++k) { 
    const int a = vst_2x2v_ser_p1_ph_v0_cmap[k]; 
    G[a] += vst_2x2v_ser_p1_ph_v0_V[j*4 + vst_2x2v_ser_p1_ph_v0_vrmap[k]]*(vst_2x2v_ser_p1_ph_v0_coefr[k]*hamil[k]); 
    Gd1[a] += vst_2x2v_ser_p1_ph_v0_Vd1[j*3 + vst_2x2v_ser_p1_ph_v0_vrd1map[k]]*(vst_2x2v_ser_p1_ph_v0_dcoefr1[k]*hamil[k]); 
  } 
  double dH_dv1 = 0.0; 
  for (int a = 0; a < 4; ++a) dH_dv1 += vst_2x2v_ser_p1_ph_v0_Cm[i*4 + a]*Gd1[a]; 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  const double *poisson_tensor_conf_x0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_x2 = &poisson_tensor_conf[8]; 
  double px0_q = 0.0; 
  double px1_q = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    px0_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_x0[a]; 
    px1_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_x2[a]; 
  } 
  double dH_dx0 = 0.0; 
  for (int a = 0; a < 4; ++a) dH_dx0 += vst_2x2v_ser_p1_ph_v0_CmDx0[i*4 + a]*G[a]; 
  double dH_dx1 = 0.0; 
  for (int a = 0; a < 4; ++a) dH_dx1 += vst_2x2v_ser_p1_ph_v0_CmDx1[i*4 + a]*G[a]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double dv11 = 2.0/dxv[3]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[3]; 
  const double vt1 = 0.7071067811865475*vmap_v0[0] - 1.224744871391589*vmap_v0[1]; 
  const double vt2 = 0.7071067811865475*vmap_v1[0] + 1.224744871391589*vmap_v1[1]*vst_2x2v_ser_p1_vel_nodes_v0[j*1 + 0]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[20]; 
  double p0_q = 0.0; 
  double p1_q = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    p0_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_0[a]; 
    p1_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_1[a]; 
  } 
  const double om12 = vt1*p0_q + vt2*p1_q; 
  return -(px0_q*dH_dx0*dx10*jacob_cx_inv + px1_q*dH_dx1*dx11*jacob_cy_inv) + om12*dH_dv1*dv11*jacob_vy_inv; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 2; ++j) alpha_quad[i*2 + j] += nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
