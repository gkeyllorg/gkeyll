#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH double nc_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_tensor_p1_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dH_dv1 = 0.0; 
  for (int s = 0; s < 5; ++s) { 
    const int b = vst_2x2v_tensor_p1_ho_vel_sparse_idx[s]; 
    dH_dv1 += vst_2x2v_tensor_p1_ho_vel_dv1_v0[j*9 + b]*hamil[b]; 
  } 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double dv11 = 2.0/dxv[3]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[4 + j]; 
  const double vt1 = -(1.8708286933869707*vmap_v0[3])+1.5811388300841895*vmap_v0[2]-1.224744871391589*vmap_v0[1]+0.7071067811865475*vmap_v0[0]; 
  const double xn1 = vst_2x2v_tensor_p1_ho_vel_nodes_v0[j*1 + 0]; 
  const double xn1_sq = xn1*xn1; 
  const double vt2 = 4.677071733467426*vmap_v1[3]*xn1*xn1_sq+2.371708245126284*vmap_v1[2]*xn1_sq-2.806243040080455*vmap_v1[3]*xn1+1.224744871391589*vmap_v1[1]*xn1-0.7905694150420947*vmap_v1[2]+0.7071067811865475*vmap_v1[0]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[20]; 
  double p0_q = 0.0; 
  double p1_q = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    p0_q += vst_2x2v_tensor_p1_ho_conf_ev[i*4 + a]*poisson_tensor_conf_0[a]; 
    p1_q += vst_2x2v_tensor_p1_ho_conf_ev[i*4 + a]*poisson_tensor_conf_1[a]; 
  } 
  const double om12 = vt1*p0_q + vt2*p1_q; 
  return om12*dH_dv1*dv11*jacob_vy_inv; 
} 

GKYL_CU_DH void nc_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += nc_hamil_vel_sparse_ho_alpha_quad_vx_2x2v_tensor_p1_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
