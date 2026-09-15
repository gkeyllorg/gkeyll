#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p1.h> 
GKYL_CU_DH double nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1_node(int i, int j, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  double dH_dv0 = 0.0; 
  for (int b = 0; b < 4; ++b) dH_dv0 += vst_2x2v_ser_p1_vel_dv0_v1[j*4 + b]*hamil[b]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double dv10 = 2.0/dxv[2]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double vt1 = 0.7071067811865475*vmap_v0[0] + 1.224744871391589*vmap_v0[1]*vst_2x2v_ser_p1_vel_nodes_v1[j*1 + 0]; 
  const double vt2 = 0.7071067811865475*vmap_v1[0] - 1.224744871391589*vmap_v1[1]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[20]; 
  double p0_q = 0.0; 
  double p1_q = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    p0_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_0[a]; 
    p1_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_1[a]; 
  } 
  const double om12 = vt1*p0_q + vt2*p1_q; 
  return -om12*dH_dv0*dv10*jacob_vx_inv; 
} 

GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf,
  const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 2; ++j) alpha_quad[i*2 + j] += nc_hamil_vel_dense_alpha_quad_vy_2x2v_ser_p1_node(i, j, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
