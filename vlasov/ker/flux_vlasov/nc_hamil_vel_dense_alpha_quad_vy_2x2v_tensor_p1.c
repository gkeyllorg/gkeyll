#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH int nc_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double dv10 = 2.0/dxv[2]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[20]; 
  for (int i = tid; i < 4; i += nthreads) { 
    double p0_q = 0.0; 
    double p1_q = 0.0; 
    for (int a = 0; a < 4; ++a) { 
      p0_q += vst_2x2v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_0[a]; 
      p1_q += vst_2x2v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_1[a]; 
    } 
    O[0*4 + i] = -p0_q; 
    O[1*4 + i] = -p1_q; 
  } 
  for (int j = tid; j < 3; j += nthreads) { 
    double dH_dv0 = 0.0; 
    for (int b = 0; b < 9; ++b) dH_dv0 += vst_2x2v_tensor_p1_vel_dv0_v1[j*9 + b]*hamil[b]; 
    const double jacob_vx_inv = 1.0/jacob_vel_surf[0 + j]; 
    const double xn0 = vst_2x2v_tensor_p1_vel_nodes_v1[j*1 + 0]; 
    const double xn0_sq = xn0*xn0; 
    const double vt1 = 4.677071733467426*vmap_v0[3]*xn0*xn0_sq+2.371708245126284*vmap_v0[2]*xn0_sq-2.806243040080455*vmap_v0[3]*xn0+1.224744871391589*vmap_v0[1]*xn0-0.7905694150420947*vmap_v0[2]+0.7071067811865475*vmap_v0[0]; 
    const double vt2 = -(1.8708286933869707*vmap_v1[3])+1.5811388300841895*vmap_v1[2]-1.224744871391589*vmap_v1[1]+0.7071067811865475*vmap_v1[0]; 
    I[0*3 + j] = vt1*dH_dv0*dv10*jacob_vx_inv; 
    I[1*3 + j] = vt2*dH_dv0*dv10*jacob_vx_inv; 
  } 
  return 2; 
} 

GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[6]; 
  nc_hamil_vel_dense_alpha_quad_vy_2x2v_tensor_p1_shared(0, 1, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*4 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
