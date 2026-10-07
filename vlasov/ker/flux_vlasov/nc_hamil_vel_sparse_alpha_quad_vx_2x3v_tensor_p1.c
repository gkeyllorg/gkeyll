#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH int nc_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 6; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[36]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[40]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[48]; 
  const double *poisson_tensor_conf_4 = &poisson_tensor_conf[52]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[60]; 
  const double *poisson_tensor_conf_7 = &poisson_tensor_conf[64]; 
  for (int i = tid; i < 4; i += nthreads) { 
    double p0_q = 0.0; 
    double p1_q = 0.0; 
    double p3_q = 0.0; 
    double p4_q = 0.0; 
    double p6_q = 0.0; 
    double p7_q = 0.0; 
    for (int a = 0; a < 4; ++a) { 
      p0_q += vst_2x3v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_0[a]; 
      p1_q += vst_2x3v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_1[a]; 
      p3_q += vst_2x3v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_3[a]; 
      p4_q += vst_2x3v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_4[a]; 
      p6_q += vst_2x3v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_6[a]; 
      p7_q += vst_2x3v_tensor_p1_conf_ev[i*4 + a]*poisson_tensor_conf_7[a]; 
    } 
    O[0*4 + i] = p0_q; 
    O[1*4 + i] = p3_q; 
    O[2*4 + i] = p6_q; 
    O[3*4 + i] = p1_q; 
    O[4*4 + i] = p4_q; 
    O[5*4 + i] = p7_q; 
  } 
  for (int j = tid; j < 9; j += nthreads) { 
    double dH_dv1 = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_2x3v_tensor_p1_vel_sparse_idx[s]; 
      dH_dv1 += vst_2x3v_tensor_p1_vel_dv1_v0[j*27 + b]*hamil[b]; 
    } 
    double dH_dv2 = 0.0; 
    for (int s = 0; s < 7; ++s) { 
      const int b = vst_2x3v_tensor_p1_vel_sparse_idx[s]; 
      dH_dv2 += vst_2x3v_tensor_p1_vel_dv2_v0[j*27 + b]*hamil[b]; 
    } 
    const double jacob_vy_inv = 1.0/jacob_vel_surf[3 + j/3]; 
    const double jacob_vz_inv = 1.0/jacob_vel_surf[6 + j%3]; 
    const double vt1 = -(1.8708286933869707*vmap_v0[3])+1.5811388300841895*vmap_v0[2]-1.224744871391589*vmap_v0[1]+0.7071067811865475*vmap_v0[0]; 
    const double xn1 = vst_2x3v_tensor_p1_vel_nodes_v0[j*2 + 0]; 
    const double xn1_sq = xn1*xn1; 
    const double vt2 = 4.677071733467426*vmap_v1[3]*xn1*xn1_sq+2.371708245126284*vmap_v1[2]*xn1_sq-2.806243040080455*vmap_v1[3]*xn1+1.224744871391589*vmap_v1[1]*xn1-0.7905694150420947*vmap_v1[2]+0.7071067811865475*vmap_v1[0]; 
    const double xn2 = vst_2x3v_tensor_p1_vel_nodes_v0[j*2 + 1]; 
    const double xn2_sq = xn2*xn2; 
    const double vt3 = 4.677071733467426*vmap_v2[3]*xn2*xn2_sq+2.371708245126284*vmap_v2[2]*xn2_sq-2.806243040080455*vmap_v2[3]*xn2+1.224744871391589*vmap_v2[1]*xn2-0.7905694150420947*vmap_v2[2]+0.7071067811865475*vmap_v2[0]; 
    I[0*9 + j] = vt1*dH_dv1*dv11*jacob_vy_inv; 
    I[1*9 + j] = vt2*dH_dv1*dv11*jacob_vy_inv; 
    I[2*9 + j] = vt3*dH_dv1*dv11*jacob_vy_inv; 
    I[3*9 + j] = vt1*dH_dv2*dv12*jacob_vz_inv; 
    I[4*9 + j] = vt2*dH_dv2*dv12*jacob_vz_inv; 
    I[5*9 + j] = vt3*dH_dv2*dv12*jacob_vz_inv; 
  } 
  return 6; 
} 

GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[24]; 
  double I[54]; 
  nc_hamil_vel_sparse_alpha_quad_vx_2x3v_tensor_p1_shared(0, 1, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 6; ++t) alpha += O[t*4 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
