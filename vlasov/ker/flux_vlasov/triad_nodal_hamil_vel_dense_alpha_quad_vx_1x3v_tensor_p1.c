#include <gkyl_vlasov_triad_nodal_kernels.h> 
#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_tensor_p1.h> 
GKYL_CU_DH int triad_nodal_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 6; 
  O += off*2; 
  I += off*9; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  if (tid < 2) { 
    const int i = tid; 
    const double p0_q = poisson_tensor_conf[i*18 + 9]; 
    const double p1_q = poisson_tensor_conf[i*18 + 10]; 
    const double p3_q = poisson_tensor_conf[i*18 + 12]; 
    const double p4_q = poisson_tensor_conf[i*18 + 13]; 
    const double p6_q = poisson_tensor_conf[i*18 + 15]; 
    const double p7_q = poisson_tensor_conf[i*18 + 16]; 
    O[0*2 + i] = p0_q; 
    O[1*2 + i] = p3_q; 
    O[2*2 + i] = p6_q; 
    O[3*2 + i] = p1_q; 
    O[4*2 + i] = p4_q; 
    O[5*2 + i] = p7_q; 
  } 
  if (tid < 9) { 
    const int j = tid; 
    const double jacob_vy_inv = 1.0/jacob_vel_surf[3 + j/3]; 
    const double jacob_vz_inv = 1.0/jacob_vel_surf[6 + j%3]; 
    const double vt1 = -(1.8708286933869707*vmap_v0[3])+1.5811388300841895*vmap_v0[2]-1.224744871391589*vmap_v0[1]+0.7071067811865475*vmap_v0[0]; 
    const double xn1 = vst_1x3v_tensor_p1_vel_nodes_v0[j*2 + 0]; 
    const double xn1_sq = xn1*xn1; 
    const double vt2 = 4.677071733467426*vmap_v1[3]*xn1*xn1_sq+2.371708245126284*vmap_v1[2]*xn1_sq-2.806243040080455*vmap_v1[3]*xn1+1.224744871391589*vmap_v1[1]*xn1-0.7905694150420947*vmap_v1[2]+0.7071067811865475*vmap_v1[0]; 
    const double xn2 = vst_1x3v_tensor_p1_vel_nodes_v0[j*2 + 1]; 
    const double xn2_sq = xn2*xn2; 
    const double vt3 = 4.677071733467426*vmap_v2[3]*xn2*xn2_sq+2.371708245126284*vmap_v2[2]*xn2_sq-2.806243040080455*vmap_v2[3]*xn2+1.224744871391589*vmap_v2[1]*xn2-0.7905694150420947*vmap_v2[2]+0.7071067811865475*vmap_v2[0]; 
    double dH_dv1 = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dv1 += vst_1x3v_tensor_p1_vel_dv1_v0[j*27 + b]*hamil[b]; 
    double dH_dv2 = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dv2 += vst_1x3v_tensor_p1_vel_dv2_v0[j*27 + b]*hamil[b]; 
    I[0*9 + j] = vt1*dH_dv1*dv11*jacob_vy_inv; 
    I[1*9 + j] = vt2*dH_dv1*dv11*jacob_vy_inv; 
    I[2*9 + j] = vt3*dH_dv1*dv11*jacob_vy_inv; 
    I[3*9 + j] = vt1*dH_dv2*dv12*jacob_vz_inv; 
    I[4*9 + j] = vt2*dH_dv2*dv12*jacob_vz_inv; 
    I[5*9 + j] = vt3*dH_dv2*dv12*jacob_vz_inv; 
  } 
  return 6; 
} 

GKYL_CU_DH void triad_nodal_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[12]; 
  double I[54]; 
  for (int tid = 0; tid < 9; ++tid) triad_nodal_hamil_vel_dense_alpha_quad_vx_1x3v_tensor_p1_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 6; ++t) alpha += O[t*2 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
