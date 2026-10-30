#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p1.h> 
GKYL_CU_DH int nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  O += off*4; 
  I += off*2; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double dv11 = 2.0/dxv[3]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[3]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[20]; 
  if (tid < 4) { 
    const int i = tid; 
    double p0_q = 0.0; 
    double p1_q = 0.0; 
    for (int a = 0; a < 4; ++a) { 
      p0_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_0[a]; 
      p1_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_1[a]; 
    } 
    O[0*4 + i] = p0_q; 
    O[1*4 + i] = p1_q; 
  } 
  if (tid < 2) { 
    const int j = tid; 
    double dH_dv1 = 0.0; 
    for (int s = 0; s < 3; ++s) { 
      const int b = vst_2x2v_ser_p1_vel_sparse_idx[s]; 
      dH_dv1 += vst_2x2v_ser_p1_vel_dv1_v0[j*4 + b]*hamil[b]; 
    } 
    const double vt1 = 0.7071067811865475*vmap_v0[0] - 1.224744871391589*vmap_v0[1]; 
    const double vt2 = 0.7071067811865475*vmap_v1[0] + 1.224744871391589*vmap_v1[1]*vst_2x2v_ser_p1_vel_nodes_v0[j*1 + 0]; 
    I[0*2 + j] = vt1*dH_dv1*dv11*jacob_vy_inv; 
    I[1*2 + j] = vt2*dH_dv1*dv11*jacob_vy_inv; 
  } 
  return 2; 
} 

GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[4]; 
  for (int tid = 0; tid < 4; ++tid) nc_hamil_vel_sparse_alpha_quad_vx_2x2v_ser_p1_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 2; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*4 + i]*I[t*2 + j]; 
      alpha_quad[i*2 + j] += alpha; 
    } 
  } 
} 
