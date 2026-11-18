#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p2.h> 
GKYL_CU_DH int nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 2; 
  O += off*9; 
  I += off*3; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double dv10 = 2.0/dxv[2]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[32]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[40]; 
  if (tid < 9) { 
    const int i = tid; 
    double p0_q = 0.0; 
    double p1_q = 0.0; 
    for (int a = 0; a < 8; ++a) { 
      p0_q += vst_2x2v_ser_p2_conf_ev[i*8 + a]*poisson_tensor_conf_0[a]; 
      p1_q += vst_2x2v_ser_p2_conf_ev[i*8 + a]*poisson_tensor_conf_1[a]; 
    } 
    O[0*9 + i] = -p0_q; 
    O[1*9 + i] = -p1_q; 
  } 
  if (tid < 3) { 
    const int j = tid; 
    double dH_dv0 = 0.0; 
    for (int s = 0; s < 5; ++s) { 
      const int b = vst_2x2v_ser_p2_vel_sparse_idx[s]; 
      dH_dv0 += vst_2x2v_ser_p2_vel_dv0_v1[j*8 + b]*hamil[b]; 
    } 
    const double vt1 = 0.7071067811865475*vmap_v0[0] + 1.224744871391589*vmap_v0[1]*vst_2x2v_ser_p2_vel_nodes_v1[j*1 + 0]; 
    const double vt2 = 0.7071067811865475*vmap_v1[0] - 1.224744871391589*vmap_v1[1]; 
    I[0*3 + j] = vt1*dH_dv0*dv10*jacob_vx_inv; 
    I[1*3 + j] = vt2*dH_dv0*dv10*jacob_vx_inv; 
  } 
  return 2; 
} 

GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[18]; 
  double I[6]; 
  for (int tid = 0; tid < 9; ++tid) nc_hamil_vel_sparse_alpha_quad_vy_2x2v_ser_p2_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 9; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 2; ++t) alpha += O[t*9 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
