#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p1.h> 
GKYL_CU_DH int nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 6; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv12 = 2.0/dxv[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double jacob_vz_inv = 1.0/jacob_vel_surf[6]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[18]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[22]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[24]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[28]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[30]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[34]; 
  for (int i = tid; i < 2; i += nthreads) { 
    double p0_q = 0.0; 
    double p2_q = 0.0; 
    double p3_q = 0.0; 
    double p5_q = 0.0; 
    double p6_q = 0.0; 
    double p8_q = 0.0; 
    for (int a = 0; a < 2; ++a) { 
      p0_q += vst_1x3v_ser_p1_conf_ev[i*2 + a]*poisson_tensor_conf_0[a]; 
      p2_q += vst_1x3v_ser_p1_conf_ev[i*2 + a]*poisson_tensor_conf_2[a]; 
      p3_q += vst_1x3v_ser_p1_conf_ev[i*2 + a]*poisson_tensor_conf_3[a]; 
      p5_q += vst_1x3v_ser_p1_conf_ev[i*2 + a]*poisson_tensor_conf_5[a]; 
      p6_q += vst_1x3v_ser_p1_conf_ev[i*2 + a]*poisson_tensor_conf_6[a]; 
      p8_q += vst_1x3v_ser_p1_conf_ev[i*2 + a]*poisson_tensor_conf_8[a]; 
    } 
    O[0*2 + i] = -p0_q; 
    O[1*2 + i] = -p3_q; 
    O[2*2 + i] = -p6_q; 
    O[3*2 + i] = p2_q; 
    O[4*2 + i] = p5_q; 
    O[5*2 + i] = p8_q; 
  } 
  for (int j = tid; j < 4; j += nthreads) { 
    double dH_dv0 = 0.0; 
    for (int s = 0; s < 4; ++s) { 
      const int b = vst_1x3v_ser_p1_vel_sparse_idx[s]; 
      dH_dv0 += vst_1x3v_ser_p1_vel_dv0_v1[j*8 + b]*hamil[b]; 
    } 
    double dH_dv2 = 0.0; 
    for (int s = 0; s < 4; ++s) { 
      const int b = vst_1x3v_ser_p1_vel_sparse_idx[s]; 
      dH_dv2 += vst_1x3v_ser_p1_vel_dv2_v1[j*8 + b]*hamil[b]; 
    } 
    const double vt1 = 0.7071067811865475*vmap_v0[0] + 1.224744871391589*vmap_v0[1]*vst_1x3v_ser_p1_vel_nodes_v1[j*2 + 0]; 
    const double vt2 = 0.7071067811865475*vmap_v1[0] - 1.224744871391589*vmap_v1[1]; 
    const double vt3 = 0.7071067811865475*vmap_v2[0] + 1.224744871391589*vmap_v2[1]*vst_1x3v_ser_p1_vel_nodes_v1[j*2 + 1]; 
    I[0*4 + j] = vt1*dH_dv0*dv10*jacob_vx_inv; 
    I[1*4 + j] = vt2*dH_dv0*dv10*jacob_vx_inv; 
    I[2*4 + j] = vt3*dH_dv0*dv10*jacob_vx_inv; 
    I[3*4 + j] = vt1*dH_dv2*dv12*jacob_vz_inv; 
    I[4*4 + j] = vt2*dH_dv2*dv12*jacob_vz_inv; 
    I[5*4 + j] = vt3*dH_dv2*dv12*jacob_vz_inv; 
  } 
  return 6; 
} 

GKYL_CU_DH void nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[12]; 
  double I[24]; 
  nc_hamil_vel_sparse_alpha_quad_vy_1x3v_ser_p1_shared(0, 1, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 6; ++t) alpha += O[t*2 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
