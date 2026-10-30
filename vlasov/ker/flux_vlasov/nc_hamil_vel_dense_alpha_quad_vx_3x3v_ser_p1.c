#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_ser_p1.h> 
GKYL_CU_DH int nc_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 6; 
  O += off*8; 
  I += off*4; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv11 = 2.0/dxv[4]; 
  const double dv12 = 2.0/dxv[5]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[3]; 
  const double jacob_vz_inv = 1.0/jacob_vel_surf[6]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[72]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[80]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[96]; 
  const double *poisson_tensor_conf_4 = &poisson_tensor_conf[104]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[120]; 
  const double *poisson_tensor_conf_7 = &poisson_tensor_conf[128]; 
  if (tid < 8) { 
    const int i = tid; 
    double p0_q = 0.0; 
    double p1_q = 0.0; 
    double p3_q = 0.0; 
    double p4_q = 0.0; 
    double p6_q = 0.0; 
    double p7_q = 0.0; 
    for (int a = 0; a < 8; ++a) { 
      p0_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_0[a]; 
      p1_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_1[a]; 
      p3_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_3[a]; 
      p4_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_4[a]; 
      p6_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_6[a]; 
      p7_q += vst_3x3v_ser_p1_conf_ev[i*8 + a]*poisson_tensor_conf_7[a]; 
    } 
    O[0*8 + i] = p0_q; 
    O[1*8 + i] = p3_q; 
    O[2*8 + i] = p6_q; 
    O[3*8 + i] = p1_q; 
    O[4*8 + i] = p4_q; 
    O[5*8 + i] = p7_q; 
  } 
  if (tid < 4) { 
    const int j = tid; 
    double dH_dv1 = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dv1 += vst_3x3v_ser_p1_vel_dv1_v0[j*8 + b]*hamil[b]; 
    double dH_dv2 = 0.0; 
    for (int b = 0; b < 8; ++b) dH_dv2 += vst_3x3v_ser_p1_vel_dv2_v0[j*8 + b]*hamil[b]; 
    const double vt1 = 0.7071067811865475*vmap_v0[0] - 1.224744871391589*vmap_v0[1]; 
    const double vt2 = 0.7071067811865475*vmap_v1[0] + 1.224744871391589*vmap_v1[1]*vst_3x3v_ser_p1_vel_nodes_v0[j*2 + 0]; 
    const double vt3 = 0.7071067811865475*vmap_v2[0] + 1.224744871391589*vmap_v2[1]*vst_3x3v_ser_p1_vel_nodes_v0[j*2 + 1]; 
    I[0*4 + j] = vt1*dH_dv1*dv11*jacob_vy_inv; 
    I[1*4 + j] = vt2*dH_dv1*dv11*jacob_vy_inv; 
    I[2*4 + j] = vt3*dH_dv1*dv11*jacob_vy_inv; 
    I[3*4 + j] = vt1*dH_dv2*dv12*jacob_vz_inv; 
    I[4*4 + j] = vt2*dH_dv2*dv12*jacob_vz_inv; 
    I[5*4 + j] = vt3*dH_dv2*dv12*jacob_vz_inv; 
  } 
  return 6; 
} 

GKYL_CU_DH void nc_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[48]; 
  double I[24]; 
  for (int tid = 0; tid < 8; ++tid) nc_hamil_vel_dense_alpha_quad_vx_3x3v_ser_p1_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 6; ++t) alpha += O[t*8 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
