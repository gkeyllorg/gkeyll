#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH int nc_hamil_vel_dense_ho_alpha_quad_vz_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 6; 
  O += off*8; 
  I += off*16; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv10 = 2.0/dxv[3]; 
  const double dv11 = 2.0/dxv[4]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[80]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[88]; 
  const double *poisson_tensor_conf_4 = &poisson_tensor_conf[104]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[112]; 
  const double *poisson_tensor_conf_7 = &poisson_tensor_conf[128]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[136]; 
  if (tid < 8) { 
    const int i = tid; 
    double p1_q = 0.0; 
    double p2_q = 0.0; 
    double p4_q = 0.0; 
    double p5_q = 0.0; 
    double p7_q = 0.0; 
    double p8_q = 0.0; 
    for (int a = 0; a < 8; ++a) { 
      p1_q += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*poisson_tensor_conf_1[a]; 
      p2_q += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*poisson_tensor_conf_2[a]; 
      p4_q += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*poisson_tensor_conf_4[a]; 
      p5_q += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*poisson_tensor_conf_5[a]; 
      p7_q += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*poisson_tensor_conf_7[a]; 
      p8_q += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*poisson_tensor_conf_8[a]; 
    } 
    O[0*8 + i] = -p1_q; 
    O[1*8 + i] = -p4_q; 
    O[2*8 + i] = -p7_q; 
    O[3*8 + i] = -p2_q; 
    O[4*8 + i] = -p5_q; 
    O[5*8 + i] = -p8_q; 
  } 
  if (tid < 16) { 
    const int j = tid; 
    double dH_dv0 = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dv0 += vst_3x3v_tensor_p1_ho_vel_dv0_v2[j*27 + b]*hamil[b]; 
    double dH_dv1 = 0.0; 
    for (int b = 0; b < 27; ++b) dH_dv1 += vst_3x3v_tensor_p1_ho_vel_dv1_v2[j*27 + b]*hamil[b]; 
    const double jacob_vx_inv = 1.0/jacob_vel_surf[0 + j/4]; 
    const double jacob_vy_inv = 1.0/jacob_vel_surf[4 + j%4]; 
    const double xn0 = vst_3x3v_tensor_p1_ho_vel_nodes_v2[j*2 + 0]; 
    const double xn0_sq = xn0*xn0; 
    const double vt1 = 4.677071733467426*vmap_v0[3]*xn0*xn0_sq+2.371708245126284*vmap_v0[2]*xn0_sq-2.806243040080455*vmap_v0[3]*xn0+1.224744871391589*vmap_v0[1]*xn0-0.7905694150420947*vmap_v0[2]+0.7071067811865475*vmap_v0[0]; 
    const double xn1 = vst_3x3v_tensor_p1_ho_vel_nodes_v2[j*2 + 1]; 
    const double xn1_sq = xn1*xn1; 
    const double vt2 = 4.677071733467426*vmap_v1[3]*xn1*xn1_sq+2.371708245126284*vmap_v1[2]*xn1_sq-2.806243040080455*vmap_v1[3]*xn1+1.224744871391589*vmap_v1[1]*xn1-0.7905694150420947*vmap_v1[2]+0.7071067811865475*vmap_v1[0]; 
    const double vt3 = -(1.8708286933869707*vmap_v2[3])+1.5811388300841895*vmap_v2[2]-1.224744871391589*vmap_v2[1]+0.7071067811865475*vmap_v2[0]; 
    I[0*16 + j] = vt1*dH_dv0*dv10*jacob_vx_inv; 
    I[1*16 + j] = vt2*dH_dv0*dv10*jacob_vx_inv; 
    I[2*16 + j] = vt3*dH_dv0*dv10*jacob_vx_inv; 
    I[3*16 + j] = vt1*dH_dv1*dv11*jacob_vy_inv; 
    I[4*16 + j] = vt2*dH_dv1*dv11*jacob_vy_inv; 
    I[5*16 + j] = vt3*dH_dv1*dv11*jacob_vy_inv; 
  } 
  return 6; 
} 

GKYL_CU_DH void nc_hamil_vel_dense_ho_alpha_quad_vz_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[48]; 
  double I[96]; 
  for (int tid = 0; tid < 16; ++tid) nc_hamil_vel_dense_ho_alpha_quad_vz_3x3v_tensor_p1_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 6; ++t) alpha += O[t*8 + i]*I[t*16 + j]; 
      alpha_quad[i*16 + j] += alpha; 
    } 
  } 
} 
