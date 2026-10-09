#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH int nc_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 56; 
  O += off*8; 
  I += off*9; 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double dx12 = 2.0/dxv[2]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  const double jacob_cz_inv = 1.0/jacob_pos[4]; 
  const double *poisson_tensor_conf_x1 = &poisson_tensor_conf[8]; 
  const double *poisson_tensor_conf_x4 = &poisson_tensor_conf[32]; 
  const double *poisson_tensor_conf_x7 = &poisson_tensor_conf[56]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv10 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[5]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[72]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[88]; 
  const double *poisson_tensor_conf_3 = &poisson_tensor_conf[96]; 
  const double *poisson_tensor_conf_5 = &poisson_tensor_conf[112]; 
  const double *poisson_tensor_conf_6 = &poisson_tensor_conf[120]; 
  const double *poisson_tensor_conf_8 = &poisson_tensor_conf[136]; 
  if (tid < 8) { 
    const int i = tid; 
    double px0_q = 0.0; 
    double px1_q = 0.0; 
    double px2_q = 0.0; 
    for (int a = 0; a < 8; ++a) { 
      px0_q += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*poisson_tensor_conf_x1[a]; 
      px1_q += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*poisson_tensor_conf_x4[a]; 
      px2_q += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*poisson_tensor_conf_x7[a]; 
    } 
    for (int a = 0; a < 8; ++a) O[a*8 + i] = -(px0_q*vst_3x3v_tensor_p1_ph_v1_CmDx0[i*8 + a]*dx10*jacob_cx_inv + px1_q*vst_3x3v_tensor_p1_ph_v1_CmDx1[i*8 + a]*dx11*jacob_cy_inv + px2_q*vst_3x3v_tensor_p1_ph_v1_CmDx2[i*8 + a]*dx12*jacob_cz_inv); 
    double p0_q = 0.0; 
    double p2_q = 0.0; 
    double p3_q = 0.0; 
    double p5_q = 0.0; 
    double p6_q = 0.0; 
    double p8_q = 0.0; 
    for (int a = 0; a < 8; ++a) { 
      p0_q += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*poisson_tensor_conf_0[a]; 
      p2_q += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*poisson_tensor_conf_2[a]; 
      p3_q += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*poisson_tensor_conf_3[a]; 
      p5_q += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*poisson_tensor_conf_5[a]; 
      p6_q += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*poisson_tensor_conf_6[a]; 
      p8_q += vst_3x3v_tensor_p1_conf_ev[i*8 + a]*poisson_tensor_conf_8[a]; 
    } 
    for (int a = 0; a < 8; ++a) O[(8 + a)*8 + i] = -p0_q*vst_3x3v_tensor_p1_ph_v1_Cm[i*8 + a]; 
    for (int a = 0; a < 8; ++a) O[(16 + a)*8 + i] = -p3_q*vst_3x3v_tensor_p1_ph_v1_Cm[i*8 + a]; 
    for (int a = 0; a < 8; ++a) O[(24 + a)*8 + i] = -p6_q*vst_3x3v_tensor_p1_ph_v1_Cm[i*8 + a]; 
    for (int a = 0; a < 8; ++a) O[(32 + a)*8 + i] = p2_q*vst_3x3v_tensor_p1_ph_v1_Cm[i*8 + a]; 
    for (int a = 0; a < 8; ++a) O[(40 + a)*8 + i] = p5_q*vst_3x3v_tensor_p1_ph_v1_Cm[i*8 + a]; 
    for (int a = 0; a < 8; ++a) O[(48 + a)*8 + i] = p8_q*vst_3x3v_tensor_p1_ph_v1_Cm[i*8 + a]; 
  } 
  if (tid < 9) { 
    const int j = tid; 
    double G[8]; 
    double Gd0[8]; 
    double Gd2[8]; 
    for (int a = 0; a < 8; ++a) { G[a] = 0.0; Gd0[a] = 0.0; Gd2[a] = 0.0; } 
    for (int k = 0; k < 216; ++k) { 
      const int a = vst_3x3v_tensor_p1_ph_v1_cmap[k]; 
      G[a] += vst_3x3v_tensor_p1_ph_v1_V[j*27 + vst_3x3v_tensor_p1_ph_v1_vrmap[k]]*(vst_3x3v_tensor_p1_ph_v1_coefr[k]*hamil[k]); 
      Gd0[a] += vst_3x3v_tensor_p1_ph_v1_Vd0[j*19 + vst_3x3v_tensor_p1_ph_v1_vrd0map[k]]*(vst_3x3v_tensor_p1_ph_v1_dcoefr0[k]*hamil[k]); 
      Gd2[a] += vst_3x3v_tensor_p1_ph_v1_Vd2[j*19 + vst_3x3v_tensor_p1_ph_v1_vrd2map[k]]*(vst_3x3v_tensor_p1_ph_v1_dcoefr2[k]*hamil[k]); 
    } 
    for (int a = 0; a < 8; ++a) I[a*9 + j] = G[a]; 
    const double jacob_vx_inv = 1.0/jacob_vel_surf[0 + j/3]; 
    const double jacob_vz_inv = 1.0/jacob_vel_surf[6 + j%3]; 
    const double xn0 = vst_3x3v_tensor_p1_vel_nodes_v1[j*2 + 0]; 
    const double xn0_sq = xn0*xn0; 
    const double vt1 = 4.677071733467426*vmap_v0[3]*xn0*xn0_sq+2.371708245126284*vmap_v0[2]*xn0_sq-2.806243040080455*vmap_v0[3]*xn0+1.224744871391589*vmap_v0[1]*xn0-0.7905694150420947*vmap_v0[2]+0.7071067811865475*vmap_v0[0]; 
    const double vt2 = -(1.8708286933869707*vmap_v1[3])+1.5811388300841895*vmap_v1[2]-1.224744871391589*vmap_v1[1]+0.7071067811865475*vmap_v1[0]; 
    const double xn2 = vst_3x3v_tensor_p1_vel_nodes_v1[j*2 + 1]; 
    const double xn2_sq = xn2*xn2; 
    const double vt3 = 4.677071733467426*vmap_v2[3]*xn2*xn2_sq+2.371708245126284*vmap_v2[2]*xn2_sq-2.806243040080455*vmap_v2[3]*xn2+1.224744871391589*vmap_v2[1]*xn2-0.7905694150420947*vmap_v2[2]+0.7071067811865475*vmap_v2[0]; 
    for (int a = 0; a < 8; ++a) I[(8 + a)*9 + j] = vt1*Gd0[a]*dv10*jacob_vx_inv; 
    for (int a = 0; a < 8; ++a) I[(16 + a)*9 + j] = vt2*Gd0[a]*dv10*jacob_vx_inv; 
    for (int a = 0; a < 8; ++a) I[(24 + a)*9 + j] = vt3*Gd0[a]*dv10*jacob_vx_inv; 
    for (int a = 0; a < 8; ++a) I[(32 + a)*9 + j] = vt1*Gd2[a]*dv12*jacob_vz_inv; 
    for (int a = 0; a < 8; ++a) I[(40 + a)*9 + j] = vt2*Gd2[a]*dv12*jacob_vz_inv; 
    for (int a = 0; a < 8; ++a) I[(48 + a)*9 + j] = vt3*Gd2[a]*dv12*jacob_vz_inv; 
  } 
  return 56; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_3x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[448]; 
  double I[504]; 
  for (int tid = 0; tid < 9; ++tid) nc_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 56; ++t) alpha += O[t*8 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
