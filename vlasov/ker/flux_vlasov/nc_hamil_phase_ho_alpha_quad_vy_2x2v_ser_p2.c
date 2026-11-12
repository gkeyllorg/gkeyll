#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p2.h> 
GKYL_CU_DH int nc_hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 24; 
  O += off*16; 
  I += off*4; 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[3]; 
  const double *poisson_tensor_conf_x1 = &poisson_tensor_conf[8]; 
  const double *poisson_tensor_conf_x3 = &poisson_tensor_conf[24]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double dv10 = 2.0/dxv[2]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[32]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[40]; 
  if (tid < 16) { 
    const int i = tid; 
    double px0_q = 0.0; 
    double px1_q = 0.0; 
    for (int a = 0; a < 8; ++a) { 
      px0_q += vst_2x2v_ser_p2_ho_conf_ev[i*8 + a]*poisson_tensor_conf_x1[a]; 
      px1_q += vst_2x2v_ser_p2_ho_conf_ev[i*8 + a]*poisson_tensor_conf_x3[a]; 
    } 
    for (int a = 0; a < 8; ++a) O[a*16 + i] = -(px0_q*vst_2x2v_ser_p2_ho_ph_v1_CmDx0[i*8 + a]*dx10*jacob_cx_inv + px1_q*vst_2x2v_ser_p2_ho_ph_v1_CmDx1[i*8 + a]*dx11*jacob_cy_inv); 
    double p0_q = 0.0; 
    double p1_q = 0.0; 
    for (int a = 0; a < 8; ++a) { 
      p0_q += vst_2x2v_ser_p2_ho_conf_ev[i*8 + a]*poisson_tensor_conf_0[a]; 
      p1_q += vst_2x2v_ser_p2_ho_conf_ev[i*8 + a]*poisson_tensor_conf_1[a]; 
    } 
    for (int a = 0; a < 8; ++a) O[(8 + a)*16 + i] = -p0_q*vst_2x2v_ser_p2_ho_ph_v1_Cm[i*8 + a]; 
    for (int a = 0; a < 8; ++a) O[(16 + a)*16 + i] = -p1_q*vst_2x2v_ser_p2_ho_ph_v1_Cm[i*8 + a]; 
  } 
  if (tid < 4) { 
    const int j = tid; 
    double G[8]; 
    double Gd0[8]; 
    for (int a = 0; a < 8; ++a) { G[a] = 0.0; Gd0[a] = 0.0; } 
    for (int k = 0; k < 48; ++k) { 
      const int a = vst_2x2v_ser_p2_ho_ph_v1_cmap[k]; 
      G[a] += vst_2x2v_ser_p2_ho_ph_v1_V[j*8 + vst_2x2v_ser_p2_ho_ph_v1_vrmap[k]]*(vst_2x2v_ser_p2_ho_ph_v1_coefr[k]*hamil[k]); 
      Gd0[a] += vst_2x2v_ser_p2_ho_ph_v1_Vd0[j*6 + vst_2x2v_ser_p2_ho_ph_v1_vrd0map[k]]*(vst_2x2v_ser_p2_ho_ph_v1_dcoefr0[k]*hamil[k]); 
    } 
    for (int a = 0; a < 8; ++a) I[a*4 + j] = G[a]; 
    const double vt1 = 0.7071067811865475*vmap_v0[0] + 1.224744871391589*vmap_v0[1]*vst_2x2v_ser_p2_ho_vel_nodes_v1[j*1 + 0]; 
    const double vt2 = 0.7071067811865475*vmap_v1[0] - 1.224744871391589*vmap_v1[1]; 
    for (int a = 0; a < 8; ++a) I[(8 + a)*4 + j] = vt1*Gd0[a]*dv10*jacob_vx_inv; 
    for (int a = 0; a < 8; ++a) I[(16 + a)*4 + j] = vt2*Gd0[a]*dv10*jacob_vx_inv; 
  } 
  return 24; 
} 

GKYL_CU_DH void nc_hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[384]; 
  double I[96]; 
  for (int tid = 0; tid < 16; ++tid) nc_hamil_phase_ho_alpha_quad_vy_2x2v_ser_p2_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 24; ++t) alpha += O[t*16 + i]*I[t*4 + j]; 
      alpha_quad[i*4 + j] += alpha; 
    } 
  } 
} 
