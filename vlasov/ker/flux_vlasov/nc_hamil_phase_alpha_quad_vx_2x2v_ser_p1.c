#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p1.h> 
GKYL_CU_DH int nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 12; 
  O += off*4; 
  I += off*2; 
  const double dx10 = 2.0/dxv[0]; 
  const double dx11 = 2.0/dxv[1]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double jacob_cy_inv = 1.0/jacob_pos[2]; 
  const double *poisson_tensor_conf_x0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_x2 = &poisson_tensor_conf[8]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double dv11 = 2.0/dxv[3]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[3]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[20]; 
  if (tid < 4) { 
    const int i = tid; 
    double px0_q = 0.0; 
    double px1_q = 0.0; 
    for (int a = 0; a < 4; ++a) { 
      px0_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_x0[a]; 
      px1_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_x2[a]; 
    } 
    for (int a = 0; a < 4; ++a) O[a*4 + i] = -(px0_q*vst_2x2v_ser_p1_ph_v0_CmDx0[i*4 + a]*dx10*jacob_cx_inv + px1_q*vst_2x2v_ser_p1_ph_v0_CmDx1[i*4 + a]*dx11*jacob_cy_inv); 
    double p0_q = 0.0; 
    double p1_q = 0.0; 
    for (int a = 0; a < 4; ++a) { 
      p0_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_0[a]; 
      p1_q += vst_2x2v_ser_p1_conf_ev[i*4 + a]*poisson_tensor_conf_1[a]; 
    } 
    for (int a = 0; a < 4; ++a) O[(4 + a)*4 + i] = p0_q*vst_2x2v_ser_p1_ph_v0_Cm[i*4 + a]; 
    for (int a = 0; a < 4; ++a) O[(8 + a)*4 + i] = p1_q*vst_2x2v_ser_p1_ph_v0_Cm[i*4 + a]; 
  } 
  if (tid < 2) { 
    const int j = tid; 
    double G[4]; 
    double Gd1[4]; 
    for (int a = 0; a < 4; ++a) { G[a] = 0.0; Gd1[a] = 0.0; } 
    for (int k = 0; k < 16; ++k) { 
      const int a = vst_2x2v_ser_p1_ph_v0_cmap[k]; 
      G[a] += vst_2x2v_ser_p1_ph_v0_V[j*4 + vst_2x2v_ser_p1_ph_v0_vrmap[k]]*(vst_2x2v_ser_p1_ph_v0_coefr[k]*hamil[k]); 
      Gd1[a] += vst_2x2v_ser_p1_ph_v0_Vd1[j*3 + vst_2x2v_ser_p1_ph_v0_vrd1map[k]]*(vst_2x2v_ser_p1_ph_v0_dcoefr1[k]*hamil[k]); 
    } 
    for (int a = 0; a < 4; ++a) I[a*2 + j] = G[a]; 
    const double vt1 = 0.7071067811865475*vmap_v0[0] - 1.224744871391589*vmap_v0[1]; 
    const double vt2 = 0.7071067811865475*vmap_v1[0] + 1.224744871391589*vmap_v1[1]*vst_2x2v_ser_p1_vel_nodes_v0[j*1 + 0]; 
    for (int a = 0; a < 4; ++a) I[(4 + a)*2 + j] = vt1*Gd1[a]*dv11*jacob_vy_inv; 
    for (int a = 0; a < 4; ++a) I[(8 + a)*2 + j] = vt2*Gd1[a]*dv11*jacob_vy_inv; 
  } 
  return 12; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[48]; 
  double I[24]; 
  for (int tid = 0; tid < 4; ++tid) nc_hamil_phase_alpha_quad_vx_2x2v_ser_p1_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 2; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 12; ++t) alpha += O[t*4 + i]*I[t*2 + j]; 
      alpha_quad[i*2 + j] += alpha; 
    } 
  } 
} 
