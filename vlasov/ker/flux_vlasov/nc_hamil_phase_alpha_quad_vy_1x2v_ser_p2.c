#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_ser_p2.h> 
GKYL_CU_DH int nc_hamil_phase_alpha_quad_vy_1x2v_ser_p2_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 9; 
  O += off*3; 
  I += off*3; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double *poisson_tensor_conf_x1 = &poisson_tensor_conf[3]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double dv10 = 2.0/dxv[1]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[12]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[15]; 
  if (tid < 3) { 
    const int i = tid; 
    double px0_q = 0.0; 
    for (int a = 0; a < 3; ++a) { 
      px0_q += vst_1x2v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_x1[a]; 
    } 
    for (int a = 0; a < 3; ++a) O[a*3 + i] = -(px0_q*vst_1x2v_ser_p2_ph_v1_CmDx0[i*3 + a]*dx10*jacob_cx_inv); 
    double p0_q = 0.0; 
    double p1_q = 0.0; 
    for (int a = 0; a < 3; ++a) { 
      p0_q += vst_1x2v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_0[a]; 
      p1_q += vst_1x2v_ser_p2_conf_ev[i*3 + a]*poisson_tensor_conf_1[a]; 
    } 
    for (int a = 0; a < 3; ++a) O[(3 + a)*3 + i] = -p0_q*vst_1x2v_ser_p2_ph_v1_Cm[i*3 + a]; 
    for (int a = 0; a < 3; ++a) O[(6 + a)*3 + i] = -p1_q*vst_1x2v_ser_p2_ph_v1_Cm[i*3 + a]; 
  } 
  if (tid < 3) { 
    const int j = tid; 
    double G[3]; 
    double Gd0[3]; 
    for (int a = 0; a < 3; ++a) { G[a] = 0.0; Gd0[a] = 0.0; } 
    for (int k = 0; k < 20; ++k) { 
      const int a = vst_1x2v_ser_p2_ph_v1_cmap[k]; 
      G[a] += vst_1x2v_ser_p2_ph_v1_V[j*8 + vst_1x2v_ser_p2_ph_v1_vrmap[k]]*(vst_1x2v_ser_p2_ph_v1_coefr[k]*hamil[k]); 
      Gd0[a] += vst_1x2v_ser_p2_ph_v1_Vd0[j*6 + vst_1x2v_ser_p2_ph_v1_vrd0map[k]]*(vst_1x2v_ser_p2_ph_v1_dcoefr0[k]*hamil[k]); 
    } 
    for (int a = 0; a < 3; ++a) I[a*3 + j] = G[a]; 
    const double vt1 = 0.7071067811865475*vmap_v0[0] + 1.224744871391589*vmap_v0[1]*vst_1x2v_ser_p2_vel_nodes_v1[j*1 + 0]; 
    const double vt2 = 0.7071067811865475*vmap_v1[0] - 1.224744871391589*vmap_v1[1]; 
    for (int a = 0; a < 3; ++a) I[(3 + a)*3 + j] = vt1*Gd0[a]*dv10*jacob_vx_inv; 
    for (int a = 0; a < 3; ++a) I[(6 + a)*3 + j] = vt2*Gd0[a]*dv10*jacob_vx_inv; 
  } 
  return 9; 
} 

GKYL_CU_DH void nc_hamil_phase_alpha_quad_vy_1x2v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[27]; 
  double I[27]; 
  for (int tid = 0; tid < 3; ++tid) nc_hamil_phase_alpha_quad_vy_1x2v_ser_p2_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 9; ++t) alpha += O[t*3 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
