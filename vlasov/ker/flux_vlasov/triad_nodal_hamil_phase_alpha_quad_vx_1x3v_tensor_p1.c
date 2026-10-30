#include <gkyl_vlasov_triad_nodal_kernels.h> 
#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_tensor_p1.h> 
GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_vx_1x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 14; 
  O += off*2; 
  I += off*9; 
  const double dx10 = 2.0/dxv[0]; 
  const double jacob_cx_inv = 1.0/jacob_pos[0]; 
  const double *vmap_v0 = &vmap[0]; 
  const double *vmap_v1 = &vmap[4]; 
  const double *vmap_v2 = &vmap[8]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  if (tid < 2) { 
    const int i = tid; 
    const double px0_q = poisson_tensor_conf[i*18 + 0]; 
    for (int a = 0; a < 2; ++a) O[a*2 + i] = -(px0_q*vst_1x3v_tensor_p1_ph_v0_CmDx0[i*2 + a]*dx10*jacob_cx_inv); 
    const double p0_q = poisson_tensor_conf[i*18 + 9]; 
    const double p1_q = poisson_tensor_conf[i*18 + 10]; 
    const double p3_q = poisson_tensor_conf[i*18 + 12]; 
    const double p4_q = poisson_tensor_conf[i*18 + 13]; 
    const double p6_q = poisson_tensor_conf[i*18 + 15]; 
    const double p7_q = poisson_tensor_conf[i*18 + 16]; 
    for (int a = 0; a < 2; ++a) O[(2 + a)*2 + i] = p0_q*vst_1x3v_tensor_p1_ph_v0_Cm[i*2 + a]; 
    for (int a = 0; a < 2; ++a) O[(4 + a)*2 + i] = p3_q*vst_1x3v_tensor_p1_ph_v0_Cm[i*2 + a]; 
    for (int a = 0; a < 2; ++a) O[(6 + a)*2 + i] = p6_q*vst_1x3v_tensor_p1_ph_v0_Cm[i*2 + a]; 
    for (int a = 0; a < 2; ++a) O[(8 + a)*2 + i] = p1_q*vst_1x3v_tensor_p1_ph_v0_Cm[i*2 + a]; 
    for (int a = 0; a < 2; ++a) O[(10 + a)*2 + i] = p4_q*vst_1x3v_tensor_p1_ph_v0_Cm[i*2 + a]; 
    for (int a = 0; a < 2; ++a) O[(12 + a)*2 + i] = p7_q*vst_1x3v_tensor_p1_ph_v0_Cm[i*2 + a]; 
  } 
  if (tid < 9) { 
    const int j = tid; 
    double G[2]; 
    double Gd1[2]; 
    double Gd2[2]; 
    for (int a = 0; a < 2; ++a) { G[a] = 0.0; Gd1[a] = 0.0; Gd2[a] = 0.0; } 
    for (int k = 0; k < 54; ++k) { 
      const int a = vst_1x3v_tensor_p1_ph_v0_cmap[k]; 
      G[a] += vst_1x3v_tensor_p1_ph_v0_V[j*27 + vst_1x3v_tensor_p1_ph_v0_vrmap[k]]*(vst_1x3v_tensor_p1_ph_v0_coefr[k]*hamil[k]); 
      Gd1[a] += vst_1x3v_tensor_p1_ph_v0_Vd1[j*19 + vst_1x3v_tensor_p1_ph_v0_vrd1map[k]]*(vst_1x3v_tensor_p1_ph_v0_dcoefr1[k]*hamil[k]); 
      Gd2[a] += vst_1x3v_tensor_p1_ph_v0_Vd2[j*19 + vst_1x3v_tensor_p1_ph_v0_vrd2map[k]]*(vst_1x3v_tensor_p1_ph_v0_dcoefr2[k]*hamil[k]); 
    } 
    for (int a = 0; a < 2; ++a) I[a*9 + j] = G[a]; 
    const double jacob_vy_inv = 1.0/jacob_vel_surf[3 + j/3]; 
    const double jacob_vz_inv = 1.0/jacob_vel_surf[6 + j%3]; 
    const double vt1 = -(1.8708286933869707*vmap_v0[3])+1.5811388300841895*vmap_v0[2]-1.224744871391589*vmap_v0[1]+0.7071067811865475*vmap_v0[0]; 
    const double xn1 = vst_1x3v_tensor_p1_vel_nodes_v0[j*2 + 0]; 
    const double xn1_sq = xn1*xn1; 
    const double vt2 = 4.677071733467426*vmap_v1[3]*xn1*xn1_sq+2.371708245126284*vmap_v1[2]*xn1_sq-2.806243040080455*vmap_v1[3]*xn1+1.224744871391589*vmap_v1[1]*xn1-0.7905694150420947*vmap_v1[2]+0.7071067811865475*vmap_v1[0]; 
    const double xn2 = vst_1x3v_tensor_p1_vel_nodes_v0[j*2 + 1]; 
    const double xn2_sq = xn2*xn2; 
    const double vt3 = 4.677071733467426*vmap_v2[3]*xn2*xn2_sq+2.371708245126284*vmap_v2[2]*xn2_sq-2.806243040080455*vmap_v2[3]*xn2+1.224744871391589*vmap_v2[1]*xn2-0.7905694150420947*vmap_v2[2]+0.7071067811865475*vmap_v2[0]; 
    for (int a = 0; a < 2; ++a) I[(2 + a)*9 + j] = vt1*Gd1[a]*dv11*jacob_vy_inv; 
    for (int a = 0; a < 2; ++a) I[(4 + a)*9 + j] = vt2*Gd1[a]*dv11*jacob_vy_inv; 
    for (int a = 0; a < 2; ++a) I[(6 + a)*9 + j] = vt3*Gd1[a]*dv11*jacob_vy_inv; 
    for (int a = 0; a < 2; ++a) I[(8 + a)*9 + j] = vt1*Gd2[a]*dv12*jacob_vz_inv; 
    for (int a = 0; a < 2; ++a) I[(10 + a)*9 + j] = vt2*Gd2[a]*dv12*jacob_vz_inv; 
    for (int a = 0; a < 2; ++a) I[(12 + a)*9 + j] = vt3*Gd2[a]*dv12*jacob_vz_inv; 
  } 
  return 14; 
} 

GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_vx_1x3v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[28]; 
  double I[126]; 
  for (int tid = 0; tid < 9; ++tid) triad_nodal_hamil_phase_alpha_quad_vx_1x3v_tensor_p1_shared(tid, 0, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 14; ++t) alpha += O[t*2 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
