#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH int hamil_phase_alpha_quad_x_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 4; 
  O += off*2; 
  I += off*9; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[4]; 
  if (tid < 2) { 
    const int i = tid; 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 4; ++a) P0 += vst_2x2v_tensor_p1_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 2; ++a) O[(0 + a)*2 + i] = P0*vst_2x2v_tensor_p1_ph_x0_Cm[i*2 + a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 4; ++a) P1 += vst_2x2v_tensor_p1_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_1[a]; 
      for (int a = 0; a < 2; ++a) O[(2 + a)*2 + i] = P1*vst_2x2v_tensor_p1_ph_x0_Cm[i*2 + a]; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 4; ++a) P0 += vst_2x2v_tensor_p1_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 2; ++a) O[(0 + a)*2 + i] = P0*vst_2x2v_tensor_p1_ph_x0_Cm[i*2 + a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 4; ++a) P1 += vst_2x2v_tensor_p1_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_1[a]; 
      for (int a = 0; a < 2; ++a) O[(2 + a)*2 + i] = P1*vst_2x2v_tensor_p1_ph_x0_Cm[i*2 + a]; 
    } 
  } 
  if (tid < 9) { 
    const int j = tid; 
    if (hamil_pt_edge == -1) { 
      double G0[2]; 
      for (int a = 0; a < 2; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 36; ++k) { 
        G0[vst_2x2v_tensor_p1_ph_x0_cmap[k]] += vst_2x2v_tensor_p1_ph_x0_Vd0[j*13 + vst_2x2v_tensor_p1_ph_x0_vrd0map[k]]*(vst_2x2v_tensor_p1_ph_x0_dcoefr0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(0 + a)*9 + j] = G0[a]*dv10*(1.0/jacob_vel_surf[0 + j/3]); 
      double G1[2]; 
      for (int a = 0; a < 2; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 36; ++k) { 
        G1[vst_2x2v_tensor_p1_ph_x0_cmap[k]] += vst_2x2v_tensor_p1_ph_x0_Vd1[j*13 + vst_2x2v_tensor_p1_ph_x0_vrd1map[k]]*(vst_2x2v_tensor_p1_ph_x0_dcoefr1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(2 + a)*9 + j] = G1[a]*dv11*(1.0/jacob_vel_surf[3 + j%3]); 
    } 
    else if (hamil_pt_edge == 1) { 
      double G0[2]; 
      for (int a = 0; a < 2; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 36; ++k) { 
        G0[vst_2x2v_tensor_p1_ph_x0_cmap[k]] += vst_2x2v_tensor_p1_ph_x0_Vd0[j*13 + vst_2x2v_tensor_p1_ph_x0_vld0map[k]]*(vst_2x2v_tensor_p1_ph_x0_dcoefl0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(0 + a)*9 + j] = G0[a]*dv10*(1.0/jacob_vel_surf[0 + j/3]); 
      double G1[2]; 
      for (int a = 0; a < 2; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 36; ++k) { 
        G1[vst_2x2v_tensor_p1_ph_x0_cmap[k]] += vst_2x2v_tensor_p1_ph_x0_Vd1[j*13 + vst_2x2v_tensor_p1_ph_x0_vld1map[k]]*(vst_2x2v_tensor_p1_ph_x0_dcoefl1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(2 + a)*9 + j] = G1[a]*dv11*(1.0/jacob_vel_surf[3 + j%3]); 
    } 
  } 
  return 4; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_x_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[36]; 
  for (int tid = 0; tid < 9; ++tid) hamil_phase_alpha_quad_x_2x2v_tensor_p1_shared(tid, 0, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 4; ++t) alpha += O[t*2 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
