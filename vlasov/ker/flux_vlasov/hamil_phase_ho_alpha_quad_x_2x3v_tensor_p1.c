#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH int hamil_phase_ho_alpha_quad_x_2x3v_tensor_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 6; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double dv12 = 2.0/dxv[4]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[4]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[8]; 
  for (int i = tid; i < 2; i += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 4; ++a) P0 += vst_2x3v_tensor_p1_ho_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 2; ++a) O[(0 + a)*2 + i] = P0*vst_2x3v_tensor_p1_ho_ph_x0_Cm[i*2 + a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 4; ++a) P1 += vst_2x3v_tensor_p1_ho_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_1[a]; 
      for (int a = 0; a < 2; ++a) O[(2 + a)*2 + i] = P1*vst_2x3v_tensor_p1_ho_ph_x0_Cm[i*2 + a]; 
      double P2 = 0.0; 
      for (int a = 0; a < 4; ++a) P2 += vst_2x3v_tensor_p1_ho_confsurf_x0_ev_r[i*4 + a]*poisson_tensor_conf_2[a]; 
      for (int a = 0; a < 2; ++a) O[(4 + a)*2 + i] = P2*vst_2x3v_tensor_p1_ho_ph_x0_Cm[i*2 + a]; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 4; ++a) P0 += vst_2x3v_tensor_p1_ho_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 2; ++a) O[(0 + a)*2 + i] = P0*vst_2x3v_tensor_p1_ho_ph_x0_Cm[i*2 + a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 4; ++a) P1 += vst_2x3v_tensor_p1_ho_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_1[a]; 
      for (int a = 0; a < 2; ++a) O[(2 + a)*2 + i] = P1*vst_2x3v_tensor_p1_ho_ph_x0_Cm[i*2 + a]; 
      double P2 = 0.0; 
      for (int a = 0; a < 4; ++a) P2 += vst_2x3v_tensor_p1_ho_confsurf_x0_ev_l[i*4 + a]*poisson_tensor_conf_2[a]; 
      for (int a = 0; a < 2; ++a) O[(4 + a)*2 + i] = P2*vst_2x3v_tensor_p1_ho_ph_x0_Cm[i*2 + a]; 
    } 
  } 
  for (int j = tid; j < 64; j += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double G0[2]; 
      for (int a = 0; a < 2; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 108; ++k) { 
        G0[vst_2x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_2x3v_tensor_p1_ho_ph_x0_Vd0[j*37 + vst_2x3v_tensor_p1_ho_ph_x0_vrd0map[k]]*(vst_2x3v_tensor_p1_ho_ph_x0_dcoefr0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(0 + a)*64 + j] = G0[a]*dv10*(1.0/jacob_vel_surf[0 + j/16]); 
      double G1[2]; 
      for (int a = 0; a < 2; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 108; ++k) { 
        G1[vst_2x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_2x3v_tensor_p1_ho_ph_x0_Vd1[j*37 + vst_2x3v_tensor_p1_ho_ph_x0_vrd1map[k]]*(vst_2x3v_tensor_p1_ho_ph_x0_dcoefr1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(2 + a)*64 + j] = G1[a]*dv11*(1.0/jacob_vel_surf[4 + j/4%4]); 
      double G2[2]; 
      for (int a = 0; a < 2; ++a) G2[a] = 0.0; 
      for (int k = 0; k < 108; ++k) { 
        G2[vst_2x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_2x3v_tensor_p1_ho_ph_x0_Vd2[j*37 + vst_2x3v_tensor_p1_ho_ph_x0_vrd2map[k]]*(vst_2x3v_tensor_p1_ho_ph_x0_dcoefr2[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(4 + a)*64 + j] = G2[a]*dv12*(1.0/jacob_vel_surf[8 + j%4]); 
    } 
    else if (hamil_pt_edge == 1) { 
      double G0[2]; 
      for (int a = 0; a < 2; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 108; ++k) { 
        G0[vst_2x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_2x3v_tensor_p1_ho_ph_x0_Vd0[j*37 + vst_2x3v_tensor_p1_ho_ph_x0_vld0map[k]]*(vst_2x3v_tensor_p1_ho_ph_x0_dcoefl0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(0 + a)*64 + j] = G0[a]*dv10*(1.0/jacob_vel_surf[0 + j/16]); 
      double G1[2]; 
      for (int a = 0; a < 2; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 108; ++k) { 
        G1[vst_2x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_2x3v_tensor_p1_ho_ph_x0_Vd1[j*37 + vst_2x3v_tensor_p1_ho_ph_x0_vld1map[k]]*(vst_2x3v_tensor_p1_ho_ph_x0_dcoefl1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(2 + a)*64 + j] = G1[a]*dv11*(1.0/jacob_vel_surf[4 + j/4%4]); 
      double G2[2]; 
      for (int a = 0; a < 2; ++a) G2[a] = 0.0; 
      for (int k = 0; k < 108; ++k) { 
        G2[vst_2x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_2x3v_tensor_p1_ho_ph_x0_Vd2[j*37 + vst_2x3v_tensor_p1_ho_ph_x0_vld2map[k]]*(vst_2x3v_tensor_p1_ho_ph_x0_dcoefl2[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(4 + a)*64 + j] = G2[a]*dv12*(1.0/jacob_vel_surf[8 + j%4]); 
    } 
  } 
  return 6; 
} 

GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_2x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[12]; 
  double I[384]; 
  hamil_phase_ho_alpha_quad_x_2x3v_tensor_p1_shared(0, 1, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 64; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 6; ++t) alpha += O[t*2 + i]*I[t*64 + j]; 
      alpha_quad[i*64 + j] += alpha; 
    } 
  } 
} 
