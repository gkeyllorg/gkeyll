#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p2.h> 
GKYL_CU_DH int hamil_phase_alpha_quad_x_2x2v_ser_p2_shared(int tid, int nthreads, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 6; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[4]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[8]; 
  for (int i = tid; i < 3; i += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_2x2v_ser_p2_confsurf_x0_ev_r[i*8 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 3; ++a) O[(0 + a)*3 + i] = P0*vst_2x2v_ser_p2_ph_x0_Cm[i*3 + a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_2x2v_ser_p2_confsurf_x0_ev_r[i*8 + a]*poisson_tensor_conf_1[a]; 
      for (int a = 0; a < 3; ++a) O[(3 + a)*3 + i] = P1*vst_2x2v_ser_p2_ph_x0_Cm[i*3 + a]; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_2x2v_ser_p2_confsurf_x0_ev_l[i*8 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 3; ++a) O[(0 + a)*3 + i] = P0*vst_2x2v_ser_p2_ph_x0_Cm[i*3 + a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_2x2v_ser_p2_confsurf_x0_ev_l[i*8 + a]*poisson_tensor_conf_1[a]; 
      for (int a = 0; a < 3; ++a) O[(3 + a)*3 + i] = P1*vst_2x2v_ser_p2_ph_x0_Cm[i*3 + a]; 
    } 
  } 
  for (int j = tid; j < 9; j += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double G0[3]; 
      for (int a = 0; a < 3; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 48; ++k) { 
        G0[vst_2x2v_ser_p2_ph_x0_cmap[k]] += vst_2x2v_ser_p2_ph_x0_Vd0[j*13 + vst_2x2v_ser_p2_ph_x0_vrd0map[k]]*(vst_2x2v_ser_p2_ph_x0_dcoefr0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 3; ++a) I[(0 + a)*9 + j] = G0[a]*dv10*jacob_vx_inv; 
      double G1[3]; 
      for (int a = 0; a < 3; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 48; ++k) { 
        G1[vst_2x2v_ser_p2_ph_x0_cmap[k]] += vst_2x2v_ser_p2_ph_x0_Vd1[j*13 + vst_2x2v_ser_p2_ph_x0_vrd1map[k]]*(vst_2x2v_ser_p2_ph_x0_dcoefr1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 3; ++a) I[(3 + a)*9 + j] = G1[a]*dv11*jacob_vy_inv; 
    } 
    else if (hamil_pt_edge == 1) { 
      double G0[3]; 
      for (int a = 0; a < 3; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 48; ++k) { 
        G0[vst_2x2v_ser_p2_ph_x0_cmap[k]] += vst_2x2v_ser_p2_ph_x0_Vd0[j*13 + vst_2x2v_ser_p2_ph_x0_vld0map[k]]*(vst_2x2v_ser_p2_ph_x0_dcoefl0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 3; ++a) I[(0 + a)*9 + j] = G0[a]*dv10*jacob_vx_inv; 
      double G1[3]; 
      for (int a = 0; a < 3; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 48; ++k) { 
        G1[vst_2x2v_ser_p2_ph_x0_cmap[k]] += vst_2x2v_ser_p2_ph_x0_Vd1[j*13 + vst_2x2v_ser_p2_ph_x0_vld1map[k]]*(vst_2x2v_ser_p2_ph_x0_dcoefl1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 3; ++a) I[(3 + a)*9 + j] = G1[a]*dv11*jacob_vy_inv; 
    } 
  } 
  return 6; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_x_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[18]; 
  double I[54]; 
  hamil_phase_alpha_quad_x_2x2v_ser_p2_shared(0, 1, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 3; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 6; ++t) alpha += O[t*3 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
