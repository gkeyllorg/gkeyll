#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p1.h> 
GKYL_CU_DH int hamil_phase_alpha_quad_x_1x3v_ser_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 3; 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double dv12 = 2.0/dxv[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[3]; 
  const double jacob_vz_inv = 1.0/jacob_vel_surf[6]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[2]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[4]; 
  for (int i = tid; i < 1; i += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 2; ++a) P0 += vst_1x3v_ser_p1_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 1; ++a) O[(0 + a)*1 + i] = P0*vst_1x3v_ser_p1_ph_x0_Cm[i*1 + a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 2; ++a) P1 += vst_1x3v_ser_p1_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_1[a]; 
      for (int a = 0; a < 1; ++a) O[(1 + a)*1 + i] = P1*vst_1x3v_ser_p1_ph_x0_Cm[i*1 + a]; 
      double P2 = 0.0; 
      for (int a = 0; a < 2; ++a) P2 += vst_1x3v_ser_p1_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_2[a]; 
      for (int a = 0; a < 1; ++a) O[(2 + a)*1 + i] = P2*vst_1x3v_ser_p1_ph_x0_Cm[i*1 + a]; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 2; ++a) P0 += vst_1x3v_ser_p1_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 1; ++a) O[(0 + a)*1 + i] = P0*vst_1x3v_ser_p1_ph_x0_Cm[i*1 + a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 2; ++a) P1 += vst_1x3v_ser_p1_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_1[a]; 
      for (int a = 0; a < 1; ++a) O[(1 + a)*1 + i] = P1*vst_1x3v_ser_p1_ph_x0_Cm[i*1 + a]; 
      double P2 = 0.0; 
      for (int a = 0; a < 2; ++a) P2 += vst_1x3v_ser_p1_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_2[a]; 
      for (int a = 0; a < 1; ++a) O[(2 + a)*1 + i] = P2*vst_1x3v_ser_p1_ph_x0_Cm[i*1 + a]; 
    } 
  } 
  for (int j = tid; j < 8; j += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double G0[1]; 
      for (int a = 0; a < 1; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 16; ++k) { 
        G0[vst_1x3v_ser_p1_ph_x0_cmap[k]] += vst_1x3v_ser_p1_ph_x0_Vd0[j*9 + vst_1x3v_ser_p1_ph_x0_vrd0map[k]]*(vst_1x3v_ser_p1_ph_x0_dcoefr0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 1; ++a) I[(0 + a)*8 + j] = G0[a]*dv10*jacob_vx_inv; 
      double G1[1]; 
      for (int a = 0; a < 1; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 16; ++k) { 
        G1[vst_1x3v_ser_p1_ph_x0_cmap[k]] += vst_1x3v_ser_p1_ph_x0_Vd1[j*9 + vst_1x3v_ser_p1_ph_x0_vrd1map[k]]*(vst_1x3v_ser_p1_ph_x0_dcoefr1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 1; ++a) I[(1 + a)*8 + j] = G1[a]*dv11*jacob_vy_inv; 
      double G2[1]; 
      for (int a = 0; a < 1; ++a) G2[a] = 0.0; 
      for (int k = 0; k < 16; ++k) { 
        G2[vst_1x3v_ser_p1_ph_x0_cmap[k]] += vst_1x3v_ser_p1_ph_x0_Vd2[j*9 + vst_1x3v_ser_p1_ph_x0_vrd2map[k]]*(vst_1x3v_ser_p1_ph_x0_dcoefr2[k]*hamil[k]); 
      } 
      for (int a = 0; a < 1; ++a) I[(2 + a)*8 + j] = G2[a]*dv12*jacob_vz_inv; 
    } 
    else if (hamil_pt_edge == 1) { 
      double G0[1]; 
      for (int a = 0; a < 1; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 16; ++k) { 
        G0[vst_1x3v_ser_p1_ph_x0_cmap[k]] += vst_1x3v_ser_p1_ph_x0_Vd0[j*9 + vst_1x3v_ser_p1_ph_x0_vld0map[k]]*(vst_1x3v_ser_p1_ph_x0_dcoefl0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 1; ++a) I[(0 + a)*8 + j] = G0[a]*dv10*jacob_vx_inv; 
      double G1[1]; 
      for (int a = 0; a < 1; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 16; ++k) { 
        G1[vst_1x3v_ser_p1_ph_x0_cmap[k]] += vst_1x3v_ser_p1_ph_x0_Vd1[j*9 + vst_1x3v_ser_p1_ph_x0_vld1map[k]]*(vst_1x3v_ser_p1_ph_x0_dcoefl1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 1; ++a) I[(1 + a)*8 + j] = G1[a]*dv11*jacob_vy_inv; 
      double G2[1]; 
      for (int a = 0; a < 1; ++a) G2[a] = 0.0; 
      for (int k = 0; k < 16; ++k) { 
        G2[vst_1x3v_ser_p1_ph_x0_cmap[k]] += vst_1x3v_ser_p1_ph_x0_Vd2[j*9 + vst_1x3v_ser_p1_ph_x0_vld2map[k]]*(vst_1x3v_ser_p1_ph_x0_dcoefl2[k]*hamil[k]); 
      } 
      for (int a = 0; a < 1; ++a) I[(2 + a)*8 + j] = G2[a]*dv12*jacob_vz_inv; 
    } 
  } 
  return 3; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_x_1x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[3]; 
  double I[24]; 
  hamil_phase_alpha_quad_x_1x3v_ser_p1_shared(0, 1, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 1; ++i) { 
    for (int j = 0; j < 8; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 3; ++t) alpha += O[t*1 + i]*I[t*8 + j]; 
      alpha_quad[i*8 + j] += alpha; 
    } 
  } 
} 
