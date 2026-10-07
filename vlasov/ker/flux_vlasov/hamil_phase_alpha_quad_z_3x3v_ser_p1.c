#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_ser_p1.h> 
GKYL_CU_DH int hamil_phase_alpha_quad_z_3x3v_ser_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 12; 
  const double dv10 = 2.0/dxv[3]; 
  const double dv11 = 2.0/dxv[4]; 
  const double dv12 = 2.0/dxv[5]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[3]; 
  const double jacob_vz_inv = 1.0/jacob_vel_surf[6]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[48]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[56]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[64]; 
  for (int i = tid; i < 4; i += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_3x3v_ser_p1_confsurf_x2_ev_r[i*8 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 4; ++a) O[(0 + a)*4 + i] = P0*vst_3x3v_ser_p1_ph_x2_Cm[i*4 + a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_3x3v_ser_p1_confsurf_x2_ev_r[i*8 + a]*poisson_tensor_conf_1[a]; 
      for (int a = 0; a < 4; ++a) O[(4 + a)*4 + i] = P1*vst_3x3v_ser_p1_ph_x2_Cm[i*4 + a]; 
      double P2 = 0.0; 
      for (int a = 0; a < 8; ++a) P2 += vst_3x3v_ser_p1_confsurf_x2_ev_r[i*8 + a]*poisson_tensor_conf_2[a]; 
      for (int a = 0; a < 4; ++a) O[(8 + a)*4 + i] = P2*vst_3x3v_ser_p1_ph_x2_Cm[i*4 + a]; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 8; ++a) P0 += vst_3x3v_ser_p1_confsurf_x2_ev_l[i*8 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 4; ++a) O[(0 + a)*4 + i] = P0*vst_3x3v_ser_p1_ph_x2_Cm[i*4 + a]; 
      double P1 = 0.0; 
      for (int a = 0; a < 8; ++a) P1 += vst_3x3v_ser_p1_confsurf_x2_ev_l[i*8 + a]*poisson_tensor_conf_1[a]; 
      for (int a = 0; a < 4; ++a) O[(4 + a)*4 + i] = P1*vst_3x3v_ser_p1_ph_x2_Cm[i*4 + a]; 
      double P2 = 0.0; 
      for (int a = 0; a < 8; ++a) P2 += vst_3x3v_ser_p1_confsurf_x2_ev_l[i*8 + a]*poisson_tensor_conf_2[a]; 
      for (int a = 0; a < 4; ++a) O[(8 + a)*4 + i] = P2*vst_3x3v_ser_p1_ph_x2_Cm[i*4 + a]; 
    } 
  } 
  for (int j = tid; j < 8; j += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double G0[4]; 
      for (int a = 0; a < 4; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 64; ++k) { 
        G0[vst_3x3v_ser_p1_ph_x2_cmap[k]] += vst_3x3v_ser_p1_ph_x2_Vd0[j*9 + vst_3x3v_ser_p1_ph_x2_vrd0map[k]]*(vst_3x3v_ser_p1_ph_x2_dcoefr0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(0 + a)*8 + j] = G0[a]*dv10*jacob_vx_inv; 
      double G1[4]; 
      for (int a = 0; a < 4; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 64; ++k) { 
        G1[vst_3x3v_ser_p1_ph_x2_cmap[k]] += vst_3x3v_ser_p1_ph_x2_Vd1[j*9 + vst_3x3v_ser_p1_ph_x2_vrd1map[k]]*(vst_3x3v_ser_p1_ph_x2_dcoefr1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(4 + a)*8 + j] = G1[a]*dv11*jacob_vy_inv; 
      double G2[4]; 
      for (int a = 0; a < 4; ++a) G2[a] = 0.0; 
      for (int k = 0; k < 64; ++k) { 
        G2[vst_3x3v_ser_p1_ph_x2_cmap[k]] += vst_3x3v_ser_p1_ph_x2_Vd2[j*9 + vst_3x3v_ser_p1_ph_x2_vrd2map[k]]*(vst_3x3v_ser_p1_ph_x2_dcoefr2[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(8 + a)*8 + j] = G2[a]*dv12*jacob_vz_inv; 
    } 
    else if (hamil_pt_edge == 1) { 
      double G0[4]; 
      for (int a = 0; a < 4; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 64; ++k) { 
        G0[vst_3x3v_ser_p1_ph_x2_cmap[k]] += vst_3x3v_ser_p1_ph_x2_Vd0[j*9 + vst_3x3v_ser_p1_ph_x2_vld0map[k]]*(vst_3x3v_ser_p1_ph_x2_dcoefl0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(0 + a)*8 + j] = G0[a]*dv10*jacob_vx_inv; 
      double G1[4]; 
      for (int a = 0; a < 4; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 64; ++k) { 
        G1[vst_3x3v_ser_p1_ph_x2_cmap[k]] += vst_3x3v_ser_p1_ph_x2_Vd1[j*9 + vst_3x3v_ser_p1_ph_x2_vld1map[k]]*(vst_3x3v_ser_p1_ph_x2_dcoefl1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(4 + a)*8 + j] = G1[a]*dv11*jacob_vy_inv; 
      double G2[4]; 
      for (int a = 0; a < 4; ++a) G2[a] = 0.0; 
      for (int k = 0; k < 64; ++k) { 
        G2[vst_3x3v_ser_p1_ph_x2_cmap[k]] += vst_3x3v_ser_p1_ph_x2_Vd2[j*9 + vst_3x3v_ser_p1_ph_x2_vld2map[k]]*(vst_3x3v_ser_p1_ph_x2_dcoefl2[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(8 + a)*8 + j] = G2[a]*dv12*jacob_vz_inv; 
    } 
  } 
  return 12; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_z_3x3v_ser_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[48]; 
  double I[96]; 
  hamil_phase_alpha_quad_z_3x3v_ser_p1_shared(0, 1, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 8; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 12; ++t) alpha += O[t*4 + i]*I[t*8 + j]; 
      alpha_quad[i*8 + j] += alpha; 
    } 
  } 
} 
