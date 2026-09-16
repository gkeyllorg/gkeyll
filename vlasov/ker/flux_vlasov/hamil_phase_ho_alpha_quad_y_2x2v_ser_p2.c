#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p2.h> 
GKYL_CU_DH double hamil_phase_ho_alpha_quad_y_2x2v_ser_p2_node(int i, int m, int hamil_pt_edge, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[4]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[16]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[24]; 
  if (hamil_pt_edge == -1) { 
    double G0[3]; 
    for (int a = 0; a < 3; ++a) G0[a] = 0.0; 
    for (int k = 0; k < 48; ++k) { 
      G0[vst_2x2v_ser_p2_ho_ph_x1_cmap[k]] += vst_2x2v_ser_p2_ho_ph_x1_Vd0[m*13 + vst_2x2v_ser_p2_ho_ph_x1_vrd0map[k]]*(vst_2x2v_ser_p2_ho_ph_x1_dcoefr0[k]*hamil[k]); 
    } 
    double G1[3]; 
    for (int a = 0; a < 3; ++a) G1[a] = 0.0; 
    for (int k = 0; k < 48; ++k) { 
      G1[vst_2x2v_ser_p2_ho_ph_x1_cmap[k]] += vst_2x2v_ser_p2_ho_ph_x1_Vd1[m*13 + vst_2x2v_ser_p2_ho_ph_x1_vrd1map[k]]*(vst_2x2v_ser_p2_ho_ph_x1_dcoefr1[k]*hamil[k]); 
    } 
    double P0 = 0.0; 
    for (int a = 0; a < 8; ++a) P0 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_0[a]; 
    double dH0 = 0.0; 
    for (int a = 0; a < 3; ++a) dH0 += vst_2x2v_ser_p2_ho_ph_x1_Cm[i*3 + a]*G0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 8; ++a) P1 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_r[i*8 + a]*poisson_tensor_conf_1[a]; 
    double dH1 = 0.0; 
    for (int a = 0; a < 3; ++a) dH1 += vst_2x2v_ser_p2_ho_ph_x1_Cm[i*3 + a]*G1[a]; 
    return P0*dH0*dv10*jacob_vx_inv + P1*dH1*dv11*jacob_vy_inv; 
  } 
  else if (hamil_pt_edge == 1) { 
    double G0[3]; 
    for (int a = 0; a < 3; ++a) G0[a] = 0.0; 
    for (int k = 0; k < 48; ++k) { 
      G0[vst_2x2v_ser_p2_ho_ph_x1_cmap[k]] += vst_2x2v_ser_p2_ho_ph_x1_Vd0[m*13 + vst_2x2v_ser_p2_ho_ph_x1_vld0map[k]]*(vst_2x2v_ser_p2_ho_ph_x1_dcoefl0[k]*hamil[k]); 
    } 
    double G1[3]; 
    for (int a = 0; a < 3; ++a) G1[a] = 0.0; 
    for (int k = 0; k < 48; ++k) { 
      G1[vst_2x2v_ser_p2_ho_ph_x1_cmap[k]] += vst_2x2v_ser_p2_ho_ph_x1_Vd1[m*13 + vst_2x2v_ser_p2_ho_ph_x1_vld1map[k]]*(vst_2x2v_ser_p2_ho_ph_x1_dcoefl1[k]*hamil[k]); 
    } 
    double P0 = 0.0; 
    for (int a = 0; a < 8; ++a) P0 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_0[a]; 
    double dH0 = 0.0; 
    for (int a = 0; a < 3; ++a) dH0 += vst_2x2v_ser_p2_ho_ph_x1_Cm[i*3 + a]*G0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 8; ++a) P1 += vst_2x2v_ser_p2_ho_confsurf_x1_ev_l[i*8 + a]*poisson_tensor_conf_1[a]; 
    double dH1 = 0.0; 
    for (int a = 0; a < 3; ++a) dH1 += vst_2x2v_ser_p2_ho_ph_x1_Cm[i*3 + a]*G1[a]; 
    return P0*dH0*dv10*jacob_vx_inv + P1*dH1*dv11*jacob_vy_inv; 
  } 
  return 0.0; 
} 

GKYL_CU_DH void hamil_phase_ho_alpha_quad_y_2x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int m = 0; m < 16; ++m) alpha_quad[i*16 + m] += hamil_phase_ho_alpha_quad_y_2x2v_ser_p2_node(i, m, hamil_pt_edge, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
