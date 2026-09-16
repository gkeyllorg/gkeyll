#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double hamil_phase_alpha_quad_x_3x3v_tensor_p1_node(int i, int m, int hamil_pt_edge, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  const double dv10 = 2.0/dxv[3]; 
  const double dv11 = 2.0/dxv[4]; 
  const double dv12 = 2.0/dxv[5]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[8]; 
  const double *poisson_tensor_conf_2 = &poisson_tensor_conf[16]; 
  if (hamil_pt_edge == -1) { 
    double G0[4]; 
    for (int a = 0; a < 4; ++a) G0[a] = 0.0; 
    for (int k = 0; k < 216; ++k) { 
      G0[vst_3x3v_tensor_p1_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ph_x0_Vd0[m*37 + vst_3x3v_tensor_p1_ph_x0_vrd0map[k]]*(vst_3x3v_tensor_p1_ph_x0_dcoefr0[k]*hamil[k]); 
    } 
    double G1[4]; 
    for (int a = 0; a < 4; ++a) G1[a] = 0.0; 
    for (int k = 0; k < 216; ++k) { 
      G1[vst_3x3v_tensor_p1_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ph_x0_Vd1[m*37 + vst_3x3v_tensor_p1_ph_x0_vrd1map[k]]*(vst_3x3v_tensor_p1_ph_x0_dcoefr1[k]*hamil[k]); 
    } 
    double G2[4]; 
    for (int a = 0; a < 4; ++a) G2[a] = 0.0; 
    for (int k = 0; k < 216; ++k) { 
      G2[vst_3x3v_tensor_p1_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ph_x0_Vd2[m*37 + vst_3x3v_tensor_p1_ph_x0_vrd2map[k]]*(vst_3x3v_tensor_p1_ph_x0_dcoefr2[k]*hamil[k]); 
    } 
    double P0 = 0.0; 
    for (int a = 0; a < 8; ++a) P0 += vst_3x3v_tensor_p1_confsurf_x0_ev_r[i*8 + a]*poisson_tensor_conf_0[a]; 
    double dH0 = 0.0; 
    for (int a = 0; a < 4; ++a) dH0 += vst_3x3v_tensor_p1_ph_x0_Cm[i*4 + a]*G0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 8; ++a) P1 += vst_3x3v_tensor_p1_confsurf_x0_ev_r[i*8 + a]*poisson_tensor_conf_1[a]; 
    double dH1 = 0.0; 
    for (int a = 0; a < 4; ++a) dH1 += vst_3x3v_tensor_p1_ph_x0_Cm[i*4 + a]*G1[a]; 
    double P2 = 0.0; 
    for (int a = 0; a < 8; ++a) P2 += vst_3x3v_tensor_p1_confsurf_x0_ev_r[i*8 + a]*poisson_tensor_conf_2[a]; 
    double dH2 = 0.0; 
    for (int a = 0; a < 4; ++a) dH2 += vst_3x3v_tensor_p1_ph_x0_Cm[i*4 + a]*G2[a]; 
    return P0*dH0*dv10*(1.0/jacob_vel_surf[0 + m/9]) + P1*dH1*dv11*(1.0/jacob_vel_surf[3 + m/3%3]) + P2*dH2*dv12*(1.0/jacob_vel_surf[6 + m%3]); 
  } 
  else if (hamil_pt_edge == 1) { 
    double G0[4]; 
    for (int a = 0; a < 4; ++a) G0[a] = 0.0; 
    for (int k = 0; k < 216; ++k) { 
      G0[vst_3x3v_tensor_p1_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ph_x0_Vd0[m*37 + vst_3x3v_tensor_p1_ph_x0_vld0map[k]]*(vst_3x3v_tensor_p1_ph_x0_dcoefl0[k]*hamil[k]); 
    } 
    double G1[4]; 
    for (int a = 0; a < 4; ++a) G1[a] = 0.0; 
    for (int k = 0; k < 216; ++k) { 
      G1[vst_3x3v_tensor_p1_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ph_x0_Vd1[m*37 + vst_3x3v_tensor_p1_ph_x0_vld1map[k]]*(vst_3x3v_tensor_p1_ph_x0_dcoefl1[k]*hamil[k]); 
    } 
    double G2[4]; 
    for (int a = 0; a < 4; ++a) G2[a] = 0.0; 
    for (int k = 0; k < 216; ++k) { 
      G2[vst_3x3v_tensor_p1_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ph_x0_Vd2[m*37 + vst_3x3v_tensor_p1_ph_x0_vld2map[k]]*(vst_3x3v_tensor_p1_ph_x0_dcoefl2[k]*hamil[k]); 
    } 
    double P0 = 0.0; 
    for (int a = 0; a < 8; ++a) P0 += vst_3x3v_tensor_p1_confsurf_x0_ev_l[i*8 + a]*poisson_tensor_conf_0[a]; 
    double dH0 = 0.0; 
    for (int a = 0; a < 4; ++a) dH0 += vst_3x3v_tensor_p1_ph_x0_Cm[i*4 + a]*G0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 8; ++a) P1 += vst_3x3v_tensor_p1_confsurf_x0_ev_l[i*8 + a]*poisson_tensor_conf_1[a]; 
    double dH1 = 0.0; 
    for (int a = 0; a < 4; ++a) dH1 += vst_3x3v_tensor_p1_ph_x0_Cm[i*4 + a]*G1[a]; 
    double P2 = 0.0; 
    for (int a = 0; a < 8; ++a) P2 += vst_3x3v_tensor_p1_confsurf_x0_ev_l[i*8 + a]*poisson_tensor_conf_2[a]; 
    double dH2 = 0.0; 
    for (int a = 0; a < 4; ++a) dH2 += vst_3x3v_tensor_p1_ph_x0_Cm[i*4 + a]*G2[a]; 
    return P0*dH0*dv10*(1.0/jacob_vel_surf[0 + m/9]) + P1*dH1*dv11*(1.0/jacob_vel_surf[3 + m/3%3]) + P2*dH2*dv12*(1.0/jacob_vel_surf[6 + m%3]); 
  } 
  return 0.0; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int m = 0; m < 27; ++m) alpha_quad[i*27 + m] += hamil_phase_alpha_quad_x_3x3v_tensor_p1_node(i, m, hamil_pt_edge, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
