#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_ser_p2.h> 
GKYL_CU_DH double hamil_phase_ho_alpha_quad_x_1x2v_ser_p2_node(int i, int m, int hamil_pt_edge, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  const double dv10 = 2.0/dxv[1]; 
  const double dv11 = 2.0/dxv[2]; 
  const double jacob_vx_inv = 1.0/jacob_vel_surf[0]; 
  const double jacob_vy_inv = 1.0/jacob_vel_surf[4]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  const double *poisson_tensor_conf_1 = &poisson_tensor_conf[3]; 
  if (hamil_pt_edge == -1) { 
    double G0[1]; 
    for (int a = 0; a < 1; ++a) G0[a] = 0.0; 
    for (int k = 0; k < 20; ++k) { 
      G0[vst_1x2v_ser_p2_ho_ph_x0_cmap[k]] += vst_1x2v_ser_p2_ho_ph_x0_Vd0[m*13 + vst_1x2v_ser_p2_ho_ph_x0_vrd0map[k]]*(vst_1x2v_ser_p2_ho_ph_x0_dcoefr0[k]*hamil[k]); 
    } 
    double G1[1]; 
    for (int a = 0; a < 1; ++a) G1[a] = 0.0; 
    for (int k = 0; k < 20; ++k) { 
      G1[vst_1x2v_ser_p2_ho_ph_x0_cmap[k]] += vst_1x2v_ser_p2_ho_ph_x0_Vd1[m*13 + vst_1x2v_ser_p2_ho_ph_x0_vrd1map[k]]*(vst_1x2v_ser_p2_ho_ph_x0_dcoefr1[k]*hamil[k]); 
    } 
    double P0 = 0.0; 
    for (int a = 0; a < 3; ++a) P0 += vst_1x2v_ser_p2_ho_confsurf_x0_ev_r[i*3 + a]*poisson_tensor_conf_0[a]; 
    double dH0 = 0.0; 
    for (int a = 0; a < 1; ++a) dH0 += vst_1x2v_ser_p2_ho_ph_x0_Cm[i*1 + a]*G0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 3; ++a) P1 += vst_1x2v_ser_p2_ho_confsurf_x0_ev_r[i*3 + a]*poisson_tensor_conf_1[a]; 
    double dH1 = 0.0; 
    for (int a = 0; a < 1; ++a) dH1 += vst_1x2v_ser_p2_ho_ph_x0_Cm[i*1 + a]*G1[a]; 
    return P0*dH0*dv10*jacob_vx_inv + P1*dH1*dv11*jacob_vy_inv; 
  } 
  else if (hamil_pt_edge == 1) { 
    double G0[1]; 
    for (int a = 0; a < 1; ++a) G0[a] = 0.0; 
    for (int k = 0; k < 20; ++k) { 
      G0[vst_1x2v_ser_p2_ho_ph_x0_cmap[k]] += vst_1x2v_ser_p2_ho_ph_x0_Vd0[m*13 + vst_1x2v_ser_p2_ho_ph_x0_vld0map[k]]*(vst_1x2v_ser_p2_ho_ph_x0_dcoefl0[k]*hamil[k]); 
    } 
    double G1[1]; 
    for (int a = 0; a < 1; ++a) G1[a] = 0.0; 
    for (int k = 0; k < 20; ++k) { 
      G1[vst_1x2v_ser_p2_ho_ph_x0_cmap[k]] += vst_1x2v_ser_p2_ho_ph_x0_Vd1[m*13 + vst_1x2v_ser_p2_ho_ph_x0_vld1map[k]]*(vst_1x2v_ser_p2_ho_ph_x0_dcoefl1[k]*hamil[k]); 
    } 
    double P0 = 0.0; 
    for (int a = 0; a < 3; ++a) P0 += vst_1x2v_ser_p2_ho_confsurf_x0_ev_l[i*3 + a]*poisson_tensor_conf_0[a]; 
    double dH0 = 0.0; 
    for (int a = 0; a < 1; ++a) dH0 += vst_1x2v_ser_p2_ho_ph_x0_Cm[i*1 + a]*G0[a]; 
    double P1 = 0.0; 
    for (int a = 0; a < 3; ++a) P1 += vst_1x2v_ser_p2_ho_confsurf_x0_ev_l[i*3 + a]*poisson_tensor_conf_1[a]; 
    double dH1 = 0.0; 
    for (int a = 0; a < 1; ++a) dH1 += vst_1x2v_ser_p2_ho_ph_x0_Cm[i*1 + a]*G1[a]; 
    return P0*dH0*dv10*jacob_vx_inv + P1*dH1*dv11*jacob_vy_inv; 
  } 
  return 0.0; 
} 

GKYL_CU_DH void hamil_phase_ho_alpha_quad_x_1x2v_ser_p2(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 1; ++i) { 
    for (int m = 0; m < 16; ++m) alpha_quad[i*16 + m] += hamil_phase_ho_alpha_quad_x_1x2v_ser_p2_node(i, m, hamil_pt_edge, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
