#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x1v_tensor_p1.h> 
GKYL_CU_DH double hamil_phase_alpha_quad_x_1x1v_tensor_p1_node(int i, int m, int hamil_pt_edge, const double *w, const double *dxv,
  const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
  const double *poisson_tensor_conf, const double *hamil) 
{ 
  const double dv10 = 2.0/dxv[1]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  if (hamil_pt_edge == -1) { 
    double G0[1]; 
    for (int a = 0; a < 1; ++a) G0[a] = 0.0; 
    for (int k = 0; k < 6; ++k) { 
      G0[vst_1x1v_tensor_p1_ph_x0_cmap[k]] += vst_1x1v_tensor_p1_ph_x0_Vd0[m*5 + vst_1x1v_tensor_p1_ph_x0_vrd0map[k]]*(vst_1x1v_tensor_p1_ph_x0_dcoefr0[k]*hamil[k]); 
    } 
    double P0 = 0.0; 
    for (int a = 0; a < 2; ++a) P0 += vst_1x1v_tensor_p1_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_0[a]; 
    double dH0 = 0.0; 
    for (int a = 0; a < 1; ++a) dH0 += vst_1x1v_tensor_p1_ph_x0_Cm[i*1 + a]*G0[a]; 
    return P0*dH0*dv10*(1.0/jacob_vel_surf[0 + m]); 
  } 
  else if (hamil_pt_edge == 1) { 
    double G0[1]; 
    for (int a = 0; a < 1; ++a) G0[a] = 0.0; 
    for (int k = 0; k < 6; ++k) { 
      G0[vst_1x1v_tensor_p1_ph_x0_cmap[k]] += vst_1x1v_tensor_p1_ph_x0_Vd0[m*5 + vst_1x1v_tensor_p1_ph_x0_vld0map[k]]*(vst_1x1v_tensor_p1_ph_x0_dcoefl0[k]*hamil[k]); 
    } 
    double P0 = 0.0; 
    for (int a = 0; a < 2; ++a) P0 += vst_1x1v_tensor_p1_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_0[a]; 
    double dH0 = 0.0; 
    for (int a = 0; a < 1; ++a) dH0 += vst_1x1v_tensor_p1_ph_x0_Cm[i*1 + a]*G0[a]; 
    return P0*dH0*dv10*(1.0/jacob_vel_surf[0 + m]); 
  } 
  return 0.0; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_x_1x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge,
    const double *vmap, const double *jacob_pos, const double *jacob_vel_surf,
    const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 1; ++i) { 
    for (int m = 0; m < 3; ++m) alpha_quad[i*3 + m] += hamil_phase_alpha_quad_x_1x1v_tensor_p1_node(i, m, hamil_pt_edge, w, dxv, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil); 
  } 
} 
