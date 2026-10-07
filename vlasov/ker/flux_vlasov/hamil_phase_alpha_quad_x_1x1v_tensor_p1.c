#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x1v_tensor_p1.h> 
GKYL_CU_DH int hamil_phase_alpha_quad_x_1x1v_tensor_p1_shared(int tid, int nthreads, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 1; 
  const double dv10 = 2.0/dxv[1]; 
  const double *poisson_tensor_conf_0 = &poisson_tensor_conf[0]; 
  for (int i = tid; i < 1; i += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 2; ++a) P0 += vst_1x1v_tensor_p1_confsurf_x0_ev_r[i*2 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 1; ++a) O[(0 + a)*1 + i] = P0*vst_1x1v_tensor_p1_ph_x0_Cm[i*1 + a]; 
    } 
    else if (hamil_pt_edge == 1) { 
      double P0 = 0.0; 
      for (int a = 0; a < 2; ++a) P0 += vst_1x1v_tensor_p1_confsurf_x0_ev_l[i*2 + a]*poisson_tensor_conf_0[a]; 
      for (int a = 0; a < 1; ++a) O[(0 + a)*1 + i] = P0*vst_1x1v_tensor_p1_ph_x0_Cm[i*1 + a]; 
    } 
  } 
  for (int j = tid; j < 3; j += nthreads) { 
    if (hamil_pt_edge == -1) { 
      double G0[1]; 
      for (int a = 0; a < 1; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 6; ++k) { 
        G0[vst_1x1v_tensor_p1_ph_x0_cmap[k]] += vst_1x1v_tensor_p1_ph_x0_Vd0[j*5 + vst_1x1v_tensor_p1_ph_x0_vrd0map[k]]*(vst_1x1v_tensor_p1_ph_x0_dcoefr0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 1; ++a) I[(0 + a)*3 + j] = G0[a]*dv10*(1.0/jacob_vel_surf[0 + j]); 
    } 
    else if (hamil_pt_edge == 1) { 
      double G0[1]; 
      for (int a = 0; a < 1; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 6; ++k) { 
        G0[vst_1x1v_tensor_p1_ph_x0_cmap[k]] += vst_1x1v_tensor_p1_ph_x0_Vd0[j*5 + vst_1x1v_tensor_p1_ph_x0_vld0map[k]]*(vst_1x1v_tensor_p1_ph_x0_dcoefl0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 1; ++a) I[(0 + a)*3 + j] = G0[a]*dv10*(1.0/jacob_vel_surf[0 + j]); 
    } 
  } 
  return 1; 
} 

GKYL_CU_DH void hamil_phase_alpha_quad_x_1x1v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[1]; 
  double I[3]; 
  hamil_phase_alpha_quad_x_1x1v_tensor_p1_shared(0, 1, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 1; ++i) { 
    for (int j = 0; j < 3; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 1; ++t) alpha += O[t*1 + i]*I[t*3 + j]; 
      alpha_quad[i*3 + j] += alpha; 
    } 
  } 
} 
