#include <gkyl_vlasov_triad_nodal_kernels.h> 
#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH int triad_nodal_hamil_phase_alpha_quad_y_2x2v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 4; 
  O += off*2; 
  I += off*9; 
  const double dv10 = 2.0/dxv[2]; 
  const double dv11 = 2.0/dxv[3]; 
  const double *pt_face = &poisson_tensor_conf[(2 + (hamil_pt_edge == -1 ? 0 : 1))*4]; 
  if (tid < 2) { 
    const int i = tid; 
    for (int a = 0; a < 2; ++a) O[(0 + a)*2 + i] = pt_face[i*2 + 0]*vst_2x2v_tensor_p1_ph_x1_Cm[i*2 + a]; 
    for (int a = 0; a < 2; ++a) O[(2 + a)*2 + i] = pt_face[i*2 + 1]*vst_2x2v_tensor_p1_ph_x1_Cm[i*2 + a]; 
  } 
  if (tid < 9) { 
    const int j = tid; 
    if (hamil_pt_edge == -1) { 
      double G0[2]; 
      for (int a = 0; a < 2; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 36; ++k) { 
        G0[vst_2x2v_tensor_p1_ph_x1_cmap[k]] += vst_2x2v_tensor_p1_ph_x1_Vd0[j*13 + vst_2x2v_tensor_p1_ph_x1_vrd0map[k]]*(vst_2x2v_tensor_p1_ph_x1_dcoefr0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(0 + a)*9 + j] = G0[a]*dv10*(1.0/jacob_vel_surf[0 + j/3]); 
      double G1[2]; 
      for (int a = 0; a < 2; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 36; ++k) { 
        G1[vst_2x2v_tensor_p1_ph_x1_cmap[k]] += vst_2x2v_tensor_p1_ph_x1_Vd1[j*13 + vst_2x2v_tensor_p1_ph_x1_vrd1map[k]]*(vst_2x2v_tensor_p1_ph_x1_dcoefr1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(2 + a)*9 + j] = G1[a]*dv11*(1.0/jacob_vel_surf[3 + j%3]); 
    } 
    else if (hamil_pt_edge == 1) { 
      double G0[2]; 
      for (int a = 0; a < 2; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 36; ++k) { 
        G0[vst_2x2v_tensor_p1_ph_x1_cmap[k]] += vst_2x2v_tensor_p1_ph_x1_Vd0[j*13 + vst_2x2v_tensor_p1_ph_x1_vld0map[k]]*(vst_2x2v_tensor_p1_ph_x1_dcoefl0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(0 + a)*9 + j] = G0[a]*dv10*(1.0/jacob_vel_surf[0 + j/3]); 
      double G1[2]; 
      for (int a = 0; a < 2; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 36; ++k) { 
        G1[vst_2x2v_tensor_p1_ph_x1_cmap[k]] += vst_2x2v_tensor_p1_ph_x1_Vd1[j*13 + vst_2x2v_tensor_p1_ph_x1_vld1map[k]]*(vst_2x2v_tensor_p1_ph_x1_dcoefl1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 2; ++a) I[(2 + a)*9 + j] = G1[a]*dv11*(1.0/jacob_vel_surf[3 + j%3]); 
    } 
  } 
  return 4; 
} 

GKYL_CU_DH void triad_nodal_hamil_phase_alpha_quad_y_2x2v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[8]; 
  double I[36]; 
  for (int tid = 0; tid < 9; ++tid) triad_nodal_hamil_phase_alpha_quad_y_2x2v_tensor_p1_shared(tid, 0, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 4; ++t) alpha += O[t*2 + i]*I[t*9 + j]; 
      alpha_quad[i*9 + j] += alpha; 
    } 
  } 
} 
