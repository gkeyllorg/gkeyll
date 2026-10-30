#include <gkyl_vlasov_triad_nodal_kernels.h> 
#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH int triad_nodal_hamil_phase_ho_alpha_quad_x_3x3v_tensor_p1_shared(int tid, int off, const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil,
  double* GKYL_RESTRICT O, double* GKYL_RESTRICT I) 
{ 
  if (O == NULL) return 12; 
  O += off*4; 
  I += off*64; 
  const double dv10 = 2.0/dxv[3]; 
  const double dv11 = 2.0/dxv[4]; 
  const double dv12 = 2.0/dxv[5]; 
  const double *pt_face = &poisson_tensor_conf[(0 + (hamil_pt_edge == -1 ? 0 : 1))*12]; 
  if (tid < 4) { 
    const int i = tid; 
    for (int a = 0; a < 4; ++a) O[(0 + a)*4 + i] = pt_face[i*3 + 0]*vst_3x3v_tensor_p1_ho_ph_x0_Cm[i*4 + a]; 
    for (int a = 0; a < 4; ++a) O[(4 + a)*4 + i] = pt_face[i*3 + 1]*vst_3x3v_tensor_p1_ho_ph_x0_Cm[i*4 + a]; 
    for (int a = 0; a < 4; ++a) O[(8 + a)*4 + i] = pt_face[i*3 + 2]*vst_3x3v_tensor_p1_ho_ph_x0_Cm[i*4 + a]; 
  } 
  if (tid < 64) { 
    const int j = tid; 
    if (hamil_pt_edge == -1) { 
      double G0[4]; 
      for (int a = 0; a < 4; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 216; ++k) { 
        G0[vst_3x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ho_ph_x0_Vd0[j*37 + vst_3x3v_tensor_p1_ho_ph_x0_vrd0map[k]]*(vst_3x3v_tensor_p1_ho_ph_x0_dcoefr0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(0 + a)*64 + j] = G0[a]*dv10*(1.0/jacob_vel_surf[0 + j/16]); 
      double G1[4]; 
      for (int a = 0; a < 4; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 216; ++k) { 
        G1[vst_3x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ho_ph_x0_Vd1[j*37 + vst_3x3v_tensor_p1_ho_ph_x0_vrd1map[k]]*(vst_3x3v_tensor_p1_ho_ph_x0_dcoefr1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(4 + a)*64 + j] = G1[a]*dv11*(1.0/jacob_vel_surf[4 + j/4%4]); 
      double G2[4]; 
      for (int a = 0; a < 4; ++a) G2[a] = 0.0; 
      for (int k = 0; k < 216; ++k) { 
        G2[vst_3x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ho_ph_x0_Vd2[j*37 + vst_3x3v_tensor_p1_ho_ph_x0_vrd2map[k]]*(vst_3x3v_tensor_p1_ho_ph_x0_dcoefr2[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(8 + a)*64 + j] = G2[a]*dv12*(1.0/jacob_vel_surf[8 + j%4]); 
    } 
    else if (hamil_pt_edge == 1) { 
      double G0[4]; 
      for (int a = 0; a < 4; ++a) G0[a] = 0.0; 
      for (int k = 0; k < 216; ++k) { 
        G0[vst_3x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ho_ph_x0_Vd0[j*37 + vst_3x3v_tensor_p1_ho_ph_x0_vld0map[k]]*(vst_3x3v_tensor_p1_ho_ph_x0_dcoefl0[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(0 + a)*64 + j] = G0[a]*dv10*(1.0/jacob_vel_surf[0 + j/16]); 
      double G1[4]; 
      for (int a = 0; a < 4; ++a) G1[a] = 0.0; 
      for (int k = 0; k < 216; ++k) { 
        G1[vst_3x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ho_ph_x0_Vd1[j*37 + vst_3x3v_tensor_p1_ho_ph_x0_vld1map[k]]*(vst_3x3v_tensor_p1_ho_ph_x0_dcoefl1[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(4 + a)*64 + j] = G1[a]*dv11*(1.0/jacob_vel_surf[4 + j/4%4]); 
      double G2[4]; 
      for (int a = 0; a < 4; ++a) G2[a] = 0.0; 
      for (int k = 0; k < 216; ++k) { 
        G2[vst_3x3v_tensor_p1_ho_ph_x0_cmap[k]] += vst_3x3v_tensor_p1_ho_ph_x0_Vd2[j*37 + vst_3x3v_tensor_p1_ho_ph_x0_vld2map[k]]*(vst_3x3v_tensor_p1_ho_ph_x0_dcoefl2[k]*hamil[k]); 
      } 
      for (int a = 0; a < 4; ++a) I[(8 + a)*64 + j] = G2[a]*dv12*(1.0/jacob_vel_surf[8 + j%4]); 
    } 
  } 
  return 12; 
} 

GKYL_CU_DH void triad_nodal_hamil_phase_ho_alpha_quad_x_3x3v_tensor_p1(const double *w, const double *dxv, const int hamil_pt_edge, const double *vmap, const double *jacob_pos, const double *jacob_vel_surf, const double *poisson_tensor_conf, const double *hamil, double* GKYL_RESTRICT alpha_quad) 
{ 
  double O[48]; 
  double I[768]; 
  for (int tid = 0; tid < 64; ++tid) triad_nodal_hamil_phase_ho_alpha_quad_x_3x3v_tensor_p1_shared(tid, 0, w, dxv, hamil_pt_edge, vmap, jacob_pos, jacob_vel_surf, poisson_tensor_conf, hamil, O, I); 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 64; ++j) { 
      double alpha = 0.0; 
      for (int t = 0; t < 12; ++t) alpha += O[t*4 + i]*I[t*64 + j]; 
      alpha_quad[i*64 + j] += alpha; 
    } 
  } 
} 
