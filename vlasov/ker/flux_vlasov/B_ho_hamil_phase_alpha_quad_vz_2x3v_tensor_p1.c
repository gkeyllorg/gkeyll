#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH double B_ho_hamil_phase_alpha_quad_vz_2x3v_tensor_p1_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvx = 0.0; 
  { 
    double G[4]; 
    for (int a = 0; a < 4; ++a) G[a] = 0.0; 
    for (int b = 0; b < 108; ++b) { 
      G[vst_2x3v_tensor_p1_ho_ph_v2_cmap[b]] += vst_2x3v_tensor_p1_ho_ph_v2_Vd0[j*19 + vst_2x3v_tensor_p1_ho_ph_v2_vrd0map[b]]*(vst_2x3v_tensor_p1_ho_ph_v2_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 4; ++a) dH_dvx += vst_2x3v_tensor_p1_ho_ph_v2_Cm[i*4 + a]*G[a]; 
  } 
  double dH_dvy = 0.0; 
  { 
    double G[4]; 
    for (int a = 0; a < 4; ++a) G[a] = 0.0; 
    for (int b = 0; b < 108; ++b) { 
      G[vst_2x3v_tensor_p1_ho_ph_v2_cmap[b]] += vst_2x3v_tensor_p1_ho_ph_v2_Vd1[j*19 + vst_2x3v_tensor_p1_ho_ph_v2_vrd1map[b]]*(vst_2x3v_tensor_p1_ho_ph_v2_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 4; ++a) dH_dvy += vst_2x3v_tensor_p1_ho_ph_v2_Cm[i*4 + a]*G[a]; 
  } 
  double dv10 = 2.0/dxv[2]; 
  double dv11 = 2.0/dxv[3]; 
  const double *Bx = &qmem[12]; 
  const double *By = &qmem[16]; 
  double Bx_quad = 0.0; 
  double By_quad = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    Bx_quad += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*Bx[a]; 
    By_quad += vst_2x3v_tensor_p1_ho_conf_ev[i*4 + a]*By[a]; 
  } 
  return dv10*dH_dvx*By_quad/jacob_vel_surf[0 + j/4] - dv11*dH_dvy*Bx_quad/jacob_vel_surf[4 + j%4]; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vz_2x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += B_ho_hamil_phase_alpha_quad_vz_2x3v_tensor_p1_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
