#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double B_ho_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvx = 0.0; 
  { 
    double G[8]; 
    for (int a = 0; a < 8; ++a) G[a] = 0.0; 
    for (int b = 0; b < 216; ++b) { 
      G[vst_3x3v_tensor_p1_ho_ph_v1_cmap[b]] += vst_3x3v_tensor_p1_ho_ph_v1_Vd0[j*19 + vst_3x3v_tensor_p1_ho_ph_v1_vrd0map[b]]*(vst_3x3v_tensor_p1_ho_ph_v1_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 8; ++a) dH_dvx += vst_3x3v_tensor_p1_ho_ph_v1_Cm[i*8 + a]*G[a]; 
  } 
  double dH_dvz = 0.0; 
  { 
    double G[8]; 
    for (int a = 0; a < 8; ++a) G[a] = 0.0; 
    for (int b = 0; b < 216; ++b) { 
      G[vst_3x3v_tensor_p1_ho_ph_v1_cmap[b]] += vst_3x3v_tensor_p1_ho_ph_v1_Vd2[j*19 + vst_3x3v_tensor_p1_ho_ph_v1_vrd2map[b]]*(vst_3x3v_tensor_p1_ho_ph_v1_dcoefr2[b]*hamil[b]); 
    } 
    for (int a = 0; a < 8; ++a) dH_dvz += vst_3x3v_tensor_p1_ho_ph_v1_Cm[i*8 + a]*G[a]; 
  } 
  double dv10 = 2.0/dxv[3]; 
  double dv12 = 2.0/dxv[5]; 
  const double *Bx = &qmem[24]; 
  const double *Bz = &qmem[40]; 
  double Bx_quad = 0.0; 
  double Bz_quad = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    Bx_quad += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*Bx[a]; 
    Bz_quad += vst_3x3v_tensor_p1_ho_conf_ev[i*8 + a]*Bz[a]; 
  } 
  return dv12*dH_dvz*Bx_quad/jacob_vel_surf[8 + j%4] - dv10*dH_dvx*Bz_quad/jacob_vel_surf[0 + j/4]; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_3x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += B_ho_hamil_phase_alpha_quad_vy_3x3v_tensor_p1_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
