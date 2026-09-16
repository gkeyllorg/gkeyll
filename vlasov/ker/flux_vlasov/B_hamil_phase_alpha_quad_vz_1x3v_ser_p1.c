#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p1.h> 
GKYL_CU_DH double B_hamil_phase_alpha_quad_vz_1x3v_ser_p1_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvx = 0.0; 
  { 
    double G[2]; 
    for (int a = 0; a < 2; ++a) G[a] = 0.0; 
    for (int b = 0; b < 16; ++b) { 
      G[vst_1x3v_ser_p1_ph_v2_cmap[b]] += vst_1x3v_ser_p1_ph_v2_Vd0[j*5 + vst_1x3v_ser_p1_ph_v2_vrd0map[b]]*(vst_1x3v_ser_p1_ph_v2_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 2; ++a) dH_dvx += vst_1x3v_ser_p1_ph_v2_Cm[i*2 + a]*G[a]; 
  } 
  double dH_dvy = 0.0; 
  { 
    double G[2]; 
    for (int a = 0; a < 2; ++a) G[a] = 0.0; 
    for (int b = 0; b < 16; ++b) { 
      G[vst_1x3v_ser_p1_ph_v2_cmap[b]] += vst_1x3v_ser_p1_ph_v2_Vd1[j*5 + vst_1x3v_ser_p1_ph_v2_vrd1map[b]]*(vst_1x3v_ser_p1_ph_v2_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 2; ++a) dH_dvy += vst_1x3v_ser_p1_ph_v2_Cm[i*2 + a]*G[a]; 
  } 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Bx = &qmem[6]; 
  const double *By = &qmem[8]; 
  double Bx_quad = 0.0; 
  double By_quad = 0.0; 
  for (int a = 0; a < 2; ++a) { 
    Bx_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*Bx[a]; 
    By_quad += vst_1x3v_ser_p1_conf_ev[i*2 + a]*By[a]; 
  } 
  return dv10*dH_dvx*By_quad/jacob_vel_surf[0] - dv11*dH_dvy*Bx_quad/jacob_vel_surf[3]; 
} 

GKYL_CU_DH void B_hamil_phase_alpha_quad_vz_1x3v_ser_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += B_hamil_phase_alpha_quad_vz_1x3v_ser_p1_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
