#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double B_ho_hamil_phase_alpha_quad_vz_1x3v_ser_p2_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvx = 0.0; 
  { 
    double G[3]; 
    for (int a = 0; a < 3; ++a) G[a] = 0.0; 
    for (int b = 0; b < 48; ++b) { 
      G[vst_1x3v_ser_p2_ho_ph_v2_cmap[b]] += vst_1x3v_ser_p2_ho_ph_v2_Vd0[j*13 + vst_1x3v_ser_p2_ho_ph_v2_vrd0map[b]]*(vst_1x3v_ser_p2_ho_ph_v2_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 3; ++a) dH_dvx += vst_1x3v_ser_p2_ho_ph_v2_Cm[i*3 + a]*G[a]; 
  } 
  double dH_dvy = 0.0; 
  { 
    double G[3]; 
    for (int a = 0; a < 3; ++a) G[a] = 0.0; 
    for (int b = 0; b < 48; ++b) { 
      G[vst_1x3v_ser_p2_ho_ph_v2_cmap[b]] += vst_1x3v_ser_p2_ho_ph_v2_Vd1[j*13 + vst_1x3v_ser_p2_ho_ph_v2_vrd1map[b]]*(vst_1x3v_ser_p2_ho_ph_v2_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 3; ++a) dH_dvy += vst_1x3v_ser_p2_ho_ph_v2_Cm[i*3 + a]*G[a]; 
  } 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Bx = &qmem[9]; 
  const double *By = &qmem[12]; 
  double Bx_quad = 0.0; 
  double By_quad = 0.0; 
  for (int a = 0; a < 3; ++a) { 
    Bx_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*Bx[a]; 
    By_quad += vst_1x3v_ser_p2_ho_conf_ev[i*3 + a]*By[a]; 
  } 
  return dv10*dH_dvx*By_quad/jacob_vel_surf[0] - dv11*dH_dvy*Bx_quad/jacob_vel_surf[4]; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) alpha_quad[i*16 + j] += B_ho_hamil_phase_alpha_quad_vz_1x3v_ser_p2_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
