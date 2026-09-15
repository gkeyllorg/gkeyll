#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p2.h> 
GKYL_CU_DH double B_ho_hamil_phase_alpha_quad_vx_2x2v_ser_p2_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvy = 0.0; 
  { 
    double G[8]; 
    for (int a = 0; a < 8; ++a) G[a] = 0.0; 
    for (int b = 0; b < 48; ++b) { 
      G[vst_2x2v_ser_p2_ho_ph_v0_cmap[b]] += vst_2x2v_ser_p2_ho_ph_v0_Vd1[j*6 + vst_2x2v_ser_p2_ho_ph_v0_vrd1map[b]]*(vst_2x2v_ser_p2_ho_ph_v0_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 8; ++a) dH_dvy += vst_2x2v_ser_p2_ho_ph_v0_Cm[i*8 + a]*G[a]; 
  } 
  const double *Bz = &qmem[40]; 
  double Bz_quad = 0.0; 
  for (int a = 0; a < 8; ++a) Bz_quad += vst_2x2v_ser_p2_ho_conf_ev[i*8 + a]*Bz[a]; 
  return 2.0/(dxv[3]*jacob_vel_surf[4])*dH_dvy*Bz_quad; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_2x2v_ser_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += B_ho_hamil_phase_alpha_quad_vx_2x2v_ser_p2_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
