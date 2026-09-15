#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_ser_p2.h> 
GKYL_CU_DH double B_ho_hamil_phase_alpha_quad_vy_1x2v_ser_p2_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvx = 0.0; 
  { 
    double G[3]; 
    for (int a = 0; a < 3; ++a) G[a] = 0.0; 
    for (int b = 0; b < 20; ++b) { 
      G[vst_1x2v_ser_p2_ho_ph_v1_cmap[b]] += vst_1x2v_ser_p2_ho_ph_v1_Vd0[j*6 + vst_1x2v_ser_p2_ho_ph_v1_vrd0map[b]]*(vst_1x2v_ser_p2_ho_ph_v1_dcoefr0[b]*hamil[b]); 
    } 
    for (int a = 0; a < 3; ++a) dH_dvx += vst_1x2v_ser_p2_ho_ph_v1_Cm[i*3 + a]*G[a]; 
  } 
  const double *Bz = &qmem[15]; 
  double Bz_quad = 0.0; 
  for (int a = 0; a < 3; ++a) Bz_quad += vst_1x2v_ser_p2_ho_conf_ev[i*3 + a]*Bz[a]; 
  return -2.0/(dxv[1]*jacob_vel_surf[0])*dH_dvx*Bz_quad; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vy_1x2v_ser_p2(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += B_ho_hamil_phase_alpha_quad_vy_1x2v_ser_p2_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
