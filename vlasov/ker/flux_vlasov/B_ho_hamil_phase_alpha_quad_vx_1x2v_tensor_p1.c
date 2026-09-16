#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_tensor_p1.h> 
GKYL_CU_DH double B_ho_hamil_phase_alpha_quad_vx_1x2v_tensor_p1_node(int i, int j, const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem) 
{ 
  double dH_dvy = 0.0; 
  { 
    double G[2]; 
    for (int a = 0; a < 2; ++a) G[a] = 0.0; 
    for (int b = 0; b < 18; ++b) { 
      G[vst_1x2v_tensor_p1_ho_ph_v0_cmap[b]] += vst_1x2v_tensor_p1_ho_ph_v0_Vd1[j*7 + vst_1x2v_tensor_p1_ho_ph_v0_vrd1map[b]]*(vst_1x2v_tensor_p1_ho_ph_v0_dcoefr1[b]*hamil[b]); 
    } 
    for (int a = 0; a < 2; ++a) dH_dvy += vst_1x2v_tensor_p1_ho_ph_v0_Cm[i*2 + a]*G[a]; 
  } 
  const double *Bz = &qmem[10]; 
  double Bz_quad = 0.0; 
  for (int a = 0; a < 2; ++a) Bz_quad += vst_1x2v_tensor_p1_ho_conf_ev[i*2 + a]*Bz[a]; 
  return 2.0/(dxv[2]*jacob_vel_surf[4 + j])*dH_dvy*Bz_quad; 
} 

GKYL_CU_DH void B_ho_hamil_phase_alpha_quad_vx_1x2v_tensor_p1(const double *dxv, const double *jacob_vel_surf,
  const double *hamil, const double *qmem, double* GKYL_RESTRICT alpha_quad) 
{ 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 4; ++j) alpha_quad[i*4 + j] += B_ho_hamil_phase_alpha_quad_vx_1x2v_tensor_p1_node(i, j, dxv, jacob_vel_surf, hamil, qmem); 
  } 
} 
