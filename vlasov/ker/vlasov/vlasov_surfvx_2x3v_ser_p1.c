#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x3v_ser_p1.h> 
GKYL_CU_DH double vlasov_surfvx_2x3v_ser_p1(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[2]; 
  for (int k = 0; k < 16; ++k) { 
  const int a = vst_2x3v_ser_p1_prj_v0_kamap[k]; 
  const int b = vst_2x3v_ser_p1_prj_v0_kbmap[k]; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    double t_l = 0.0; 
    double t_r = 0.0; 
    for (int j = 0; j < 4; ++j) { 
      t_l += vst_2x3v_ser_p1_prj_v0_Vw[j*4 + b]*Fhat_l_nodal[0 + i*4 + j]; 
      t_r += vst_2x3v_ser_p1_prj_v0_Vw[j*4 + b]*Fhat_r_nodal[0 + i*4 + j]; 
    } 
    g_l += vst_2x3v_ser_p1_prj_v0_Cw[i*4 + a]*t_l; 
    g_r += vst_2x3v_ser_p1_prj_v0_Cw[i*4 + a]*t_r; 
  } 
  for (int q = vst_2x3v_ser_p1_prj_v0_out_off[k]; q < vst_2x3v_ser_p1_prj_v0_out_off[k+1]; ++q) { 
    out[vst_2x3v_ser_p1_prj_v0_out_mode[q]] += dv10*(vst_2x3v_ser_p1_prj_v0_out_cl[q]*g_l + vst_2x3v_ser_p1_prj_v0_out_cr[q]*g_r); 
  } 
  } 
  return 0.0;
} 
