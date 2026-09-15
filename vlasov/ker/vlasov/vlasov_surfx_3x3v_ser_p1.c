#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_3x3v_ser_p1.h> 
GKYL_CU_DH double vlasov_surfx_3x3v_ser_p1(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  for (int k = 0; k < 32; ++k) { 
  const int a = vst_3x3v_ser_p1_prj_x0_kamap[k]; 
  const int b = vst_3x3v_ser_p1_prj_x0_kbmap[k]; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    double t_l = 0.0; 
    double t_r = 0.0; 
    for (int j = 0; j < 8; ++j) { 
      t_l += vst_3x3v_ser_p1_prj_x0_Vw[j*8 + b]*Fhat_l_nodal[0 + i*8 + j]; 
      t_r += vst_3x3v_ser_p1_prj_x0_Vw[j*8 + b]*Fhat_r_nodal[0 + i*8 + j]; 
    } 
    g_l += vst_3x3v_ser_p1_prj_x0_Cw[i*4 + a]*t_l; 
    g_r += vst_3x3v_ser_p1_prj_x0_Cw[i*4 + a]*t_r; 
  } 
  for (int q = vst_3x3v_ser_p1_prj_x0_out_off[k]; q < vst_3x3v_ser_p1_prj_x0_out_off[k+1]; ++q) { 
    out[vst_3x3v_ser_p1_prj_x0_out_mode[q]] += dx10*(vst_3x3v_ser_p1_prj_x0_out_cl[q]*g_l + vst_3x3v_ser_p1_prj_x0_out_cr[q]*g_r); 
  } 
  } 
  return 0.0;
} 
