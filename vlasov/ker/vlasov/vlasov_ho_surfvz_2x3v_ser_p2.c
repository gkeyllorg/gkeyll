#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x3v_ser_p2.h> 
GKYL_CU_DH double vlasov_ho_surfvz_2x3v_ser_p2(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dv12 = 2.0/dxv[4]; 
  for (int k = 0; k < 48; ++k) { 
  const int a = vst_2x3v_ser_p2_ho_prj_v2_kamap[k]; 
  const int b = vst_2x3v_ser_p2_ho_prj_v2_kbmap[k]; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int i = 0; i < 16; ++i) { 
    double t_l = 0.0; 
    double t_r = 0.0; 
    for (int j = 0; j < 16; ++j) { 
      t_l += vst_2x3v_ser_p2_ho_prj_v2_Vw[j*8 + b]*Fhat_l_nodal[512 + i*16 + j]; 
      t_r += vst_2x3v_ser_p2_ho_prj_v2_Vw[j*8 + b]*Fhat_r_nodal[512 + i*16 + j]; 
    } 
    g_l += vst_2x3v_ser_p2_ho_prj_v2_Cw[i*8 + a]*t_l; 
    g_r += vst_2x3v_ser_p2_ho_prj_v2_Cw[i*8 + a]*t_r; 
  } 
  for (int q = vst_2x3v_ser_p2_ho_prj_v2_out_off[k]; q < vst_2x3v_ser_p2_ho_prj_v2_out_off[k+1]; ++q) { 
    out[vst_2x3v_ser_p2_ho_prj_v2_out_mode[q]] += dv12*(vst_2x3v_ser_p2_ho_prj_v2_out_cl[q]*g_l + vst_2x3v_ser_p2_ho_prj_v2_out_cr[q]*g_r); 
  } 
  } 
  return 0.0;
} 
