#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double vlasov_surfvz_3x3v_tensor_p1(const double *w, const double *dxv,
  const double *Fhat_l_nodal, const double *Fhat_r_nodal, double* GKYL_RESTRICT out) 
{ 
  double dv12 = 2.0/dxv[5]; 
  for (int k = 0; k < 72; ++k) { 
  const int a = vst_3x3v_tensor_p1_prj_v2_kamap[k]; 
  const int b = vst_3x3v_tensor_p1_prj_v2_kbmap[k]; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int i = 0; i < 8; ++i) { 
    double t_l = 0.0; 
    double t_r = 0.0; 
    for (int j = 0; j < 9; ++j) { 
      t_l += vst_3x3v_tensor_p1_prj_v2_Vw[j*9 + b]*Fhat_l_nodal[144 + i*9 + j]; 
      t_r += vst_3x3v_tensor_p1_prj_v2_Vw[j*9 + b]*Fhat_r_nodal[144 + i*9 + j]; 
    } 
    g_l += vst_3x3v_tensor_p1_prj_v2_Cw[i*8 + a]*t_l; 
    g_r += vst_3x3v_tensor_p1_prj_v2_Cw[i*8 + a]*t_r; 
  } 
  for (int q = vst_3x3v_tensor_p1_prj_v2_out_off[k]; q < vst_3x3v_tensor_p1_prj_v2_out_off[k+1]; ++q) { 
    out[vst_3x3v_tensor_p1_prj_v2_out_mode[q]] += dv12*(vst_3x3v_tensor_p1_prj_v2_out_cl[q]*g_l + vst_3x3v_tensor_p1_prj_v2_out_cr[q]*g_r); 
  } 
  } 
  return 0.0;
} 
