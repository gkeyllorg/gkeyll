#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double vlasov_ho_surfvz_3x3v_tensor_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  double dv12 = 2.0/dxv[5]; 
  for (int k = 0; k < 72; ++k) { 
  const double g_l = flux_l[144 + k]; 
  const double g_r = flux_r[144 + k]; 
  for (int q = vst_3x3v_tensor_p1_ho_prj_v2_out_off[k]; q < vst_3x3v_tensor_p1_ho_prj_v2_out_off[k+1]; ++q) { 
    out[vst_3x3v_tensor_p1_ho_prj_v2_out_mode[q]] += dv12*(vst_3x3v_tensor_p1_ho_prj_v2_out_cl[q]*g_l + vst_3x3v_tensor_p1_ho_prj_v2_out_cr[q]*g_r); 
  } 
  } 
  return 0.0;
} 
