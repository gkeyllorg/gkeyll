#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double vlasov_surfz_3x3v_tensor_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  double dx12 = 2.0/dxv[2]; 
  for (int k = 0; k < 108; ++k) { 
  const double g_l = flux_l[216 + k]; 
  const double g_r = flux_r[216 + k]; 
  for (int q = vst_3x3v_tensor_p1_prj_x2_out_off[k]; q < vst_3x3v_tensor_p1_prj_x2_out_off[k+1]; ++q) { 
    out[vst_3x3v_tensor_p1_prj_x2_out_mode[q]] += dx12*(vst_3x3v_tensor_p1_prj_x2_out_cl[q]*g_l + vst_3x3v_tensor_p1_prj_x2_out_cr[q]*g_r); 
  } 
  } 
  return 0.0;
} 
