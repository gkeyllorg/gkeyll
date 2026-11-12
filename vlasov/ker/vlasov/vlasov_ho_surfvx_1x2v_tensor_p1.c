#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x2v_tensor_p1.h> 
GKYL_CU_DH double vlasov_ho_surfvx_1x2v_tensor_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  for (int k = 0; k < 6; ++k) { 
  const double g_l = flux_l[0 + k]; 
  const double g_r = flux_r[0 + k]; 
  for (int q = vst_1x2v_tensor_p1_ho_prj_v0_out_off[k]; q < vst_1x2v_tensor_p1_ho_prj_v0_out_off[k+1]; ++q) { 
    out[vst_1x2v_tensor_p1_ho_prj_v0_out_mode[q]] += dv10*(vst_1x2v_tensor_p1_ho_prj_v0_out_cl[q]*g_l + vst_1x2v_tensor_p1_ho_prj_v0_out_cr[q]*g_r); 
  } 
  } 
  return 0.0;
} 
