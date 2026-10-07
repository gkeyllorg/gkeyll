#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double vlasov_ho_surfvy_1x3v_ser_p2(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  double dv11 = 2.0/dxv[2]; 
  for (int k = 0; k < 20; ++k) { 
  const double g_l = flux_l[20 + k]; 
  const double g_r = flux_r[20 + k]; 
  for (int q = vst_1x3v_ser_p2_ho_prj_v1_out_off[k]; q < vst_1x3v_ser_p2_ho_prj_v1_out_off[k+1]; ++q) { 
    out[vst_1x3v_ser_p2_ho_prj_v1_out_mode[q]] += dv11*(vst_1x3v_ser_p2_ho_prj_v1_out_cl[q]*g_l + vst_1x3v_ser_p2_ho_prj_v1_out_cr[q]*g_r); 
  } 
  } 
  return 0.0;
} 
