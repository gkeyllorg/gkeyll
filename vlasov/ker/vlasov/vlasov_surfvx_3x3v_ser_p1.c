#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_3x3v_ser_p1.h> 
GKYL_CU_DH double vlasov_surfvx_3x3v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[3]; 
  for (int k = 0; k < 32; ++k) { 
  const double g_l = flux_l[0 + k]; 
  const double g_r = flux_r[0 + k]; 
  for (int q = vst_3x3v_ser_p1_prj_v0_out_off[k]; q < vst_3x3v_ser_p1_prj_v0_out_off[k+1]; ++q) { 
    out[vst_3x3v_ser_p1_prj_v0_out_mode[q]] += dv10*(vst_3x3v_ser_p1_prj_v0_out_cl[q]*g_l + vst_3x3v_ser_p1_prj_v0_out_cr[q]*g_r); 
  } 
  } 
  return 0.0;
} 
