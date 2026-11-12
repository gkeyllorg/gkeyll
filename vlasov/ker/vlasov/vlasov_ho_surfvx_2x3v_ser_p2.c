#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x3v_ser_p2.h> 
GKYL_CU_DH double vlasov_ho_surfvx_2x3v_ser_p2(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[2]; 
  for (int k = 0; k < 48; ++k) { 
  const double g_l = flux_l[0 + k]; 
  const double g_r = flux_r[0 + k]; 
  for (int q = vst_2x3v_ser_p2_ho_prj_v0_out_off[k]; q < vst_2x3v_ser_p2_ho_prj_v0_out_off[k+1]; ++q) { 
    out[vst_2x3v_ser_p2_ho_prj_v0_out_mode[q]] += dv10*(vst_2x3v_ser_p2_ho_prj_v0_out_cl[q]*g_l + vst_2x3v_ser_p2_ho_prj_v0_out_cr[q]*g_r); 
  } 
  } 
  return 0.0;
} 
