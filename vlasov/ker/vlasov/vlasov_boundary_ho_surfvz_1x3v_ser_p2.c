#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double vlasov_boundary_ho_surfvz_1x3v_ser_p2(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv12 = 2.0/dxv[3]; 
  for (int k = 0; k < 20; ++k) { 
  const double g = flux[40 + k]; 
  if (edge == -1) { 
    for (int q = vst_1x3v_ser_p2_ho_prj_v2_out_off[k]; q < vst_1x3v_ser_p2_ho_prj_v2_out_off[k+1]; ++q) { 
      out[vst_1x3v_ser_p2_ho_prj_v2_out_mode[q]] += dv12*vst_1x3v_ser_p2_ho_prj_v2_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_1x3v_ser_p2_ho_prj_v2_out_off[k]; q < vst_1x3v_ser_p2_ho_prj_v2_out_off[k+1]; ++q) { 
      out[vst_1x3v_ser_p2_ho_prj_v2_out_mode[q]] += dv12*vst_1x3v_ser_p2_ho_prj_v2_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
