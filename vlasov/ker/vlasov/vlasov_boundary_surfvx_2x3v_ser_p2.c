#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x3v_ser_p2.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_2x3v_ser_p2(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[2]; 
  for (int k = 0; k < 48; ++k) { 
  const double g = flux[0 + k]; 
  if (edge == -1) { 
    for (int q = vst_2x3v_ser_p2_prj_v0_out_off[k]; q < vst_2x3v_ser_p2_prj_v0_out_off[k+1]; ++q) { 
      out[vst_2x3v_ser_p2_prj_v0_out_mode[q]] += dv10*vst_2x3v_ser_p2_prj_v0_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_2x3v_ser_p2_prj_v0_out_off[k]; q < vst_2x3v_ser_p2_prj_v0_out_off[k+1]; ++q) { 
      out[vst_2x3v_ser_p2_prj_v0_out_mode[q]] += dv10*vst_2x3v_ser_p2_prj_v0_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
