#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x1v_ser_p2.h> 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_1x1v_ser_p2(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  for (int k = 0; k < 3; ++k) { 
  const double g = flux[0 + k]; 
  if (edge == -1) { 
    for (int q = vst_1x1v_ser_p2_ho_prj_v0_out_off[k]; q < vst_1x1v_ser_p2_ho_prj_v0_out_off[k+1]; ++q) { 
      out[vst_1x1v_ser_p2_ho_prj_v0_out_mode[q]] += dv10*vst_1x1v_ser_p2_ho_prj_v0_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_1x1v_ser_p2_ho_prj_v0_out_off[k]; q < vst_1x1v_ser_p2_ho_prj_v0_out_off[k+1]; ++q) { 
      out[vst_1x1v_ser_p2_ho_prj_v0_out_mode[q]] += dv10*vst_1x1v_ser_p2_ho_prj_v0_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
