#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x1v_ser_p1.h> 
GKYL_CU_DH double vlasov_boundary_surfx_1x1v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  for (int k = 0; k < 2; ++k) { 
  const double g = flux[0 + k]; 
  if (edge == -1) { 
    for (int q = vst_1x1v_ser_p1_prj_x0_out_off[k]; q < vst_1x1v_ser_p1_prj_x0_out_off[k+1]; ++q) { 
      out[vst_1x1v_ser_p1_prj_x0_out_mode[q]] += dx10*vst_1x1v_ser_p1_prj_x0_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_1x1v_ser_p1_prj_x0_out_off[k]; q < vst_1x1v_ser_p1_prj_x0_out_off[k+1]; ++q) { 
      out[vst_1x1v_ser_p1_prj_x0_out_mode[q]] += dx10*vst_1x1v_ser_p1_prj_x0_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
