#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH double vlasov_boundary_ho_surfvy_2x3v_tensor_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv11 = 2.0/dxv[3]; 
  for (int k = 0; k < 36; ++k) { 
  const double g = flux[36 + k]; 
  if (edge == -1) { 
    for (int q = vst_2x3v_tensor_p1_ho_prj_v1_out_off[k]; q < vst_2x3v_tensor_p1_ho_prj_v1_out_off[k+1]; ++q) { 
      out[vst_2x3v_tensor_p1_ho_prj_v1_out_mode[q]] += dv11*vst_2x3v_tensor_p1_ho_prj_v1_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_2x3v_tensor_p1_ho_prj_v1_out_off[k]; q < vst_2x3v_tensor_p1_ho_prj_v1_out_off[k+1]; ++q) { 
      out[vst_2x3v_tensor_p1_ho_prj_v1_out_mode[q]] += dv11*vst_2x3v_tensor_p1_ho_prj_v1_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
