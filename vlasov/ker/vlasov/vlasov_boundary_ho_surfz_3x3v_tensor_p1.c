#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double vlasov_boundary_ho_surfz_3x3v_tensor_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dx12 = 2.0/dxv[2]; 
  for (int k = 0; k < 108; ++k) { 
  const double g = flux[216 + k]; 
  if (edge == -1) { 
    for (int q = vst_3x3v_tensor_p1_ho_prj_x2_out_off[k]; q < vst_3x3v_tensor_p1_ho_prj_x2_out_off[k+1]; ++q) { 
      out[vst_3x3v_tensor_p1_ho_prj_x2_out_mode[q]] += dx12*vst_3x3v_tensor_p1_ho_prj_x2_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_3x3v_tensor_p1_ho_prj_x2_out_off[k]; q < vst_3x3v_tensor_p1_ho_prj_x2_out_off[k+1]; ++q) { 
      out[vst_3x3v_tensor_p1_ho_prj_x2_out_mode[q]] += dx12*vst_3x3v_tensor_p1_ho_prj_x2_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
