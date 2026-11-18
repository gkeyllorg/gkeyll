#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH double vlasov_boundary_ho_surfy_2x3v_tensor_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dx11 = 2.0/dxv[1]; 
  for (int k = 0; k < 54; ++k) { 
  const double g = flux[54 + k]; 
  if (edge == -1) { 
    for (int q = vst_2x3v_tensor_p1_ho_prj_x1_out_off[k]; q < vst_2x3v_tensor_p1_ho_prj_x1_out_off[k+1]; ++q) { 
      out[vst_2x3v_tensor_p1_ho_prj_x1_out_mode[q]] += dx11*vst_2x3v_tensor_p1_ho_prj_x1_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_2x3v_tensor_p1_ho_prj_x1_out_off[k]; q < vst_2x3v_tensor_p1_ho_prj_x1_out_off[k+1]; ++q) { 
      out[vst_2x3v_tensor_p1_ho_prj_x1_out_mode[q]] += dx11*vst_2x3v_tensor_p1_ho_prj_x1_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
