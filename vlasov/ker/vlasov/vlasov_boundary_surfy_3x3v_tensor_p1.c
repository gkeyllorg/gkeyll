#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double vlasov_boundary_surfy_3x3v_tensor_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dx11 = 2.0/dxv[1]; 
  for (int k = 0; k < 108; ++k) { 
  const int a = vst_3x3v_tensor_p1_prj_x1_kamap[k]; 
  const int b = vst_3x3v_tensor_p1_prj_x1_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 27; ++j) { 
      t += vst_3x3v_tensor_p1_prj_x1_Vw[j*27 + b]*flux[108 + i*27 + j]; 
    } 
    g += vst_3x3v_tensor_p1_prj_x1_Cw[i*4 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_3x3v_tensor_p1_prj_x1_out_off[k]; q < vst_3x3v_tensor_p1_prj_x1_out_off[k+1]; ++q) { 
      out[vst_3x3v_tensor_p1_prj_x1_out_mode[q]] += dx11*vst_3x3v_tensor_p1_prj_x1_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_3x3v_tensor_p1_prj_x1_out_off[k]; q < vst_3x3v_tensor_p1_prj_x1_out_off[k+1]; ++q) { 
      out[vst_3x3v_tensor_p1_prj_x1_out_mode[q]] += dx11*vst_3x3v_tensor_p1_prj_x1_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
