#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH double vlasov_boundary_ho_surfx_2x3v_tensor_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  for (int k = 0; k < 54; ++k) { 
  const int a = vst_2x3v_tensor_p1_ho_prj_x0_kamap[k]; 
  const int b = vst_2x3v_tensor_p1_ho_prj_x0_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 2; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 64; ++j) { 
      t += vst_2x3v_tensor_p1_ho_prj_x0_Vw[j*27 + b]*flux[0 + i*64 + j]; 
    } 
    g += vst_2x3v_tensor_p1_ho_prj_x0_Cw[i*2 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_2x3v_tensor_p1_ho_prj_x0_out_off[k]; q < vst_2x3v_tensor_p1_ho_prj_x0_out_off[k+1]; ++q) { 
      out[vst_2x3v_tensor_p1_ho_prj_x0_out_mode[q]] += dx10*vst_2x3v_tensor_p1_ho_prj_x0_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_2x3v_tensor_p1_ho_prj_x0_out_off[k]; q < vst_2x3v_tensor_p1_ho_prj_x0_out_off[k+1]; ++q) { 
      out[vst_2x3v_tensor_p1_ho_prj_x0_out_mode[q]] += dx10*vst_2x3v_tensor_p1_ho_prj_x0_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
