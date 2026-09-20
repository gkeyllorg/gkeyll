#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x2v_tensor_p2.h> 
GKYL_CU_DH double vlasov_boundary_surfvy_1x2v_tensor_p2(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv11 = 2.0/dxv[2]; 
  for (int k = 0; k < 9; ++k) { 
  const int a = vst_1x2v_tensor_p2_prj_v1_kamap[k]; 
  const int b = vst_1x2v_tensor_p2_prj_v1_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 4; ++j) { 
      t += vst_1x2v_tensor_p2_prj_v1_Vw[j*3 + b]*flux[16 + i*4 + j]; 
    } 
    g += vst_1x2v_tensor_p2_prj_v1_Cw[i*3 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_1x2v_tensor_p2_prj_v1_out_off[k]; q < vst_1x2v_tensor_p2_prj_v1_out_off[k+1]; ++q) { 
      out[vst_1x2v_tensor_p2_prj_v1_out_mode[q]] += dv11*vst_1x2v_tensor_p2_prj_v1_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_1x2v_tensor_p2_prj_v1_out_off[k]; q < vst_1x2v_tensor_p2_prj_v1_out_off[k+1]; ++q) { 
      out[vst_1x2v_tensor_p2_prj_v1_out_mode[q]] += dv11*vst_1x2v_tensor_p2_prj_v1_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
