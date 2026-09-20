#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double vlasov_boundary_surfvx_3x3v_tensor_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[3]; 
  for (int k = 0; k < 72; ++k) { 
  const int a = vst_3x3v_tensor_p1_prj_v0_kamap[k]; 
  const int b = vst_3x3v_tensor_p1_prj_v0_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 8; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 9; ++j) { 
      t += vst_3x3v_tensor_p1_prj_v0_Vw[j*9 + b]*flux[0 + i*9 + j]; 
    } 
    g += vst_3x3v_tensor_p1_prj_v0_Cw[i*8 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_3x3v_tensor_p1_prj_v0_out_off[k]; q < vst_3x3v_tensor_p1_prj_v0_out_off[k+1]; ++q) { 
      out[vst_3x3v_tensor_p1_prj_v0_out_mode[q]] += dv10*vst_3x3v_tensor_p1_prj_v0_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_3x3v_tensor_p1_prj_v0_out_off[k]; q < vst_3x3v_tensor_p1_prj_v0_out_off[k+1]; ++q) { 
      out[vst_3x3v_tensor_p1_prj_v0_out_mode[q]] += dv10*vst_3x3v_tensor_p1_prj_v0_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
