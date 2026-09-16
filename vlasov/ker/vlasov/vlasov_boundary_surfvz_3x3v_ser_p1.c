#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_3x3v_ser_p1.h> 
GKYL_CU_DH double vlasov_boundary_surfvz_3x3v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv12 = 2.0/dxv[5]; 
  for (int k = 0; k < 32; ++k) { 
  const int a = vst_3x3v_ser_p1_prj_v2_kamap[k]; 
  const int b = vst_3x3v_ser_p1_prj_v2_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 8; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 4; ++j) { 
      t += vst_3x3v_ser_p1_prj_v2_Vw[j*4 + b]*flux[64 + i*4 + j]; 
    } 
    g += vst_3x3v_ser_p1_prj_v2_Cw[i*8 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_3x3v_ser_p1_prj_v2_out_off[k]; q < vst_3x3v_ser_p1_prj_v2_out_off[k+1]; ++q) { 
      out[vst_3x3v_ser_p1_prj_v2_out_mode[q]] += dv12*vst_3x3v_ser_p1_prj_v2_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_3x3v_ser_p1_prj_v2_out_off[k]; q < vst_3x3v_ser_p1_prj_v2_out_off[k+1]; ++q) { 
      out[vst_3x3v_ser_p1_prj_v2_out_mode[q]] += dv12*vst_3x3v_ser_p1_prj_v2_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
