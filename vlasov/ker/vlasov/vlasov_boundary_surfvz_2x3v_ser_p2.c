#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x3v_ser_p2.h> 
GKYL_CU_DH double vlasov_boundary_surfvz_2x3v_ser_p2(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv12 = 2.0/dxv[4]; 
  for (int k = 0; k < 48; ++k) { 
  const int a = vst_2x3v_ser_p2_prj_v2_kamap[k]; 
  const int b = vst_2x3v_ser_p2_prj_v2_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 9; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 9; ++j) { 
      t += vst_2x3v_ser_p2_prj_v2_Vw[j*8 + b]*flux[162 + i*9 + j]; 
    } 
    g += vst_2x3v_ser_p2_prj_v2_Cw[i*8 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_2x3v_ser_p2_prj_v2_out_off[k]; q < vst_2x3v_ser_p2_prj_v2_out_off[k+1]; ++q) { 
      out[vst_2x3v_ser_p2_prj_v2_out_mode[q]] += dv12*vst_2x3v_ser_p2_prj_v2_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_2x3v_ser_p2_prj_v2_out_off[k]; q < vst_2x3v_ser_p2_prj_v2_out_off[k+1]; ++q) { 
      out[vst_2x3v_ser_p2_prj_v2_out_mode[q]] += dv12*vst_2x3v_ser_p2_prj_v2_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
