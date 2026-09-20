#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double vlasov_boundary_ho_surfvz_1x3v_ser_p2(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv12 = 2.0/dxv[3]; 
  for (int k = 0; k < 20; ++k) { 
  const int a = vst_1x3v_ser_p2_ho_prj_v2_kamap[k]; 
  const int b = vst_1x3v_ser_p2_ho_prj_v2_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 16; ++j) { 
      t += vst_1x3v_ser_p2_ho_prj_v2_Vw[j*8 + b]*flux[128 + i*16 + j]; 
    } 
    g += vst_1x3v_ser_p2_ho_prj_v2_Cw[i*3 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_1x3v_ser_p2_ho_prj_v2_out_off[k]; q < vst_1x3v_ser_p2_ho_prj_v2_out_off[k+1]; ++q) { 
      out[vst_1x3v_ser_p2_ho_prj_v2_out_mode[q]] += dv12*vst_1x3v_ser_p2_ho_prj_v2_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_1x3v_ser_p2_ho_prj_v2_out_off[k]; q < vst_1x3v_ser_p2_ho_prj_v2_out_off[k+1]; ++q) { 
      out[vst_1x3v_ser_p2_ho_prj_v2_out_mode[q]] += dv12*vst_1x3v_ser_p2_ho_prj_v2_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
