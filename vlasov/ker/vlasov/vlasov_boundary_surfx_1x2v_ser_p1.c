#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x2v_ser_p1.h> 
GKYL_CU_DH double vlasov_boundary_surfx_1x2v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  for (int k = 0; k < 4; ++k) { 
  const int a = vst_1x2v_ser_p1_prj_x0_kamap[k]; 
  const int b = vst_1x2v_ser_p1_prj_x0_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 1; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 4; ++j) { 
      t += vst_1x2v_ser_p1_prj_x0_Vw[j*4 + b]*flux[0 + i*4 + j]; 
    } 
    g += vst_1x2v_ser_p1_prj_x0_Cw[i*1 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_1x2v_ser_p1_prj_x0_out_off[k]; q < vst_1x2v_ser_p1_prj_x0_out_off[k+1]; ++q) { 
      out[vst_1x2v_ser_p1_prj_x0_out_mode[q]] += dx10*vst_1x2v_ser_p1_prj_x0_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_1x2v_ser_p1_prj_x0_out_off[k]; q < vst_1x2v_ser_p1_prj_x0_out_off[k+1]; ++q) { 
      out[vst_1x2v_ser_p1_prj_x0_out_mode[q]] += dx10*vst_1x2v_ser_p1_prj_x0_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
