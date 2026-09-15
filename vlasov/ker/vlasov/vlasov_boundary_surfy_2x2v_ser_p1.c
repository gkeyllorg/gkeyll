#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x2v_ser_p1.h> 
GKYL_CU_DH double vlasov_boundary_surfy_2x2v_ser_p1(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dx11 = 2.0/dxv[1]; 
  for (int k = 0; k < 8; ++k) { 
  const int a = vst_2x2v_ser_p1_prj_x1_kamap[k]; 
  const int b = vst_2x2v_ser_p1_prj_x1_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 2; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 4; ++j) { 
      t += vst_2x2v_ser_p1_prj_x1_Vw[j*4 + b]*flux[8 + i*4 + j]; 
    } 
    g += vst_2x2v_ser_p1_prj_x1_Cw[i*2 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_2x2v_ser_p1_prj_x1_out_off[k]; q < vst_2x2v_ser_p1_prj_x1_out_off[k+1]; ++q) { 
      out[vst_2x2v_ser_p1_prj_x1_out_mode[q]] += dx11*vst_2x2v_ser_p1_prj_x1_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_2x2v_ser_p1_prj_x1_out_off[k]; q < vst_2x2v_ser_p1_prj_x1_out_off[k+1]; ++q) { 
      out[vst_2x2v_ser_p1_prj_x1_out_mode[q]] += dx11*vst_2x2v_ser_p1_prj_x1_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
