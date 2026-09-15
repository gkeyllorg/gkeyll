#include <gkyl_vlasov_kernels.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x1v_ser_p3.h> 
GKYL_CU_DH double vlasov_boundary_ho_surfvx_1x1v_ser_p3(const double *w, const double *dxv,
  const int edge, const double *flux, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  for (int k = 0; k < 4; ++k) { 
  const int a = vst_1x1v_ser_p3_ho_prj_v0_kamap[k]; 
  const int b = vst_1x1v_ser_p3_ho_prj_v0_kbmap[k]; 
  double g = 0.0; 
  for (int i = 0; i < 5; ++i) { 
    double t = 0.0; 
    for (int j = 0; j < 1; ++j) { 
      t += vst_1x1v_ser_p3_ho_prj_v0_Vw[j*1 + b]*flux[0 + i*1 + j]; 
    } 
    g += vst_1x1v_ser_p3_ho_prj_v0_Cw[i*4 + a]*t; 
  } 
  if (edge == -1) { 
    for (int q = vst_1x1v_ser_p3_ho_prj_v0_out_off[k]; q < vst_1x1v_ser_p3_ho_prj_v0_out_off[k+1]; ++q) { 
      out[vst_1x1v_ser_p3_ho_prj_v0_out_mode[q]] += dv10*vst_1x1v_ser_p3_ho_prj_v0_out_cr[q]*g; 
    } 
  } else { 
    for (int q = vst_1x1v_ser_p3_ho_prj_v0_out_off[k]; q < vst_1x1v_ser_p3_ho_prj_v0_out_off[k+1]; ++q) { 
      out[vst_1x1v_ser_p3_ho_prj_v0_out_mode[q]] += dv10*vst_1x1v_ser_p3_ho_prj_v0_out_cl[q]*g; 
    } 
  } 
  } 
  return 0.0;
} 
