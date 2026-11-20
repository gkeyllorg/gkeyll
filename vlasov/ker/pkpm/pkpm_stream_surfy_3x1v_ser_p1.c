#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_flux_surf_tables_3x1v_p1.h> 
GKYL_CU_DH double pkpm_stream_surfy_3x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  // flux_l/flux_r: Stored modal fluxes [F, G, M0_perp, N_perp] on the lower/upper face of the cell. 
  // out:           Incremented output [F, G, M0_perp, N_perp]. 
  const double dx1 = 2.0/dxv[1]; 
  for (int c = 0; c < 8; ++c) { 
    const double *fl = &flux_l[c*12]; 
    const double *fr = &flux_r[c*12]; 
    double *o = &out[c*24]; 
    for (int k = 0; k < 12; ++k) { 
      const double g_l = fl[k]; 
      const double g_r = fr[k]; 
      for (int q = vst_pkpm_3x1v_p1_prj_x1_out_off[k]; q < vst_pkpm_3x1v_p1_prj_x1_out_off[k+1]; ++q) { 
        o[vst_pkpm_3x1v_p1_prj_x1_out_mode[q]] += dx1*(vst_pkpm_3x1v_p1_prj_x1_out_cl[q]*g_l + vst_pkpm_3x1v_p1_prj_x1_out_cr[q]*g_r); 
      } 
    } 
  } 
  return 0.0; 
} 
