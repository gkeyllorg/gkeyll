#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_flux_surf_tables_1x1v_p1.h> 
GKYL_CU_DH double pkpm_stream_surfx_1x1v_ser_p1(const double *w, const double *dxv,
  const double *flux_l, const double *flux_r, double* GKYL_RESTRICT out) 
{ 
  // flux_l/flux_r: Stored modal fluxes [F, G, M0_perp, N_perp] on the lower/upper face of the cell. 
  // out:           Incremented output [F, G, M0_perp, N_perp]. 
  const double dx1 = 2.0/dxv[0]; 
  for (int c = 0; c < 8; ++c) { 
    const double *fl = &flux_l[c*3]; 
    const double *fr = &flux_r[c*3]; 
    double *o = &out[c*6]; 
    for (int k = 0; k < 3; ++k) { 
      const double g_l = fl[k]; 
      const double g_r = fr[k]; 
      for (int q = vst_pkpm_1x1v_p1_prj_x0_out_off[k]; q < vst_pkpm_1x1v_p1_prj_x0_out_off[k+1]; ++q) { 
        o[vst_pkpm_1x1v_p1_prj_x0_out_mode[q]] += dx1*(vst_pkpm_1x1v_p1_prj_x0_out_cl[q]*g_l + vst_pkpm_1x1v_p1_prj_x0_out_cr[q]*g_r); 
      } 
    } 
  } 
  return 0.0; 
} 
