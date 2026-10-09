#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH void ho_lax_flux_nodal_x_1x3v_ser_p2_g(int item, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (item >= 64) return; 
  const int a = item/64; 
  const int j = item - a*64; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_1x3v_ser_p2_ho_ph_x0_aoff[a]; q < vst_1x3v_ser_p2_ho_ph_x0_aoff[a+1]; ++q) { 
    const int k = vst_1x3v_ser_p2_ho_ph_x0_aks[q]; 
    g_l += vst_1x3v_ser_p2_ho_ph_x0_Wl[k*64 + j]*f_l[k]; 
    g_r += vst_1x3v_ser_p2_ho_ph_x0_Wr[k*64 + j]*f_r[k]; 
  } 
  G_l[j*1 + a] = g_l; 
  G_r[j*1 + a] = g_r; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_x_1x3v_ser_p2_node(int i, int j, const double *jacob_pos_l, const double *jacob_pos_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 1; ++a) { 
    f_l_quad += vst_1x3v_ser_p2_ho_ph_x0_Cm[i*1 + a]*G_l[j*1 + a]; 
    f_r_quad += vst_1x3v_ser_p2_ho_ph_x0_Cm[i*1 + a]*G_r[j*1 + a]; 
  } 
  f_l_quad *= 1.0/jacob_pos_l[0]; 
  f_r_quad *= 1.0/jacob_pos_r[0]; 
  const int n = i*64 + j; 
  Fhat_nodal[n] = 0.5*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void ho_lax_flux_nodal_x_1x3v_ser_p2_prj(int unit, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (unit >= 20) return; 
  const int b = unit; 
  double t[1]; 
  for (int i = 0; i < 1; ++i) t[i] = 0.0; 
  for (int j = 0; j < 64; ++j) { 
    const double w = vst_1x3v_ser_p2_ho_prj_x0_Vw[j*20 + b]; 
    for (int i = 0; i < 1; ++i) t[i] += w*Fhat_nodal[i*64 + j]; 
  } 
  for (int q = vst_1x3v_ser_p2_ho_prj_x0_boff[b]; q < vst_1x3v_ser_p2_ho_prj_x0_boff[b+1]; ++q) { 
    const int k = vst_1x3v_ser_p2_ho_prj_x0_bks[q]; 
    const int a = vst_1x3v_ser_p2_ho_prj_x0_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 1; ++i) g += vst_1x3v_ser_p2_ho_prj_x0_Cw[i*1 + a]*t[i]; 
    flux[0 + k] = g; 
  } 
} 

GKYL_CU_DH double ho_lax_flux_nodal_x_1x3v_ser_p2_cfl(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r, double alpha_max) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double jacob_pos_min = fmin(jacob_pos_l[0], jacob_pos_r[0]); 
  return 2.5*dx10*alpha_max/jacob_pos_min;
} 

GKYL_CU_DH double ho_lax_flux_nodal_x_1x3v_ser_p2(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[64]; 
  double G_r[64]; 
  double Fhat_nodal[64]; 
  for (int item = 0; item < 64; ++item) ho_lax_flux_nodal_x_1x3v_ser_p2_g(item, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 1; ++i) { 
    for (int j = 0; j < 64; ++j) { 
      alpha_max = fmax(alpha_max, ho_lax_flux_nodal_x_1x3v_ser_p2_node(i, j, jacob_pos_l, jacob_pos_r, alpha_quad[i*64 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int unit = 0; unit < 64; ++unit) ho_lax_flux_nodal_x_1x3v_ser_p2_prj(unit, Fhat_nodal, flux); 
  return ho_lax_flux_nodal_x_1x3v_ser_p2_cfl(dxv, jacob_pos_l, jacob_pos_r, alpha_max); 
} 
