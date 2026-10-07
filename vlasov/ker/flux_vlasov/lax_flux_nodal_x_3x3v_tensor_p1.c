#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH void lax_flux_nodal_x_3x3v_tensor_p1_g(int item, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (item >= 108) return; 
  const int j = item/4; 
  const int a = item - j*4; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_3x3v_tensor_p1_ph_x0_aoff[a]; q < vst_3x3v_tensor_p1_ph_x0_aoff[a+1]; ++q) { 
    const int k = vst_3x3v_tensor_p1_ph_x0_aks[q]; 
    g_l += vst_3x3v_tensor_p1_ph_x0_V[j*54 + vst_3x3v_tensor_p1_ph_x0_vlmap[k]]*(vst_3x3v_tensor_p1_ph_x0_coefl[k]*f_l[k]); 
    g_r += vst_3x3v_tensor_p1_ph_x0_V[j*54 + vst_3x3v_tensor_p1_ph_x0_vrmap[k]]*(vst_3x3v_tensor_p1_ph_x0_coefr[k]*f_r[k]); 
  } 
  G_l[item] = g_l; 
  G_r[item] = g_r; 
} 

GKYL_CU_DH double lax_flux_nodal_x_3x3v_tensor_p1_node(int i, int j, const double *jacob_pos_l, const double *jacob_pos_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    f_l_quad += vst_3x3v_tensor_p1_ph_x0_Cm[i*4 + a]*G_l[j*4 + a]; 
    f_r_quad += vst_3x3v_tensor_p1_ph_x0_Cm[i*4 + a]*G_r[j*4 + a]; 
  } 
  f_l_quad *= 1.0/jacob_pos_l[0]; 
  f_r_quad *= 1.0/jacob_pos_r[0]; 
  const int n = i*27 + j; 
  Fhat_nodal[n] = 0.5*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void lax_flux_nodal_x_3x3v_tensor_p1_prj(int k, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (k >= 108) return; 
  const int a = vst_3x3v_tensor_p1_prj_x0_kamap[k]; 
  const int b = vst_3x3v_tensor_p1_prj_x0_kbmap[k]; 
  double t[4]; 
  for (int i = 0; i < 4; ++i) t[i] = 0.0; 
  for (int j = 0; j < 27; ++j) { 
    const double w = vst_3x3v_tensor_p1_prj_x0_Vw[j*27 + b]; 
    for (int i = 0; i < 4; ++i) t[i] += w*Fhat_nodal[i*27 + j]; 
  } 
  double g = 0.0; 
  for (int i = 0; i < 4; ++i) g += vst_3x3v_tensor_p1_prj_x0_Cw[i*4 + a]*t[i]; 
  flux[0 + k] = g; 
} 

GKYL_CU_DH double lax_flux_nodal_x_3x3v_tensor_p1_cfl(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r, double alpha_max) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double jacob_pos_min = fmin(jacob_pos_l[0], jacob_pos_r[0]); 
  return 1.5*dx10*alpha_max/jacob_pos_min;
} 

GKYL_CU_DH double lax_flux_nodal_x_3x3v_tensor_p1(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[108]; 
  double G_r[108]; 
  double Fhat_nodal[108]; 
  for (int item = 0; item < 108; ++item) lax_flux_nodal_x_3x3v_tensor_p1_g(item, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 27; ++j) { 
      alpha_max = fmax(alpha_max, lax_flux_nodal_x_3x3v_tensor_p1_node(i, j, jacob_pos_l, jacob_pos_r, alpha_quad[i*27 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int k = 0; k < 108; ++k) lax_flux_nodal_x_3x3v_tensor_p1_prj(k, Fhat_nodal, flux); 
  return lax_flux_nodal_x_3x3v_tensor_p1_cfl(dxv, jacob_pos_l, jacob_pos_r, alpha_max); 
} 
