#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH void ho_lax_flux_nodal_y_2x2v_tensor_p1_g(int item, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (item >= 32) return; 
  const int a = item/16; 
  const int j = item - a*16; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_2x2v_tensor_p1_ho_ph_x1_aoff[a]; q < vst_2x2v_tensor_p1_ho_ph_x1_aoff[a+1]; ++q) { 
    const int k = vst_2x2v_tensor_p1_ho_ph_x1_aks[q]; 
    g_l += vst_2x2v_tensor_p1_ho_ph_x1_Wl[k*16 + j]*f_l[k]; 
    g_r += vst_2x2v_tensor_p1_ho_ph_x1_Wr[k*16 + j]*f_r[k]; 
  } 
  G_l[j*2 + a] = g_l; 
  G_r[j*2 + a] = g_r; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_y_2x2v_tensor_p1_node(int i, int j, const double *jacob_pos_l, const double *jacob_pos_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 2; ++a) { 
    f_l_quad += vst_2x2v_tensor_p1_ho_ph_x1_Cm[i*2 + a]*G_l[j*2 + a]; 
    f_r_quad += vst_2x2v_tensor_p1_ho_ph_x1_Cm[i*2 + a]*G_r[j*2 + a]; 
  } 
  f_l_quad *= 1.0/jacob_pos_l[2]; 
  f_r_quad *= 1.0/jacob_pos_r[2]; 
  const int n = i*16 + j; 
  Fhat_nodal[n] = 0.5*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void ho_lax_flux_nodal_y_2x2v_tensor_p1_prj(int unit, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (unit >= 9) return; 
  const int b = unit; 
  double t[2]; 
  for (int i = 0; i < 2; ++i) t[i] = 0.0; 
  for (int j = 0; j < 16; ++j) { 
    const double w = vst_2x2v_tensor_p1_ho_prj_x1_Vw[j*9 + b]; 
    for (int i = 0; i < 2; ++i) t[i] += w*Fhat_nodal[i*16 + j]; 
  } 
  for (int q = vst_2x2v_tensor_p1_ho_prj_x1_boff[b]; q < vst_2x2v_tensor_p1_ho_prj_x1_boff[b+1]; ++q) { 
    const int k = vst_2x2v_tensor_p1_ho_prj_x1_bks[q]; 
    const int a = vst_2x2v_tensor_p1_ho_prj_x1_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 2; ++i) g += vst_2x2v_tensor_p1_ho_prj_x1_Cw[i*2 + a]*t[i]; 
    flux[18 + k] = g; 
  } 
} 

GKYL_CU_DH double ho_lax_flux_nodal_y_2x2v_tensor_p1_cfl(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r, double alpha_max) 
{ 
  double dx11 = 2.0/dxv[1]; 
  const double jacob_pos_min = fmin(jacob_pos_l[2], jacob_pos_r[2]); 
  return 1.5*dx11*alpha_max/jacob_pos_min;
} 

GKYL_CU_DH double ho_lax_flux_nodal_y_2x2v_tensor_p1(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[32]; 
  double G_r[32]; 
  double Fhat_nodal[32]; 
  for (int item = 0; item < 32; ++item) ho_lax_flux_nodal_y_2x2v_tensor_p1_g(item, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      alpha_max = fmax(alpha_max, ho_lax_flux_nodal_y_2x2v_tensor_p1_node(i, j, jacob_pos_l, jacob_pos_r, alpha_quad[i*16 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int unit = 0; unit < 32; ++unit) ho_lax_flux_nodal_y_2x2v_tensor_p1_prj(unit, Fhat_nodal, flux); 
  return ho_lax_flux_nodal_y_2x2v_tensor_p1_cfl(dxv, jacob_pos_l, jacob_pos_r, alpha_max); 
} 
