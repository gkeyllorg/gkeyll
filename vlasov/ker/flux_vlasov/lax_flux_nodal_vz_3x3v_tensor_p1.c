#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH void lax_flux_nodal_vz_3x3v_tensor_p1_g(int item, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (item >= 72) return; 
  const int a = item/9; 
  const int j = item - a*9; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_3x3v_tensor_p1_ph_v2_aoff[a]; q < vst_3x3v_tensor_p1_ph_v2_aoff[a+1]; ++q) { 
    const int k = vst_3x3v_tensor_p1_ph_v2_aks[q]; 
    g_l += vst_3x3v_tensor_p1_ph_v2_Wl[k*9 + j]*f_l[k]; 
    g_r += vst_3x3v_tensor_p1_ph_v2_Wr[k*9 + j]*f_r[k]; 
  } 
  G_l[j*8 + a] = g_l; 
  G_r[j*8 + a] = g_r; 
} 

GKYL_CU_DH double lax_flux_nodal_vz_3x3v_tensor_p1_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    f_l_quad += vst_3x3v_tensor_p1_ph_v2_Cm[i*8 + a]*G_l[j*8 + a]; 
    f_r_quad += vst_3x3v_tensor_p1_ph_v2_Cm[i*8 + a]*G_r[j*8 + a]; 
  } 
  const int n = i*9 + j; 
  const double jac = jacob_vel_surf_r[0 + j/3]*jacob_vel_surf_r[3 + j%3]; 
  Fhat_nodal[n] = 0.5*jac*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void lax_flux_nodal_vz_3x3v_tensor_p1_prj(int unit, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (unit >= 72) return; 
  const int k = unit; 
  const int a = vst_3x3v_tensor_p1_prj_v2_kamap[k]; 
  const int b = vst_3x3v_tensor_p1_prj_v2_kbmap[k]; 
  double t[8]; 
  for (int i = 0; i < 8; ++i) t[i] = 0.0; 
  for (int j = 0; j < 9; ++j) { 
    const double w = vst_3x3v_tensor_p1_prj_v2_Vw[j*9 + b]; 
    for (int i = 0; i < 8; ++i) t[i] += w*Fhat_nodal[i*9 + j]; 
  } 
  double g = 0.0; 
  for (int i = 0; i < 8; ++i) g += vst_3x3v_tensor_p1_prj_v2_Cw[i*8 + a]*t[i]; 
  flux[144 + k] = g; 
} 

GKYL_CU_DH double lax_flux_nodal_vz_3x3v_tensor_p1_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv12 = 2.0/dxv[5]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf_r[6]; 
  return 1.5*dv12*alpha_max/(0.18783610896543046*jacob_vel_surf_vz[2]-0.6666666666666666*jacob_vel_surf_vz[1]+1.4788305577012362*jacob_vel_surf_vz[0]);
} 

GKYL_CU_DH double lax_flux_nodal_vz_3x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[72]; 
  double G_r[72]; 
  double Fhat_nodal[72]; 
  for (int item = 0; item < 72; ++item) lax_flux_nodal_vz_3x3v_tensor_p1_g(item, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      alpha_max = fmax(alpha_max, lax_flux_nodal_vz_3x3v_tensor_p1_node(i, j, jacob_vel_surf_r, alpha_quad[i*9 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int unit = 0; unit < 72; ++unit) lax_flux_nodal_vz_3x3v_tensor_p1_prj(unit, Fhat_nodal, flux); 
  return lax_flux_nodal_vz_3x3v_tensor_p1_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
