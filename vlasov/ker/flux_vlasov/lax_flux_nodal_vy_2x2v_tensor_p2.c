#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p2.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x2v_tensor_p2.h> 
GKYL_CU_DH void lax_flux_nodal_vy_2x2v_tensor_p2_g(int item, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (item >= 36) return; 
  const int j = item/9; 
  const int a = item - j*9; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_2x2v_tensor_p2_ph_v1_aoff[a]; q < vst_2x2v_tensor_p2_ph_v1_aoff[a+1]; ++q) { 
    const int k = vst_2x2v_tensor_p2_ph_v1_aks[q]; 
    g_l += vst_2x2v_tensor_p2_ph_v1_V[j*9 + vst_2x2v_tensor_p2_ph_v1_vlmap[k]]*(vst_2x2v_tensor_p2_ph_v1_coefl[k]*f_l[k]); 
    g_r += vst_2x2v_tensor_p2_ph_v1_V[j*9 + vst_2x2v_tensor_p2_ph_v1_vrmap[k]]*(vst_2x2v_tensor_p2_ph_v1_coefr[k]*f_r[k]); 
  } 
  G_l[item] = g_l; 
  G_r[item] = g_r; 
} 

GKYL_CU_DH double lax_flux_nodal_vy_2x2v_tensor_p2_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 9; ++a) { 
    f_l_quad += vst_2x2v_tensor_p2_ph_v1_Cm[i*9 + a]*G_l[j*9 + a]; 
    f_r_quad += vst_2x2v_tensor_p2_ph_v1_Cm[i*9 + a]*G_r[j*9 + a]; 
  } 
  const int n = i*4 + j; 
  const double jac = jacob_vel_surf_r[0 + j]; 
  Fhat_nodal[n] = 0.5*jac*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void lax_flux_nodal_vy_2x2v_tensor_p2_prj(int k, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (k >= 27) return; 
  const int a = vst_2x2v_tensor_p2_prj_v1_kamap[k]; 
  const int b = vst_2x2v_tensor_p2_prj_v1_kbmap[k]; 
  double t[16]; 
  for (int i = 0; i < 16; ++i) t[i] = 0.0; 
  for (int j = 0; j < 4; ++j) { 
    const double w = vst_2x2v_tensor_p2_prj_v1_Vw[j*3 + b]; 
    for (int i = 0; i < 16; ++i) t[i] += w*Fhat_nodal[i*4 + j]; 
  } 
  double g = 0.0; 
  for (int i = 0; i < 16; ++i) g += vst_2x2v_tensor_p2_prj_v1_Cw[i*9 + a]*t[i]; 
  flux[27 + k] = g; 
} 

GKYL_CU_DH double lax_flux_nodal_vy_2x2v_tensor_p2_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vel_surf_vy = &jacob_vel_surf_r[4]; 
  return 2.5*dv11*alpha_max/(-(0.11391719628198968*jacob_vel_surf_vy[3])+0.40076152031165013*jacob_vel_surf_vy[2]-0.8136324494869249*jacob_vel_surf_vy[1]+1.5267881254572662*jacob_vel_surf_vy[0]);
} 

GKYL_CU_DH double lax_flux_nodal_vy_2x2v_tensor_p2(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[36]; 
  double G_r[36]; 
  double Fhat_nodal[64]; 
  for (int item = 0; item < 36; ++item) lax_flux_nodal_vy_2x2v_tensor_p2_g(item, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      alpha_max = fmax(alpha_max, lax_flux_nodal_vy_2x2v_tensor_p2_node(i, j, jacob_vel_surf_r, alpha_quad[i*4 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int k = 0; k < 27; ++k) lax_flux_nodal_vy_2x2v_tensor_p2_prj(k, Fhat_nodal, flux); 
  return lax_flux_nodal_vy_2x2v_tensor_p2_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
