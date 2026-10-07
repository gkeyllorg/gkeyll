#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x2v_ser_p1.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x2v_ser_p1.h> 
GKYL_CU_DH void lax_flux_nodal_vy_1x2v_ser_p1_g(int item, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (item >= 4) return; 
  const int j = item/2; 
  const int a = item - j*2; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_1x2v_ser_p1_ph_v1_aoff[a]; q < vst_1x2v_ser_p1_ph_v1_aoff[a+1]; ++q) { 
    const int k = vst_1x2v_ser_p1_ph_v1_aks[q]; 
    g_l += vst_1x2v_ser_p1_ph_v1_V[j*4 + vst_1x2v_ser_p1_ph_v1_vlmap[k]]*(vst_1x2v_ser_p1_ph_v1_coefl[k]*f_l[k]); 
    g_r += vst_1x2v_ser_p1_ph_v1_V[j*4 + vst_1x2v_ser_p1_ph_v1_vrmap[k]]*(vst_1x2v_ser_p1_ph_v1_coefr[k]*f_r[k]); 
  } 
  G_l[item] = g_l; 
  G_r[item] = g_r; 
} 

GKYL_CU_DH double lax_flux_nodal_vy_1x2v_ser_p1_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 2; ++a) { 
    f_l_quad += vst_1x2v_ser_p1_ph_v1_Cm[i*2 + a]*G_l[j*2 + a]; 
    f_r_quad += vst_1x2v_ser_p1_ph_v1_Cm[i*2 + a]*G_r[j*2 + a]; 
  } 
  const int n = i*2 + j; 
  const double jac = jacob_vel_surf_r[0]; 
  Fhat_nodal[n] = 0.5*jac*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void lax_flux_nodal_vy_1x2v_ser_p1_prj(int k, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (k >= 4) return; 
  const int a = vst_1x2v_ser_p1_prj_v1_kamap[k]; 
  const int b = vst_1x2v_ser_p1_prj_v1_kbmap[k]; 
  double t[2]; 
  for (int i = 0; i < 2; ++i) t[i] = 0.0; 
  for (int j = 0; j < 2; ++j) { 
    const double w = vst_1x2v_ser_p1_prj_v1_Vw[j*2 + b]; 
    for (int i = 0; i < 2; ++i) t[i] += w*Fhat_nodal[i*2 + j]; 
  } 
  double g = 0.0; 
  for (int i = 0; i < 2; ++i) g += vst_1x2v_ser_p1_prj_v1_Cw[i*2 + a]*t[i]; 
  flux[4 + k] = g; 
} 

GKYL_CU_DH double lax_flux_nodal_vy_1x2v_ser_p1_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv11 = 2.0/dxv[2]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[3], jacob_vel_surf_r[3]); 
  return 1.5*dv11*alpha_max/jacob_vel_surf_min;
} 

GKYL_CU_DH double lax_flux_nodal_vy_1x2v_ser_p1(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[4]; 
  double G_r[4]; 
  double Fhat_nodal[4]; 
  for (int item = 0; item < 4; ++item) lax_flux_nodal_vy_1x2v_ser_p1_g(item, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 2; ++j) { 
      alpha_max = fmax(alpha_max, lax_flux_nodal_vy_1x2v_ser_p1_node(i, j, jacob_vel_surf_r, alpha_quad[i*2 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int k = 0; k < 4; ++k) lax_flux_nodal_vy_1x2v_ser_p1_prj(k, Fhat_nodal, flux); 
  return lax_flux_nodal_vy_1x2v_ser_p1_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
