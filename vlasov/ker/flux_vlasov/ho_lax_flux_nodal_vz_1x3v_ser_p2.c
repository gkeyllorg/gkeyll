#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH void ho_lax_flux_nodal_vz_1x3v_ser_p2_g(int item, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (item >= 48) return; 
  const int j = item/3; 
  const int a = item - j*3; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_1x3v_ser_p2_ho_ph_v2_aoff[a]; q < vst_1x3v_ser_p2_ho_ph_v2_aoff[a+1]; ++q) { 
    const int k = vst_1x3v_ser_p2_ho_ph_v2_aks[q]; 
    g_l += vst_1x3v_ser_p2_ho_ph_v2_V[j*20 + vst_1x3v_ser_p2_ho_ph_v2_vlmap[k]]*(vst_1x3v_ser_p2_ho_ph_v2_coefl[k]*f_l[k]); 
    g_r += vst_1x3v_ser_p2_ho_ph_v2_V[j*20 + vst_1x3v_ser_p2_ho_ph_v2_vrmap[k]]*(vst_1x3v_ser_p2_ho_ph_v2_coefr[k]*f_r[k]); 
  } 
  G_l[item] = g_l; 
  G_r[item] = g_r; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vz_1x3v_ser_p2_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 3; ++a) { 
    f_l_quad += vst_1x3v_ser_p2_ho_ph_v2_Cm[i*3 + a]*G_l[j*3 + a]; 
    f_r_quad += vst_1x3v_ser_p2_ho_ph_v2_Cm[i*3 + a]*G_r[j*3 + a]; 
  } 
  const int n = i*16 + j; 
  const double jac = jacob_vel_surf_r[0]*jacob_vel_surf_r[4]; 
  Fhat_nodal[n] = 0.5*jac*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void ho_lax_flux_nodal_vz_1x3v_ser_p2_prj(int k, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (k >= 20) return; 
  const int a = vst_1x3v_ser_p2_ho_prj_v2_kamap[k]; 
  const int b = vst_1x3v_ser_p2_ho_prj_v2_kbmap[k]; 
  double t[4]; 
  for (int i = 0; i < 4; ++i) t[i] = 0.0; 
  for (int j = 0; j < 16; ++j) { 
    const double w = vst_1x3v_ser_p2_ho_prj_v2_Vw[j*8 + b]; 
    for (int i = 0; i < 4; ++i) t[i] += w*Fhat_nodal[i*16 + j]; 
  } 
  double g = 0.0; 
  for (int i = 0; i < 4; ++i) g += vst_1x3v_ser_p2_ho_prj_v2_Cw[i*3 + a]*t[i]; 
  flux[40 + k] = g; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vz_1x3v_ser_p2_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv12 = 2.0/dxv[3]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[8], jacob_vel_surf_r[8]); 
  return 2.5*dv12*alpha_max/jacob_vel_surf_min;
} 

GKYL_CU_DH double ho_lax_flux_nodal_vz_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[48]; 
  double G_r[48]; 
  double Fhat_nodal[64]; 
  for (int item = 0; item < 48; ++item) ho_lax_flux_nodal_vz_1x3v_ser_p2_g(item, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      alpha_max = fmax(alpha_max, ho_lax_flux_nodal_vz_1x3v_ser_p2_node(i, j, jacob_vel_surf_r, alpha_quad[i*16 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int k = 0; k < 20; ++k) ho_lax_flux_nodal_vz_1x3v_ser_p2_prj(k, Fhat_nodal, flux); 
  return ho_lax_flux_nodal_vz_1x3v_ser_p2_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
