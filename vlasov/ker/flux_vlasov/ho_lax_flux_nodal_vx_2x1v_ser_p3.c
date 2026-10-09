#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x1v_ser_p3.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x1v_ser_p3.h> 
GKYL_CU_DH void ho_lax_flux_nodal_vx_2x1v_ser_p3_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 12) return; 
  const int a = tid/1; 
  const int j = tid - a*1; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_2x1v_ser_p3_ho_ph_v0_aoff[a]; q < vst_2x1v_ser_p3_ho_ph_v0_aoff[a+1]; ++q) { 
    const int k = vst_2x1v_ser_p3_ho_ph_v0_aks[q]; 
    g_l += vst_2x1v_ser_p3_ho_ph_v0_Wl[k*1 + j]*f_l[k]; 
    g_r += vst_2x1v_ser_p3_ho_ph_v0_Wr[k*1 + j]*f_r[k]; 
  } 
  G_l[j*12 + a] = g_l; 
  G_r[j*12 + a] = g_r; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vx_2x1v_ser_p3_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 12; ++a) { 
    f_l_quad += vst_2x1v_ser_p3_ho_ph_v0_Cm[i*12 + a]*G_l[j*12 + a]; 
    f_r_quad += vst_2x1v_ser_p3_ho_ph_v0_Cm[i*12 + a]*G_r[j*12 + a]; 
  } 
  const int n = i*1 + j; 
  Fhat_nodal[n] = 0.5*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void ho_lax_flux_nodal_vx_2x1v_ser_p3_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 12) return; 
  const int k = tid; 
  const int a = vst_2x1v_ser_p3_ho_prj_v0_kamap[k]; 
  const int b = vst_2x1v_ser_p3_ho_prj_v0_kbmap[k]; 
  double t[25]; 
  for (int i = 0; i < 25; ++i) t[i] = 0.0; 
  for (int j = 0; j < 1; ++j) { 
    const double w = vst_2x1v_ser_p3_ho_prj_v0_Vw[j*1 + b]; 
    for (int i = 0; i < 25; ++i) t[i] += w*Fhat_nodal[i*1 + j]; 
  } 
  double g = 0.0; 
  for (int i = 0; i < 25; ++i) g += vst_2x1v_ser_p3_ho_prj_v0_Cw[i*12 + a]*t[i]; 
  flux[0 + k] = g; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vx_2x1v_ser_p3_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv10 = 2.0/dxv[2]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[0], jacob_vel_surf_r[0]); 
  return 3.5*dv10*alpha_max/jacob_vel_surf_min;
} 

GKYL_CU_DH double ho_lax_flux_nodal_vx_2x1v_ser_p3(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[12]; 
  double G_r[12]; 
  double Fhat_nodal[25]; 
  for (int tid = 0; tid < 12; ++tid) ho_lax_flux_nodal_vx_2x1v_ser_p3_g(tid, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 25; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      alpha_max = fmax(alpha_max, ho_lax_flux_nodal_vx_2x1v_ser_p3_node(i, j, jacob_vel_surf_r, alpha_quad[i*1 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int tid = 0; tid < 25; ++tid) ho_lax_flux_nodal_vx_2x1v_ser_p3_prj(tid, Fhat_nodal, flux); 
  return ho_lax_flux_nodal_vx_2x1v_ser_p3_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
