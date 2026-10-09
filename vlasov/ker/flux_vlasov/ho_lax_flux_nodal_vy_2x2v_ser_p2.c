#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p2.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x2v_ser_p2.h> 
GKYL_CU_DH void ho_lax_flux_nodal_vy_2x2v_ser_p2_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 32) return; 
  const int a = tid/4; 
  const int j = tid - a*4; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_2x2v_ser_p2_ho_ph_v1_aoff[a]; q < vst_2x2v_ser_p2_ho_ph_v1_aoff[a+1]; ++q) { 
    const int k = vst_2x2v_ser_p2_ho_ph_v1_aks[q]; 
    g_l += vst_2x2v_ser_p2_ho_ph_v1_Wl[k*4 + j]*f_l[k]; 
    g_r += vst_2x2v_ser_p2_ho_ph_v1_Wr[k*4 + j]*f_r[k]; 
  } 
  G_l[j*8 + a] = g_l; 
  G_r[j*8 + a] = g_r; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vy_2x2v_ser_p2_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    f_l_quad += vst_2x2v_ser_p2_ho_ph_v1_Cm[i*8 + a]*G_l[j*8 + a]; 
    f_r_quad += vst_2x2v_ser_p2_ho_ph_v1_Cm[i*8 + a]*G_r[j*8 + a]; 
  } 
  const int n = i*4 + j; 
  const double jac = jacob_vel_surf_r[0]; 
  Fhat_nodal[n] = 0.5*jac*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void ho_lax_flux_nodal_vy_2x2v_ser_p2_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 20) return; 
  const int k = tid; 
  const int a = vst_2x2v_ser_p2_ho_prj_v1_kamap[k]; 
  const int b = vst_2x2v_ser_p2_ho_prj_v1_kbmap[k]; 
  double t[16]; 
  for (int i = 0; i < 16; ++i) t[i] = 0.0; 
  for (int j = 0; j < 4; ++j) { 
    const double w = vst_2x2v_ser_p2_ho_prj_v1_Vw[j*3 + b]; 
    for (int i = 0; i < 16; ++i) t[i] += w*Fhat_nodal[i*4 + j]; 
  } 
  double g = 0.0; 
  for (int i = 0; i < 16; ++i) g += vst_2x2v_ser_p2_ho_prj_v1_Cw[i*8 + a]*t[i]; 
  flux[20 + k] = g; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vy_2x2v_ser_p2_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv11 = 2.0/dxv[3]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[4], jacob_vel_surf_r[4]); 
  return 2.5*dv11*alpha_max/jacob_vel_surf_min;
} 

GKYL_CU_DH double ho_lax_flux_nodal_vy_2x2v_ser_p2(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[32]; 
  double G_r[32]; 
  double Fhat_nodal[64]; 
  for (int tid = 0; tid < 32; ++tid) ho_lax_flux_nodal_vy_2x2v_ser_p2_g(tid, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 16; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      alpha_max = fmax(alpha_max, ho_lax_flux_nodal_vy_2x2v_ser_p2_node(i, j, jacob_vel_surf_r, alpha_quad[i*4 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int tid = 0; tid < 64; ++tid) ho_lax_flux_nodal_vy_2x2v_ser_p2_prj(tid, Fhat_nodal, flux); 
  return ho_lax_flux_nodal_vy_2x2v_ser_p2_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
