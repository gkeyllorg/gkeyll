#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH void lax_flux_nodal_vx_2x3v_tensor_p1_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 36) return; 
  const int a = tid/9; 
  const int j = tid - a*9; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_2x3v_tensor_p1_ph_v0_aoff[a]; q < vst_2x3v_tensor_p1_ph_v0_aoff[a+1]; ++q) { 
    const int k = vst_2x3v_tensor_p1_ph_v0_aks[q]; 
    g_l += vst_2x3v_tensor_p1_ph_v0_Wl[k*9 + j]*f_l[k]; 
    g_r += vst_2x3v_tensor_p1_ph_v0_Wr[k*9 + j]*f_r[k]; 
  } 
  G_l[j*4 + a] = g_l; 
  G_r[j*4 + a] = g_r; 
} 

GKYL_CU_DH double lax_flux_nodal_vx_2x3v_tensor_p1_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    f_l_quad += vst_2x3v_tensor_p1_ph_v0_Cm[i*4 + a]*G_l[j*4 + a]; 
    f_r_quad += vst_2x3v_tensor_p1_ph_v0_Cm[i*4 + a]*G_r[j*4 + a]; 
  } 
  const int n = i*9 + j; 
  const double jac = jacob_vel_surf_r[3 + j/3]*jacob_vel_surf_r[6 + j%3]; 
  Fhat_nodal[n] = 0.5*jac*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void lax_flux_nodal_vx_2x3v_tensor_p1_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 9) return; 
  const int b = tid; 
  double t[4]; 
  for (int i = 0; i < 4; ++i) t[i] = 0.0; 
  for (int j = 0; j < 9; ++j) { 
    const double w = vst_2x3v_tensor_p1_prj_v0_Vw[j*9 + b]; 
    for (int i = 0; i < 4; ++i) t[i] += w*Fhat_nodal[i*9 + j]; 
  } 
  for (int q = vst_2x3v_tensor_p1_prj_v0_boff[b]; q < vst_2x3v_tensor_p1_prj_v0_boff[b+1]; ++q) { 
    const int k = vst_2x3v_tensor_p1_prj_v0_bks[q]; 
    const int a = vst_2x3v_tensor_p1_prj_v0_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 4; ++i) g += vst_2x3v_tensor_p1_prj_v0_Cw[i*4 + a]*t[i]; 
    flux[0 + k] = g; 
  } 
} 

GKYL_CU_DH double lax_flux_nodal_vx_2x3v_tensor_p1_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf_r[0]; 
  return 1.5*dv10*alpha_max/(0.18783610896543046*jacob_vel_surf_vx[2]-0.6666666666666666*jacob_vel_surf_vx[1]+1.4788305577012362*jacob_vel_surf_vx[0]);
} 

GKYL_CU_DH double lax_flux_nodal_vx_2x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[36]; 
  double G_r[36]; 
  double Fhat_nodal[36]; 
  for (int tid = 0; tid < 36; ++tid) lax_flux_nodal_vx_2x3v_tensor_p1_g(tid, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      alpha_max = fmax(alpha_max, lax_flux_nodal_vx_2x3v_tensor_p1_node(i, j, jacob_vel_surf_r, alpha_quad[i*9 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int tid = 0; tid < 36; ++tid) lax_flux_nodal_vx_2x3v_tensor_p1_prj(tid, Fhat_nodal, flux); 
  return lax_flux_nodal_vx_2x3v_tensor_p1_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
