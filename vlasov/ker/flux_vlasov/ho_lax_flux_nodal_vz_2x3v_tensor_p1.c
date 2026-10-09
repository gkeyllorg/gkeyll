#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x3v_tensor_p1.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x3v_tensor_p1.h> 
GKYL_CU_DH void ho_lax_flux_nodal_vz_2x3v_tensor_p1_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 64) return; 
  const int a = tid/16; 
  const int j = tid - a*16; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_2x3v_tensor_p1_ho_ph_v2_aoff[a]; q < vst_2x3v_tensor_p1_ho_ph_v2_aoff[a+1]; ++q) { 
    const int k = vst_2x3v_tensor_p1_ho_ph_v2_aks[q]; 
    g_l += vst_2x3v_tensor_p1_ho_ph_v2_Wl[k*16 + j]*f_l[k]; 
    g_r += vst_2x3v_tensor_p1_ho_ph_v2_Wr[k*16 + j]*f_r[k]; 
  } 
  G_l[j*4 + a] = g_l; 
  G_r[j*4 + a] = g_r; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vz_2x3v_tensor_p1_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    f_l_quad += vst_2x3v_tensor_p1_ho_ph_v2_Cm[i*4 + a]*G_l[j*4 + a]; 
    f_r_quad += vst_2x3v_tensor_p1_ho_ph_v2_Cm[i*4 + a]*G_r[j*4 + a]; 
  } 
  const int n = i*16 + j; 
  const double jac = jacob_vel_surf_r[0 + j/4]*jacob_vel_surf_r[4 + j%4]; 
  Fhat_nodal[n] = 0.5*jac*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void ho_lax_flux_nodal_vz_2x3v_tensor_p1_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 9) return; 
  const int b = tid; 
  double t[4]; 
  for (int i = 0; i < 4; ++i) t[i] = 0.0; 
  for (int j = 0; j < 16; ++j) { 
    const double w = vst_2x3v_tensor_p1_ho_prj_v2_Vw[j*9 + b]; 
    for (int i = 0; i < 4; ++i) t[i] += w*Fhat_nodal[i*16 + j]; 
  } 
  for (int q = vst_2x3v_tensor_p1_ho_prj_v2_boff[b]; q < vst_2x3v_tensor_p1_ho_prj_v2_boff[b+1]; ++q) { 
    const int k = vst_2x3v_tensor_p1_ho_prj_v2_bks[q]; 
    const int a = vst_2x3v_tensor_p1_ho_prj_v2_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 4; ++i) g += vst_2x3v_tensor_p1_ho_prj_v2_Cw[i*4 + a]*t[i]; 
    flux[72 + k] = g; 
  } 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vz_2x3v_tensor_p1_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv12 = 2.0/dxv[4]; 
  const double *jacob_vel_surf_vz = &jacob_vel_surf_r[8]; 
  return 1.5*dv12*alpha_max/(-(0.11391719628198968*jacob_vel_surf_vz[3])+0.40076152031165013*jacob_vel_surf_vz[2]-0.8136324494869249*jacob_vel_surf_vz[1]+1.5267881254572662*jacob_vel_surf_vz[0]);
} 

GKYL_CU_DH double ho_lax_flux_nodal_vz_2x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[64]; 
  double G_r[64]; 
  double Fhat_nodal[64]; 
  for (int tid = 0; tid < 64; ++tid) ho_lax_flux_nodal_vz_2x3v_tensor_p1_g(tid, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      alpha_max = fmax(alpha_max, ho_lax_flux_nodal_vz_2x3v_tensor_p1_node(i, j, jacob_vel_surf_r, alpha_quad[i*16 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int tid = 0; tid < 64; ++tid) ho_lax_flux_nodal_vz_2x3v_tensor_p1_prj(tid, Fhat_nodal, flux); 
  return ho_lax_flux_nodal_vz_2x3v_tensor_p1_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
