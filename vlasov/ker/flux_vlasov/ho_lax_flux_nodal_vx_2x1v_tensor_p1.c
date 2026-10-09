#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x1v_tensor_p1.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x1v_tensor_p1.h> 
GKYL_CU_DH void ho_lax_flux_nodal_vx_2x1v_tensor_p1_g(int item, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (item >= 4) return; 
  const int a = item/1; 
  const int j = item - a*1; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_2x1v_tensor_p1_ho_ph_v0_aoff[a]; q < vst_2x1v_tensor_p1_ho_ph_v0_aoff[a+1]; ++q) { 
    const int k = vst_2x1v_tensor_p1_ho_ph_v0_aks[q]; 
    g_l += vst_2x1v_tensor_p1_ho_ph_v0_Wl[k*1 + j]*f_l[k]; 
    g_r += vst_2x1v_tensor_p1_ho_ph_v0_Wr[k*1 + j]*f_r[k]; 
  } 
  G_l[j*4 + a] = g_l; 
  G_r[j*4 + a] = g_r; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vx_2x1v_tensor_p1_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    f_l_quad += vst_2x1v_tensor_p1_ho_ph_v0_Cm[i*4 + a]*G_l[j*4 + a]; 
    f_r_quad += vst_2x1v_tensor_p1_ho_ph_v0_Cm[i*4 + a]*G_r[j*4 + a]; 
  } 
  const int n = i*1 + j; 
  Fhat_nodal[n] = 0.5*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void ho_lax_flux_nodal_vx_2x1v_tensor_p1_prj(int unit, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (unit >= 4) return; 
  const int k = unit; 
  const int a = vst_2x1v_tensor_p1_ho_prj_v0_kamap[k]; 
  const int b = vst_2x1v_tensor_p1_ho_prj_v0_kbmap[k]; 
  double t[4]; 
  for (int i = 0; i < 4; ++i) t[i] = 0.0; 
  for (int j = 0; j < 1; ++j) { 
    const double w = vst_2x1v_tensor_p1_ho_prj_v0_Vw[j*1 + b]; 
    for (int i = 0; i < 4; ++i) t[i] += w*Fhat_nodal[i*1 + j]; 
  } 
  double g = 0.0; 
  for (int i = 0; i < 4; ++i) g += vst_2x1v_tensor_p1_ho_prj_v0_Cw[i*4 + a]*t[i]; 
  flux[0 + k] = g; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vx_2x1v_tensor_p1_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf_r[0]; 
  return 1.5*dv10*alpha_max/(-(0.11391719628198968*jacob_vel_surf_vx[3])+0.40076152031165013*jacob_vel_surf_vx[2]-0.8136324494869249*jacob_vel_surf_vx[1]+1.5267881254572662*jacob_vel_surf_vx[0]);
} 

GKYL_CU_DH double ho_lax_flux_nodal_vx_2x1v_tensor_p1(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[4]; 
  double G_r[4]; 
  double Fhat_nodal[4]; 
  for (int item = 0; item < 4; ++item) ho_lax_flux_nodal_vx_2x1v_tensor_p1_g(item, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 1; ++j) { 
      alpha_max = fmax(alpha_max, ho_lax_flux_nodal_vx_2x1v_tensor_p1_node(i, j, jacob_vel_surf_r, alpha_quad[i*1 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int unit = 0; unit < 4; ++unit) ho_lax_flux_nodal_vx_2x1v_tensor_p1_prj(unit, Fhat_nodal, flux); 
  return ho_lax_flux_nodal_vx_2x1v_tensor_p1_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
