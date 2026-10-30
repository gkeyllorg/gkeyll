#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_tensor_p2.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_1x3v_tensor_p2.h> 
GKYL_CU_DH void lax_flux_nodal_vx_1x3v_tensor_p2_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 48) return; 
  const int a = tid/16; 
  const int j = tid - a*16; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_1x3v_tensor_p2_ph_v0_aoff[a]; q < vst_1x3v_tensor_p2_ph_v0_aoff[a+1]; ++q) { 
    const int k = vst_1x3v_tensor_p2_ph_v0_aks[q]; 
    g_l += vst_1x3v_tensor_p2_ph_v0_Wl[k*16 + j]*f_l[k]; 
    g_r += vst_1x3v_tensor_p2_ph_v0_Wr[k*16 + j]*f_r[k]; 
  } 
  G_l[j*3 + a] = g_l; 
  G_r[j*3 + a] = g_r; 
} 

GKYL_CU_DH double lax_flux_nodal_vx_1x3v_tensor_p2_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 3; ++a) { 
    f_l_quad += vst_1x3v_tensor_p2_ph_v0_Cm[i*3 + a]*G_l[j*3 + a]; 
    f_r_quad += vst_1x3v_tensor_p2_ph_v0_Cm[i*3 + a]*G_r[j*3 + a]; 
  } 
  const int n = i*16 + j; 
  const double jac = jacob_vel_surf_r[4 + j/4]*jacob_vel_surf_r[8 + j%4]; 
  Fhat_nodal[n] = 0.5*jac*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void lax_flux_nodal_vx_1x3v_tensor_p2_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 9) return; 
  const int b = tid; 
  double t[4]; 
  for (int i = 0; i < 4; ++i) t[i] = 0.0; 
  for (int j = 0; j < 16; ++j) { 
    const double w = vst_1x3v_tensor_p2_prj_v0_Vw[j*9 + b]; 
    for (int i = 0; i < 4; ++i) t[i] += w*Fhat_nodal[i*16 + j]; 
  } 
  for (int q = vst_1x3v_tensor_p2_prj_v0_boff[b]; q < vst_1x3v_tensor_p2_prj_v0_boff[b+1]; ++q) { 
    const int k = vst_1x3v_tensor_p2_prj_v0_bks[q]; 
    const int a = vst_1x3v_tensor_p2_prj_v0_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 4; ++i) g += vst_1x3v_tensor_p2_prj_v0_Cw[i*3 + a]*t[i]; 
    flux[0 + k] = g; 
  } 
} 

GKYL_CU_DH double lax_flux_nodal_vx_1x3v_tensor_p2_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf_r[0]; 
  return 2.5*dv10*alpha_max/(-(0.11391719628198968*jacob_vel_surf_vx[3])+0.40076152031165013*jacob_vel_surf_vx[2]-0.8136324494869249*jacob_vel_surf_vx[1]+1.5267881254572662*jacob_vel_surf_vx[0]);
} 

GKYL_CU_DH double lax_flux_nodal_vx_1x3v_tensor_p2(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[48]; 
  double G_r[48]; 
  double Fhat_nodal[64]; 
  for (int tid = 0; tid < 48; ++tid) lax_flux_nodal_vx_1x3v_tensor_p2_g(tid, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      alpha_max = fmax(alpha_max, lax_flux_nodal_vx_1x3v_tensor_p2_node(i, j, jacob_vel_surf_r, alpha_quad[i*16 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int tid = 0; tid < 64; ++tid) lax_flux_nodal_vx_1x3v_tensor_p2_prj(tid, Fhat_nodal, flux); 
  return lax_flux_nodal_vx_1x3v_tensor_p2_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
