#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_ser_p2.h> 
#include <gkyl_vlasov_surf_nod2mod_tables_2x2v_ser_p2.h> 
GKYL_CU_DH void ho_lax_flux_nodal_y_2x2v_ser_p2_g(int tid, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  if (tid >= 48) return; 
  const int a = tid/16; 
  const int j = tid - a*16; 
  double g_l = 0.0; 
  double g_r = 0.0; 
  for (int q = vst_2x2v_ser_p2_ho_ph_x1_aoff[a]; q < vst_2x2v_ser_p2_ho_ph_x1_aoff[a+1]; ++q) { 
    const int k = vst_2x2v_ser_p2_ho_ph_x1_aks[q]; 
    g_l += vst_2x2v_ser_p2_ho_ph_x1_Wl[k*16 + j]*f_l[k]; 
    g_r += vst_2x2v_ser_p2_ho_ph_x1_Wr[k*16 + j]*f_r[k]; 
  } 
  G_l[j*3 + a] = g_l; 
  G_r[j*3 + a] = g_r; 
} 

GKYL_CU_DH double ho_lax_flux_nodal_y_2x2v_ser_p2_node(int i, int j, const double *jacob_pos_l, const double *jacob_pos_r,
  double alpha, const double *G_l, const double *G_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 3; ++a) { 
    f_l_quad += vst_2x2v_ser_p2_ho_ph_x1_Cm[i*3 + a]*G_l[j*3 + a]; 
    f_r_quad += vst_2x2v_ser_p2_ho_ph_x1_Cm[i*3 + a]*G_r[j*3 + a]; 
  } 
  f_l_quad *= 1.0/jacob_pos_l[3]; 
  f_r_quad *= 1.0/jacob_pos_r[3]; 
  const int n = i*16 + j; 
  Fhat_nodal[n] = 0.5*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH void ho_lax_flux_nodal_y_2x2v_ser_p2_prj(int tid, const double *Fhat_nodal, double* GKYL_RESTRICT flux) 
{ 
  if (tid >= 8) return; 
  const int b = tid; 
  double t[4]; 
  for (int i = 0; i < 4; ++i) t[i] = 0.0; 
  for (int j = 0; j < 16; ++j) { 
    const double w = vst_2x2v_ser_p2_ho_prj_x1_Vw[j*8 + b]; 
    for (int i = 0; i < 4; ++i) t[i] += w*Fhat_nodal[i*16 + j]; 
  } 
  for (int q = vst_2x2v_ser_p2_ho_prj_x1_boff[b]; q < vst_2x2v_ser_p2_ho_prj_x1_boff[b+1]; ++q) { 
    const int k = vst_2x2v_ser_p2_ho_prj_x1_bks[q]; 
    const int a = vst_2x2v_ser_p2_ho_prj_x1_kamap[k]; 
    double g = 0.0; 
    for (int i = 0; i < 4; ++i) g += vst_2x2v_ser_p2_ho_prj_x1_Cw[i*3 + a]*t[i]; 
    flux[20 + k] = g; 
  } 
} 

GKYL_CU_DH double ho_lax_flux_nodal_y_2x2v_ser_p2_cfl(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r, double alpha_max) 
{ 
  double dx11 = 2.0/dxv[1]; 
  const double jacob_pos_min = fmin(jacob_pos_l[3], jacob_pos_r[3]); 
  return 2.5*dx11*alpha_max/jacob_pos_min;
} 

GKYL_CU_DH double ho_lax_flux_nodal_y_2x2v_ser_p2(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT flux) 
{ 
  double G_l[48]; 
  double G_r[48]; 
  double Fhat_nodal[64]; 
  for (int tid = 0; tid < 48; ++tid) ho_lax_flux_nodal_y_2x2v_ser_p2_g(tid, f_l, f_r, G_l, G_r); 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      alpha_max = fmax(alpha_max, ho_lax_flux_nodal_y_2x2v_ser_p2_node(i, j, jacob_pos_l, jacob_pos_r, alpha_quad[i*16 + j], G_l, G_r, Fhat_nodal)); 
    } 
  } 
  for (int tid = 0; tid < 64; ++tid) ho_lax_flux_nodal_y_2x2v_ser_p2_prj(tid, Fhat_nodal, flux); 
  return ho_lax_flux_nodal_y_2x2v_ser_p2_cfl(dxv, jacob_pos_l, jacob_pos_r, alpha_max); 
} 
