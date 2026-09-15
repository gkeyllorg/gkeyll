#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_ser_p1.h> 
GKYL_CU_DH double lax_flux_nodal_y_3x3v_ser_p1_node(int i, int j, const double *jacob_pos_l, const double *jacob_pos_r,
  double alpha, const double *f_l, const double *f_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double G_l[4]; 
  double G_r[4]; 
  for (int a = 0; a < 4; ++a) { G_l[a] = 0.0; G_r[a] = 0.0; } 
  for (int k = 0; k < 64; ++k) { 
    const int a = vst_3x3v_ser_p1_ph_x1_cmap[k]; 
    G_l[a] += vst_3x3v_ser_p1_ph_x1_V[j*16 + vst_3x3v_ser_p1_ph_x1_vlmap[k]]*(vst_3x3v_ser_p1_ph_x1_coefl[k]*f_l[k]); 
    G_r[a] += vst_3x3v_ser_p1_ph_x1_V[j*16 + vst_3x3v_ser_p1_ph_x1_vrmap[k]]*(vst_3x3v_ser_p1_ph_x1_coefr[k]*f_r[k]); 
  } 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    f_l_quad += vst_3x3v_ser_p1_ph_x1_Cm[i*4 + a]*G_l[a]; 
    f_r_quad += vst_3x3v_ser_p1_ph_x1_Cm[i*4 + a]*G_r[a]; 
  } 
  f_l_quad *= 1.0/jacob_pos_l[2]; 
  f_r_quad *= 1.0/jacob_pos_r[2]; 
  const int n = i*8 + j; 
  Fhat_nodal[32 + n] = 0.5*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH double lax_flux_nodal_y_3x3v_ser_p1_cfl(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r, double alpha_max) 
{ 
  double dx11 = 2.0/dxv[1]; 
  const double jacob_pos_min = fmin(jacob_pos_l[2], jacob_pos_r[2]); 
  return 1.5*dx11*alpha_max/jacob_pos_min;
} 

GKYL_CU_DH double lax_flux_nodal_y_3x3v_ser_p1(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 8; ++j) { 
      alpha_max = fmax(alpha_max, lax_flux_nodal_y_3x3v_ser_p1_node(i, j, jacob_pos_l, jacob_pos_r, alpha_quad[i*8 + j], f_l, f_r, Fhat_nodal)); 
    } 
  } 
  return lax_flux_nodal_y_3x3v_ser_p1_cfl(dxv, jacob_pos_l, jacob_pos_r, alpha_max); 
} 
