#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH double ho_lax_flux_nodal_x_2x2v_tensor_p1_node(int i, int j, const double *jacob_pos_l, const double *jacob_pos_r,
  double alpha, const double *f_l, const double *f_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double G_l[2]; 
  double G_r[2]; 
  for (int a = 0; a < 2; ++a) { G_l[a] = 0.0; G_r[a] = 0.0; } 
  for (int k = 0; k < 36; ++k) { 
    const int a = vst_2x2v_tensor_p1_ho_ph_x0_cmap[k]; 
    G_l[a] += vst_2x2v_tensor_p1_ho_ph_x0_V[j*18 + vst_2x2v_tensor_p1_ho_ph_x0_vlmap[k]]*(vst_2x2v_tensor_p1_ho_ph_x0_coefl[k]*f_l[k]); 
    G_r[a] += vst_2x2v_tensor_p1_ho_ph_x0_V[j*18 + vst_2x2v_tensor_p1_ho_ph_x0_vrmap[k]]*(vst_2x2v_tensor_p1_ho_ph_x0_coefr[k]*f_r[k]); 
  } 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 2; ++a) { 
    f_l_quad += vst_2x2v_tensor_p1_ho_ph_x0_Cm[i*2 + a]*G_l[a]; 
    f_r_quad += vst_2x2v_tensor_p1_ho_ph_x0_Cm[i*2 + a]*G_r[a]; 
  } 
  f_l_quad *= 1.0/jacob_pos_l[0]; 
  f_r_quad *= 1.0/jacob_pos_r[0]; 
  const int n = i*16 + j; 
  Fhat_nodal[0 + n] = 0.5*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH double ho_lax_flux_nodal_x_2x2v_tensor_p1_cfl(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r, double alpha_max) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double jacob_pos_min = fmin(jacob_pos_l[0], jacob_pos_r[0]); 
  return 1.5*dx10*alpha_max/jacob_pos_min;
} 

GKYL_CU_DH double ho_lax_flux_nodal_x_2x2v_tensor_p1(const double *dxv, const double *jacob_pos_l, const double *jacob_pos_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double alpha_max = 0.0; 
  for (int i = 0; i < 2; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      alpha_max = fmax(alpha_max, ho_lax_flux_nodal_x_2x2v_tensor_p1_node(i, j, jacob_pos_l, jacob_pos_r, alpha_quad[i*16 + j], f_l, f_r, Fhat_nodal)); 
    } 
  } 
  return ho_lax_flux_nodal_x_2x2v_tensor_p1_cfl(dxv, jacob_pos_l, jacob_pos_r, alpha_max); 
} 
