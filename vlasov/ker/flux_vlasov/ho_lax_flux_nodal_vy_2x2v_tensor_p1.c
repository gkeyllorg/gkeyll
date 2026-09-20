#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_2x2v_tensor_p1.h> 
GKYL_CU_DH double ho_lax_flux_nodal_vy_2x2v_tensor_p1_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *f_l, const double *f_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double G_l[4]; 
  double G_r[4]; 
  for (int a = 0; a < 4; ++a) { G_l[a] = 0.0; G_r[a] = 0.0; } 
  for (int k = 0; k < 36; ++k) { 
    const int a = vst_2x2v_tensor_p1_ho_ph_v1_cmap[k]; 
    G_l[a] += vst_2x2v_tensor_p1_ho_ph_v1_V[j*9 + vst_2x2v_tensor_p1_ho_ph_v1_vlmap[k]]*(vst_2x2v_tensor_p1_ho_ph_v1_coefl[k]*f_l[k]); 
    G_r[a] += vst_2x2v_tensor_p1_ho_ph_v1_V[j*9 + vst_2x2v_tensor_p1_ho_ph_v1_vrmap[k]]*(vst_2x2v_tensor_p1_ho_ph_v1_coefr[k]*f_r[k]); 
  } 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 4; ++a) { 
    f_l_quad += vst_2x2v_tensor_p1_ho_ph_v1_Cm[i*4 + a]*G_l[a]; 
    f_r_quad += vst_2x2v_tensor_p1_ho_ph_v1_Cm[i*4 + a]*G_r[a]; 
  } 
  const int n = i*4 + j; 
  const double jac = jacob_vel_surf_r[0 + j]; 
  Fhat_nodal[16 + n] = 0.5*jac*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vy_2x2v_tensor_p1_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vel_surf_vy = &jacob_vel_surf_r[4]; 
  return 1.5*dv11*alpha_max/(-(0.11391719628198968*jacob_vel_surf_vy[3])+0.40076152031165013*jacob_vel_surf_vy[2]-0.8136324494869249*jacob_vel_surf_vy[1]+1.5267881254572662*jacob_vel_surf_vy[0]);
} 

GKYL_CU_DH double ho_lax_flux_nodal_vy_2x2v_tensor_p1(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 4; ++j) { 
      alpha_max = fmax(alpha_max, ho_lax_flux_nodal_vy_2x2v_tensor_p1_node(i, j, jacob_vel_surf_r, alpha_quad[i*4 + j], f_l, f_r, Fhat_nodal)); 
    } 
  } 
  return ho_lax_flux_nodal_vy_2x2v_tensor_p1_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
