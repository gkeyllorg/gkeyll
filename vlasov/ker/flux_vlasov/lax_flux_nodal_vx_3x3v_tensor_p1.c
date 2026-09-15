#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_3x3v_tensor_p1.h> 
GKYL_CU_DH double lax_flux_nodal_vx_3x3v_tensor_p1_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *f_l, const double *f_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double G_l[8]; 
  double G_r[8]; 
  for (int a = 0; a < 8; ++a) { G_l[a] = 0.0; G_r[a] = 0.0; } 
  for (int k = 0; k < 216; ++k) { 
    const int a = vst_3x3v_tensor_p1_ph_v0_cmap[k]; 
    G_l[a] += vst_3x3v_tensor_p1_ph_v0_V[j*27 + vst_3x3v_tensor_p1_ph_v0_vlmap[k]]*(vst_3x3v_tensor_p1_ph_v0_coefl[k]*f_l[k]); 
    G_r[a] += vst_3x3v_tensor_p1_ph_v0_V[j*27 + vst_3x3v_tensor_p1_ph_v0_vrmap[k]]*(vst_3x3v_tensor_p1_ph_v0_coefr[k]*f_r[k]); 
  } 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 8; ++a) { 
    f_l_quad += vst_3x3v_tensor_p1_ph_v0_Cm[i*8 + a]*G_l[a]; 
    f_r_quad += vst_3x3v_tensor_p1_ph_v0_Cm[i*8 + a]*G_r[a]; 
  } 
  const int n = i*9 + j; 
  const double jac = jacob_vel_surf_r[3 + j/3]*jacob_vel_surf_r[6 + j%3]; 
  Fhat_nodal[0 + n] = 0.5*jac*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH double lax_flux_nodal_vx_3x3v_tensor_p1_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv10 = 2.0/dxv[3]; 
  const double *jacob_vel_surf_vx = &jacob_vel_surf_r[0]; 
  return 1.5*dv10*alpha_max/(0.18783610896543046*jacob_vel_surf_vx[2]-0.6666666666666666*jacob_vel_surf_vx[1]+1.4788305577012362*jacob_vel_surf_vx[0]);
} 

GKYL_CU_DH double lax_flux_nodal_vx_3x3v_tensor_p1(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double alpha_max = 0.0; 
  for (int i = 0; i < 8; ++i) { 
    for (int j = 0; j < 9; ++j) { 
      alpha_max = fmax(alpha_max, lax_flux_nodal_vx_3x3v_tensor_p1_node(i, j, jacob_vel_surf_r, alpha_quad[i*9 + j], f_l, f_r, Fhat_nodal)); 
    } 
  } 
  return lax_flux_nodal_vx_3x3v_tensor_p1_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
