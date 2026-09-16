#include <gkyl_flux_vlasov_kernels.h> 
#include <gkyl_vlasov_flux_surf_mod2nod_tables_1x3v_ser_p2.h> 
GKYL_CU_DH double ho_lax_flux_nodal_vx_1x3v_ser_p2_node(int i, int j, const double *jacob_vel_surf_r,
  double alpha, const double *f_l, const double *f_r, double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double G_l[3]; 
  double G_r[3]; 
  for (int a = 0; a < 3; ++a) { G_l[a] = 0.0; G_r[a] = 0.0; } 
  for (int k = 0; k < 48; ++k) { 
    const int a = vst_1x3v_ser_p2_ho_ph_v0_cmap[k]; 
    G_l[a] += vst_1x3v_ser_p2_ho_ph_v0_V[j*20 + vst_1x3v_ser_p2_ho_ph_v0_vlmap[k]]*(vst_1x3v_ser_p2_ho_ph_v0_coefl[k]*f_l[k]); 
    G_r[a] += vst_1x3v_ser_p2_ho_ph_v0_V[j*20 + vst_1x3v_ser_p2_ho_ph_v0_vrmap[k]]*(vst_1x3v_ser_p2_ho_ph_v0_coefr[k]*f_r[k]); 
  } 
  double f_l_quad = 0.0; 
  double f_r_quad = 0.0; 
  for (int a = 0; a < 3; ++a) { 
    f_l_quad += vst_1x3v_ser_p2_ho_ph_v0_Cm[i*3 + a]*G_l[a]; 
    f_r_quad += vst_1x3v_ser_p2_ho_ph_v0_Cm[i*3 + a]*G_r[a]; 
  } 
  const int n = i*16 + j; 
  const double jac = jacob_vel_surf_r[4]*jacob_vel_surf_r[8]; 
  Fhat_nodal[0 + n] = 0.5*jac*(alpha*(f_r_quad + f_l_quad) - fabs(alpha)*(f_r_quad - f_l_quad)); 
  return fabs(alpha); 
} 

GKYL_CU_DH double ho_lax_flux_nodal_vx_1x3v_ser_p2_cfl(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r, double alpha_max) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double jacob_vel_surf_min = fmin(jacob_vel_surf_l[0], jacob_vel_surf_r[0]); 
  return 2.5*dv10*alpha_max/jacob_vel_surf_min;
} 

GKYL_CU_DH double ho_lax_flux_nodal_vx_1x3v_ser_p2(const double *dxv, const double *jacob_vel_surf_l, const double *jacob_vel_surf_r,
  const double *alpha_quad, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT Fhat_nodal) 
{ 
  double alpha_max = 0.0; 
  for (int i = 0; i < 4; ++i) { 
    for (int j = 0; j < 16; ++j) { 
      alpha_max = fmax(alpha_max, ho_lax_flux_nodal_vx_1x3v_ser_p2_node(i, j, jacob_vel_surf_r, alpha_quad[i*16 + j], f_l, f_r, Fhat_nodal)); 
    } 
  } 
  return ho_lax_flux_nodal_vx_1x3v_ser_p2_cfl(dxv, jacob_vel_surf_l, jacob_vel_surf_r, alpha_max); 
} 
