#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH double lax_flux_y_2x_tensor_p1(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
       double* GKYL_RESTRICT conf_flux_surf) 
{ 
  double dx11 = 2.0/dxv[1]; 
  const double jacob_pos_l_inv = 1.0/jacob_pos_l[2]; 
  const double jacob_pos_r_inv = 1.0/jacob_pos_r[2]; 

  double *out = &conf_flux_surf[16]; 
  double alpha_max = 0.0; 
  const double *em_l; 
  const double *em_r; 
  const double *em_no_J_l; 
  const double *em_no_J_r; 
  const double *flux_l_quad; 
  const double *flux_r_quad; 
  for (int i=0; i<8; ++i) {
    em_l = &field_con_l[i*4]; 
    em_r = &field_con_r[i*4]; 
    em_no_J_l = &field_no_J_con_l[i*4]; 
    em_no_J_r = &field_no_J_con_r[i*4]; 
    flux_l_quad = &flux_l[i*2]; 
    flux_r_quad = &flux_r[i*2]; 
    double Q_l_quad = 0.0; 
    double Q_r_quad = 0.0; 
    double U_l_quad = 0.0; 
    double U_r_quad = 0.0; 
    // At the pole max_alpha_quad includes Jc multiplication already to get the correct asymptotic limit from pre-computed values 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[0])); 
    Q_l_quad = jacob_pos_l_inv*(-(0.8660254037844386*em_l[3])+0.8660254037844386*em_l[2]-0.5*em_l[1]+0.5*em_l[0]); 
    Q_r_quad = jacob_pos_r_inv*(0.8660254037844386*em_r[3]-0.8660254037844386*em_r[2]-0.5*em_r[1]+0.5*em_r[0]); 
    U_l_quad = jacob_pos_l_inv*(-(0.8660254037844386*em_no_J_l[3])+0.8660254037844386*em_no_J_l[2]-0.5*em_no_J_l[1]+0.5*em_no_J_l[0]); 
    U_r_quad = jacob_pos_r_inv*(0.8660254037844386*em_no_J_r[3]-0.8660254037844386*em_no_J_r[2]-0.5*em_no_J_r[1]+0.5*em_no_J_r[0]); 
    if (theta_pole) { 
      out[0+i*2] = 0.5*((flux_r_quad[0] + flux_l_quad[0]) - fabs(max_alpha_quad[0])*(U_r_quad - U_l_quad)); 
    }
    else {
      out[0+i*2] = 0.5*((flux_r_quad[0] + flux_l_quad[0]) - fabs(max_alpha_quad[0])*(Q_r_quad - Q_l_quad)); 
    }
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[1])); 
    Q_l_quad = jacob_pos_l_inv*(0.8660254037844386*em_l[3]+0.8660254037844386*em_l[2]+0.5*em_l[1]+0.5*em_l[0]); 
    Q_r_quad = jacob_pos_r_inv*(-(0.8660254037844386*em_r[3])-0.8660254037844386*em_r[2]+0.5*em_r[1]+0.5*em_r[0]); 
    U_l_quad = jacob_pos_l_inv*(0.8660254037844386*em_no_J_l[3]+0.8660254037844386*em_no_J_l[2]+0.5*em_no_J_l[1]+0.5*em_no_J_l[0]); 
    U_r_quad = jacob_pos_r_inv*(-(0.8660254037844386*em_no_J_r[3])-0.8660254037844386*em_no_J_r[2]+0.5*em_no_J_r[1]+0.5*em_no_J_r[0]); 
    if (theta_pole) { 
      out[1+i*2] = 0.5*((flux_r_quad[1] + flux_l_quad[1]) - fabs(max_alpha_quad[1])*(U_r_quad - U_l_quad)); 
    }
    else {
      out[1+i*2] = 0.5*((flux_r_quad[1] + flux_l_quad[1]) - fabs(max_alpha_quad[1])*(Q_r_quad - Q_l_quad)); 
    }
  }
  
  return 1.5*dx11*fmax(jacob_pos_l_inv,jacob_pos_r_inv)*alpha_max;

} 
