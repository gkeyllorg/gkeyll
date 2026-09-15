#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH double lax_flux_x_2x_tensor_p2(const double *dxv, const int theta_pole,
      const double *jacob_pos_l, const double *jacob_pos_r,
      const double *J_c, const double *flux_l, const double *flux_r, const double *max_alpha_quad, 
      const double *field_con_l, const double *field_con_r, const double *field_no_J_con_l, const double *field_no_J_con_r,
       double* GKYL_RESTRICT conf_flux_surf) 
{ 
  double dx10 = 2.0/dxv[0]; 
  const double jacob_pos_l_inv = 1.0/jacob_pos_l[0]; 
  const double jacob_pos_r_inv = 1.0/jacob_pos_r[0]; 

  double *out = &conf_flux_surf[0]; 
  double alpha_max = 0.0; 
  const double *em_l; 
  const double *em_r; 
  const double *em_no_J_l; 
  const double *em_no_J_r; 
  const double *flux_l_quad; 
  const double *flux_r_quad; 
  for (int i=0; i<8; ++i) {
    em_l = &field_con_l[i*9]; 
    em_r = &field_con_r[i*9]; 
    em_no_J_l = &field_no_J_con_l[i*9]; 
    em_no_J_r = &field_no_J_con_r[i*9]; 
    flux_l_quad = &flux_l[i*3]; 
    flux_r_quad = &flux_r[i*3]; 
    double Q_l_quad = 0.0; 
    double Q_r_quad = 0.0; 
    double U_l_quad = 0.0; 
    double U_r_quad = 0.0; 
    // At the pole max_alpha_quad includes Jc multiplication already to get the correct asymptotic limit from pre-computed values 
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[0])); 
    Q_l_quad = jacob_pos_l_inv*(em_l[8]+0.7745966692414833*em_l[7]-1.5*em_l[6]+0.4472135954999579*em_l[5]+1.118033988749895*em_l[4]-1.1618950038622249*em_l[3]-0.6708203932499369*em_l[2]+0.8660254037844386*em_l[1]+0.5*em_l[0]); 
    Q_r_quad = jacob_pos_r_inv*(em_r[8]-0.7745966692414833*em_r[7]-1.5*em_r[6]+0.4472135954999579*em_r[5]+1.118033988749895*em_r[4]+1.1618950038622249*em_r[3]-0.6708203932499369*em_r[2]-0.8660254037844386*em_r[1]+0.5*em_r[0]); 
    U_l_quad = jacob_pos_l_inv*(em_no_J_l[8]+0.7745966692414833*em_no_J_l[7]-1.5*em_no_J_l[6]+0.4472135954999579*em_no_J_l[5]+1.118033988749895*em_no_J_l[4]-1.1618950038622249*em_no_J_l[3]-0.6708203932499369*em_no_J_l[2]+0.8660254037844386*em_no_J_l[1]+0.5*em_no_J_l[0]); 
    U_r_quad = jacob_pos_r_inv*(em_no_J_r[8]-0.7745966692414833*em_no_J_r[7]-1.5*em_no_J_r[6]+0.4472135954999579*em_no_J_r[5]+1.118033988749895*em_no_J_r[4]+1.1618950038622249*em_no_J_r[3]-0.6708203932499369*em_no_J_r[2]-0.8660254037844386*em_no_J_r[1]+0.5*em_no_J_r[0]); 
    if (theta_pole) { 
      out[0+i*3] = 0.5*((flux_r_quad[0] + flux_l_quad[0]) - fabs(max_alpha_quad[0])*(U_r_quad - U_l_quad)); 
    }
    else {
      out[0+i*3] = 0.5*((flux_r_quad[0] + flux_l_quad[0]) - fabs(max_alpha_quad[0])*(Q_r_quad - Q_l_quad)); 
    }
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[1])); 
    Q_l_quad = jacob_pos_l_inv*(-(1.25*em_l[8])-0.9682458365518543*em_l[7]-0.5590169943749475*em_l[5]+1.118033988749895*em_l[4]+0.8660254037844386*em_l[1]+0.5*em_l[0]); 
    Q_r_quad = jacob_pos_r_inv*(-(1.25*em_r[8])+0.9682458365518543*em_r[7]-0.5590169943749475*em_r[5]+1.118033988749895*em_r[4]-0.8660254037844386*em_r[1]+0.5*em_r[0]); 
    U_l_quad = jacob_pos_l_inv*(-(1.25*em_no_J_l[8])-0.9682458365518543*em_no_J_l[7]-0.5590169943749475*em_no_J_l[5]+1.118033988749895*em_no_J_l[4]+0.8660254037844386*em_no_J_l[1]+0.5*em_no_J_l[0]); 
    U_r_quad = jacob_pos_r_inv*(-(1.25*em_no_J_r[8])+0.9682458365518543*em_no_J_r[7]-0.5590169943749475*em_no_J_r[5]+1.118033988749895*em_no_J_r[4]-0.8660254037844386*em_no_J_r[1]+0.5*em_no_J_r[0]); 
    if (theta_pole) { 
      out[1+i*3] = 0.5*((flux_r_quad[1] + flux_l_quad[1]) - fabs(max_alpha_quad[1])*(U_r_quad - U_l_quad)); 
    }
    else {
      out[1+i*3] = 0.5*((flux_r_quad[1] + flux_l_quad[1]) - fabs(max_alpha_quad[1])*(Q_r_quad - Q_l_quad)); 
    }
    alpha_max = fmax(alpha_max, fabs(max_alpha_quad[2])); 
    Q_l_quad = jacob_pos_l_inv*(em_l[8]+0.7745966692414833*em_l[7]+1.5*em_l[6]+0.4472135954999579*em_l[5]+1.118033988749895*em_l[4]+1.1618950038622249*em_l[3]+0.6708203932499369*em_l[2]+0.8660254037844386*em_l[1]+0.5*em_l[0]); 
    Q_r_quad = jacob_pos_r_inv*(em_r[8]-0.7745966692414833*em_r[7]+1.5*em_r[6]+0.4472135954999579*em_r[5]+1.118033988749895*em_r[4]-1.1618950038622249*em_r[3]+0.6708203932499369*em_r[2]-0.8660254037844386*em_r[1]+0.5*em_r[0]); 
    U_l_quad = jacob_pos_l_inv*(em_no_J_l[8]+0.7745966692414833*em_no_J_l[7]+1.5*em_no_J_l[6]+0.4472135954999579*em_no_J_l[5]+1.118033988749895*em_no_J_l[4]+1.1618950038622249*em_no_J_l[3]+0.6708203932499369*em_no_J_l[2]+0.8660254037844386*em_no_J_l[1]+0.5*em_no_J_l[0]); 
    U_r_quad = jacob_pos_r_inv*(em_no_J_r[8]-0.7745966692414833*em_no_J_r[7]+1.5*em_no_J_r[6]+0.4472135954999579*em_no_J_r[5]+1.118033988749895*em_no_J_r[4]-1.1618950038622249*em_no_J_r[3]+0.6708203932499369*em_no_J_r[2]-0.8660254037844386*em_no_J_r[1]+0.5*em_no_J_r[0]); 
    if (theta_pole) { 
      out[2+i*3] = 0.5*((flux_r_quad[2] + flux_l_quad[2]) - fabs(max_alpha_quad[2])*(U_r_quad - U_l_quad)); 
    }
    else {
      out[2+i*3] = 0.5*((flux_r_quad[2] + flux_l_quad[2]) - fabs(max_alpha_quad[2])*(Q_r_quad - Q_l_quad)); 
    }
  }
  
  return 2.5*dx10*fmax(jacob_pos_l_inv,jacob_pos_r_inv)*alpha_max;

} 
