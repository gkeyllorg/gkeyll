#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_alpha_quad_x_2x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dxv, 
      const int theta_pole, const double *jacob_pos_l, const double *jacob_pos_r,
      const double *lapse_nodal, const double *shift_nodal, const double *h_ij_nodal, 
      const double *h_ij_inv_nodal, const double *J_c, const double *field_con_l, const double *field_con_r, 
      const double *field_no_J_con_l, const double *field_no_J_con_r, 
      double* GKYL_RESTRICT flux_l, double* GKYL_RESTRICT flux_r, double* GKYL_RESTRICT max_alpha_quad) 
{ 
  const double chi = meq->chi, gamma = meq->gamma; 
  const double jacob_pos_l_inv = 1.0/jacob_pos_l[0]; 
  const double jacob_pos_r_inv = 1.0/jacob_pos_r[0]; 
  const double *h_xx_nodal = &h_ij_nodal[0]; 
  const double *h_xy_nodal = &h_ij_nodal[3]; 
  const double *h_xz_nodal = &h_ij_nodal[6]; 
  const double *h_yy_nodal = &h_ij_nodal[9]; 
  const double *h_yz_nodal = &h_ij_nodal[12]; 
  const double *h_zz_nodal = &h_ij_nodal[15]; 
  const double *h_xx_inv_nodal = &h_ij_inv_nodal[0]; 
  const double *h_xy_inv_nodal = &h_ij_inv_nodal[3]; 
  const double *h_xz_inv_nodal = &h_ij_inv_nodal[6]; 
  const double *h_yy_inv_nodal = &h_ij_inv_nodal[9]; 
  const double *h_yz_inv_nodal = &h_ij_inv_nodal[12]; 
  const double *h_zz_inv_nodal = &h_ij_inv_nodal[15]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[3]; 
  const double *shift_nodal_z = &shift_nodal[6]; 
  
  double JDx_con_l_n[3] = {0.0};
  double JDy_con_l_n[3] = {0.0};
  double JDz_con_l_n[3] = {0.0};
  double JBx_con_l_n[3] = {0.0};
  double JBy_con_l_n[3] = {0.0};
  double JBz_con_l_n[3] = {0.0};
  double Jphi_l_n[3] = {0.0};
  double Jpsi_l_n[3] = {0.0};
  double Dx_con_l_n[3] = {0.0};
  double Dy_con_l_n[3] = {0.0};
  double Dz_con_l_n[3] = {0.0};
  double Bx_con_l_n[3] = {0.0};
  double By_con_l_n[3] = {0.0};
  double Bz_con_l_n[3] = {0.0};
  double phi_l_n[3] = {0.0};
  double psi_l_n[3] = {0.0};
  double Ex_l_n[3] = {0.0};
  double Ey_l_n[3] = {0.0};
  double Ez_l_n[3] = {0.0};
  double Hx_l_n[3] = {0.0};
  double Hy_l_n[3] = {0.0};
  double Hz_l_n[3] = {0.0};
  double FD11_l_n[3] = {0.0};
  double FD12_l_n[3] = {0.0};
  double FD13_l_n[3] = {0.0};
  double FD21_l_n[3] = {0.0};
  double FD22_l_n[3] = {0.0};
  double FD23_l_n[3] = {0.0};
  double FD31_l_n[3] = {0.0};
  double FD32_l_n[3] = {0.0};
  double FD33_l_n[3] = {0.0};
  double FB11_l_n[3] = {0.0};
  double FB12_l_n[3] = {0.0};
  double FB13_l_n[3] = {0.0};
  double FB21_l_n[3] = {0.0};
  double FB22_l_n[3] = {0.0};
  double FB23_l_n[3] = {0.0};
  double FB31_l_n[3] = {0.0};
  double FB32_l_n[3] = {0.0};
  double FB33_l_n[3] = {0.0};
  double FPhi1_l_n[3] = {0.0};
  double FPhi2_l_n[3] = {0.0};
  double FPhi3_l_n[3] = {0.0};
  double FPsi1_l_n[3] = {0.0};
  double FPsi2_l_n[3] = {0.0};
  double FPsi3_l_n[3] = {0.0};
  
  double JDx_con_r_n[3] = {0.0};
  double JDy_con_r_n[3] = {0.0};
  double JDz_con_r_n[3] = {0.0};
  double JBx_con_r_n[3] = {0.0};
  double JBy_con_r_n[3] = {0.0};
  double JBz_con_r_n[3] = {0.0};
  double Jphi_r_n[3] = {0.0};
  double Jpsi_r_n[3] = {0.0};
  double Dx_con_r_n[3] = {0.0};
  double Dy_con_r_n[3] = {0.0};
  double Dz_con_r_n[3] = {0.0};
  double Bx_con_r_n[3] = {0.0};
  double By_con_r_n[3] = {0.0};
  double Bz_con_r_n[3] = {0.0};
  double phi_r_n[3] = {0.0};
  double psi_r_n[3] = {0.0};
  double Ex_r_n[3] = {0.0};
  double Ey_r_n[3] = {0.0};
  double Ez_r_n[3] = {0.0};
  double Hx_r_n[3] = {0.0};
  double Hy_r_n[3] = {0.0};
  double Hz_r_n[3] = {0.0};
  double FD11_r_n[3] = {0.0};
  double FD12_r_n[3] = {0.0};
  double FD13_r_n[3] = {0.0};
  double FD21_r_n[3] = {0.0};
  double FD22_r_n[3] = {0.0};
  double FD23_r_n[3] = {0.0};
  double FD31_r_n[3] = {0.0};
  double FD32_r_n[3] = {0.0};
  double FD33_r_n[3] = {0.0};
  double FB11_r_n[3] = {0.0};
  double FB12_r_n[3] = {0.0};
  double FB13_r_n[3] = {0.0};
  double FB21_r_n[3] = {0.0};
  double FB22_r_n[3] = {0.0};
  double FB23_r_n[3] = {0.0};
  double FB31_r_n[3] = {0.0};
  double FB32_r_n[3] = {0.0};
  double FB33_r_n[3] = {0.0};
  double FPhi1_r_n[3] = {0.0};
  double FPhi2_r_n[3] = {0.0};
  double FPhi3_r_n[3] = {0.0};
  double FPsi1_r_n[3] = {0.0};
  double FPsi2_r_n[3] = {0.0};
  double FPsi3_r_n[3] = {0.0};
  
  const double *JDx_l = &field_con_l[0]; 
  const double *JDy_l = &field_con_l[8]; 
  const double *JDz_l = &field_con_l[16]; 
  const double *JBx_l = &field_con_l[24]; 
  const double *JBy_l = &field_con_l[32]; 
  const double *JBz_l = &field_con_l[40]; 
  const double *Jphi_l = &field_con_l[48]; 
  const double *Jpsi_l = &field_con_l[56]; 
  
  JDx_con_l_n[0] = 0.7745966692414833*JDx_l[7]-1.5*JDx_l[6]+0.4472135954999579*JDx_l[5]+1.118033988749895*JDx_l[4]-1.1618950038622249*JDx_l[3]-0.6708203932499369*JDx_l[2]+0.8660254037844386*JDx_l[1]+0.5*JDx_l[0];
  JDy_con_l_n[0] = 0.7745966692414833*JDy_l[7]-1.5*JDy_l[6]+0.4472135954999579*JDy_l[5]+1.118033988749895*JDy_l[4]-1.1618950038622249*JDy_l[3]-0.6708203932499369*JDy_l[2]+0.8660254037844386*JDy_l[1]+0.5*JDy_l[0];
  JDz_con_l_n[0] = 0.7745966692414833*JDz_l[7]-1.5*JDz_l[6]+0.4472135954999579*JDz_l[5]+1.118033988749895*JDz_l[4]-1.1618950038622249*JDz_l[3]-0.6708203932499369*JDz_l[2]+0.8660254037844386*JDz_l[1]+0.5*JDz_l[0];
  JBx_con_l_n[0] = 0.7745966692414833*JBx_l[7]-1.5*JBx_l[6]+0.4472135954999579*JBx_l[5]+1.118033988749895*JBx_l[4]-1.1618950038622249*JBx_l[3]-0.6708203932499369*JBx_l[2]+0.8660254037844386*JBx_l[1]+0.5*JBx_l[0];
  JBy_con_l_n[0] = 0.7745966692414833*JBy_l[7]-1.5*JBy_l[6]+0.4472135954999579*JBy_l[5]+1.118033988749895*JBy_l[4]-1.1618950038622249*JBy_l[3]-0.6708203932499369*JBy_l[2]+0.8660254037844386*JBy_l[1]+0.5*JBy_l[0];
  JBz_con_l_n[0] = 0.7745966692414833*JBz_l[7]-1.5*JBz_l[6]+0.4472135954999579*JBz_l[5]+1.118033988749895*JBz_l[4]-1.1618950038622249*JBz_l[3]-0.6708203932499369*JBz_l[2]+0.8660254037844386*JBz_l[1]+0.5*JBz_l[0];
  Jphi_l_n[0] = 0.7745966692414833*Jphi_l[7]-1.5*Jphi_l[6]+0.4472135954999579*Jphi_l[5]+1.118033988749895*Jphi_l[4]-1.1618950038622249*Jphi_l[3]-0.6708203932499369*Jphi_l[2]+0.8660254037844386*Jphi_l[1]+0.5*Jphi_l[0];
  Jpsi_l_n[0] = 0.7745966692414833*Jpsi_l[7]-1.5*Jpsi_l[6]+0.4472135954999579*Jpsi_l[5]+1.118033988749895*Jpsi_l[4]-1.1618950038622249*Jpsi_l[3]-0.6708203932499369*Jpsi_l[2]+0.8660254037844386*Jpsi_l[1]+0.5*Jpsi_l[0];
  
  JDx_con_l_n[1] = -(0.9682458365518543*JDx_l[7])-0.5590169943749475*JDx_l[5]+1.118033988749895*JDx_l[4]+0.8660254037844386*JDx_l[1]+0.5*JDx_l[0];
  JDy_con_l_n[1] = -(0.9682458365518543*JDy_l[7])-0.5590169943749475*JDy_l[5]+1.118033988749895*JDy_l[4]+0.8660254037844386*JDy_l[1]+0.5*JDy_l[0];
  JDz_con_l_n[1] = -(0.9682458365518543*JDz_l[7])-0.5590169943749475*JDz_l[5]+1.118033988749895*JDz_l[4]+0.8660254037844386*JDz_l[1]+0.5*JDz_l[0];
  JBx_con_l_n[1] = -(0.9682458365518543*JBx_l[7])-0.5590169943749475*JBx_l[5]+1.118033988749895*JBx_l[4]+0.8660254037844386*JBx_l[1]+0.5*JBx_l[0];
  JBy_con_l_n[1] = -(0.9682458365518543*JBy_l[7])-0.5590169943749475*JBy_l[5]+1.118033988749895*JBy_l[4]+0.8660254037844386*JBy_l[1]+0.5*JBy_l[0];
  JBz_con_l_n[1] = -(0.9682458365518543*JBz_l[7])-0.5590169943749475*JBz_l[5]+1.118033988749895*JBz_l[4]+0.8660254037844386*JBz_l[1]+0.5*JBz_l[0];
  Jphi_l_n[1] = -(0.9682458365518543*Jphi_l[7])-0.5590169943749475*Jphi_l[5]+1.118033988749895*Jphi_l[4]+0.8660254037844386*Jphi_l[1]+0.5*Jphi_l[0];
  Jpsi_l_n[1] = -(0.9682458365518543*Jpsi_l[7])-0.5590169943749475*Jpsi_l[5]+1.118033988749895*Jpsi_l[4]+0.8660254037844386*Jpsi_l[1]+0.5*Jpsi_l[0];
  
  JDx_con_l_n[2] = 0.7745966692414833*JDx_l[7]+1.5*JDx_l[6]+0.4472135954999579*JDx_l[5]+1.118033988749895*JDx_l[4]+1.1618950038622249*JDx_l[3]+0.6708203932499369*JDx_l[2]+0.8660254037844386*JDx_l[1]+0.5*JDx_l[0];
  JDy_con_l_n[2] = 0.7745966692414833*JDy_l[7]+1.5*JDy_l[6]+0.4472135954999579*JDy_l[5]+1.118033988749895*JDy_l[4]+1.1618950038622249*JDy_l[3]+0.6708203932499369*JDy_l[2]+0.8660254037844386*JDy_l[1]+0.5*JDy_l[0];
  JDz_con_l_n[2] = 0.7745966692414833*JDz_l[7]+1.5*JDz_l[6]+0.4472135954999579*JDz_l[5]+1.118033988749895*JDz_l[4]+1.1618950038622249*JDz_l[3]+0.6708203932499369*JDz_l[2]+0.8660254037844386*JDz_l[1]+0.5*JDz_l[0];
  JBx_con_l_n[2] = 0.7745966692414833*JBx_l[7]+1.5*JBx_l[6]+0.4472135954999579*JBx_l[5]+1.118033988749895*JBx_l[4]+1.1618950038622249*JBx_l[3]+0.6708203932499369*JBx_l[2]+0.8660254037844386*JBx_l[1]+0.5*JBx_l[0];
  JBy_con_l_n[2] = 0.7745966692414833*JBy_l[7]+1.5*JBy_l[6]+0.4472135954999579*JBy_l[5]+1.118033988749895*JBy_l[4]+1.1618950038622249*JBy_l[3]+0.6708203932499369*JBy_l[2]+0.8660254037844386*JBy_l[1]+0.5*JBy_l[0];
  JBz_con_l_n[2] = 0.7745966692414833*JBz_l[7]+1.5*JBz_l[6]+0.4472135954999579*JBz_l[5]+1.118033988749895*JBz_l[4]+1.1618950038622249*JBz_l[3]+0.6708203932499369*JBz_l[2]+0.8660254037844386*JBz_l[1]+0.5*JBz_l[0];
  Jphi_l_n[2] = 0.7745966692414833*Jphi_l[7]+1.5*Jphi_l[6]+0.4472135954999579*Jphi_l[5]+1.118033988749895*Jphi_l[4]+1.1618950038622249*Jphi_l[3]+0.6708203932499369*Jphi_l[2]+0.8660254037844386*Jphi_l[1]+0.5*Jphi_l[0];
  Jpsi_l_n[2] = 0.7745966692414833*Jpsi_l[7]+1.5*Jpsi_l[6]+0.4472135954999579*Jpsi_l[5]+1.118033988749895*Jpsi_l[4]+1.1618950038622249*Jpsi_l[3]+0.6708203932499369*Jpsi_l[2]+0.8660254037844386*Jpsi_l[1]+0.5*Jpsi_l[0];
  
  const double *Dx_l = &field_no_J_con_l[0]; 
  const double *Dy_l = &field_no_J_con_l[8]; 
  const double *Dz_l = &field_no_J_con_l[16]; 
  const double *Bx_l = &field_no_J_con_l[24]; 
  const double *By_l = &field_no_J_con_l[32]; 
  const double *Bz_l = &field_no_J_con_l[40]; 
  const double *phi_l = &field_no_J_con_l[48]; 
  const double *psi_l = &field_no_J_con_l[56]; 
  
  Dx_con_l_n[0] = 0.7745966692414833*Dx_l[7]-1.5*Dx_l[6]+0.4472135954999579*Dx_l[5]+1.118033988749895*Dx_l[4]-1.1618950038622249*Dx_l[3]-0.6708203932499369*Dx_l[2]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
  Dy_con_l_n[0] = 0.7745966692414833*Dy_l[7]-1.5*Dy_l[6]+0.4472135954999579*Dy_l[5]+1.118033988749895*Dy_l[4]-1.1618950038622249*Dy_l[3]-0.6708203932499369*Dy_l[2]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
  Dz_con_l_n[0] = 0.7745966692414833*Dz_l[7]-1.5*Dz_l[6]+0.4472135954999579*Dz_l[5]+1.118033988749895*Dz_l[4]-1.1618950038622249*Dz_l[3]-0.6708203932499369*Dz_l[2]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
  Bx_con_l_n[0] = 0.7745966692414833*Bx_l[7]-1.5*Bx_l[6]+0.4472135954999579*Bx_l[5]+1.118033988749895*Bx_l[4]-1.1618950038622249*Bx_l[3]-0.6708203932499369*Bx_l[2]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
  By_con_l_n[0] = 0.7745966692414833*By_l[7]-1.5*By_l[6]+0.4472135954999579*By_l[5]+1.118033988749895*By_l[4]-1.1618950038622249*By_l[3]-0.6708203932499369*By_l[2]+0.8660254037844386*By_l[1]+0.5*By_l[0];
  Bz_con_l_n[0] = 0.7745966692414833*Bz_l[7]-1.5*Bz_l[6]+0.4472135954999579*Bz_l[5]+1.118033988749895*Bz_l[4]-1.1618950038622249*Bz_l[3]-0.6708203932499369*Bz_l[2]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];
  phi_l_n[0] = 0.7745966692414833*phi_l[7]-1.5*phi_l[6]+0.4472135954999579*phi_l[5]+1.118033988749895*phi_l[4]-1.1618950038622249*phi_l[3]-0.6708203932499369*phi_l[2]+0.8660254037844386*phi_l[1]+0.5*phi_l[0];
  psi_l_n[0] = 0.7745966692414833*psi_l[7]-1.5*psi_l[6]+0.4472135954999579*psi_l[5]+1.118033988749895*psi_l[4]-1.1618950038622249*psi_l[3]-0.6708203932499369*psi_l[2]+0.8660254037844386*psi_l[1]+0.5*psi_l[0];
  
  Dx_con_l_n[1] = -(0.9682458365518543*Dx_l[7])-0.5590169943749475*Dx_l[5]+1.118033988749895*Dx_l[4]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
  Dy_con_l_n[1] = -(0.9682458365518543*Dy_l[7])-0.5590169943749475*Dy_l[5]+1.118033988749895*Dy_l[4]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
  Dz_con_l_n[1] = -(0.9682458365518543*Dz_l[7])-0.5590169943749475*Dz_l[5]+1.118033988749895*Dz_l[4]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
  Bx_con_l_n[1] = -(0.9682458365518543*Bx_l[7])-0.5590169943749475*Bx_l[5]+1.118033988749895*Bx_l[4]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
  By_con_l_n[1] = -(0.9682458365518543*By_l[7])-0.5590169943749475*By_l[5]+1.118033988749895*By_l[4]+0.8660254037844386*By_l[1]+0.5*By_l[0];
  Bz_con_l_n[1] = -(0.9682458365518543*Bz_l[7])-0.5590169943749475*Bz_l[5]+1.118033988749895*Bz_l[4]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];
  phi_l_n[1] = -(0.9682458365518543*phi_l[7])-0.5590169943749475*phi_l[5]+1.118033988749895*phi_l[4]+0.8660254037844386*phi_l[1]+0.5*phi_l[0];
  psi_l_n[1] = -(0.9682458365518543*psi_l[7])-0.5590169943749475*psi_l[5]+1.118033988749895*psi_l[4]+0.8660254037844386*psi_l[1]+0.5*psi_l[0];
  
  Dx_con_l_n[2] = 0.7745966692414833*Dx_l[7]+1.5*Dx_l[6]+0.4472135954999579*Dx_l[5]+1.118033988749895*Dx_l[4]+1.1618950038622249*Dx_l[3]+0.6708203932499369*Dx_l[2]+0.8660254037844386*Dx_l[1]+0.5*Dx_l[0];
  Dy_con_l_n[2] = 0.7745966692414833*Dy_l[7]+1.5*Dy_l[6]+0.4472135954999579*Dy_l[5]+1.118033988749895*Dy_l[4]+1.1618950038622249*Dy_l[3]+0.6708203932499369*Dy_l[2]+0.8660254037844386*Dy_l[1]+0.5*Dy_l[0];
  Dz_con_l_n[2] = 0.7745966692414833*Dz_l[7]+1.5*Dz_l[6]+0.4472135954999579*Dz_l[5]+1.118033988749895*Dz_l[4]+1.1618950038622249*Dz_l[3]+0.6708203932499369*Dz_l[2]+0.8660254037844386*Dz_l[1]+0.5*Dz_l[0];
  Bx_con_l_n[2] = 0.7745966692414833*Bx_l[7]+1.5*Bx_l[6]+0.4472135954999579*Bx_l[5]+1.118033988749895*Bx_l[4]+1.1618950038622249*Bx_l[3]+0.6708203932499369*Bx_l[2]+0.8660254037844386*Bx_l[1]+0.5*Bx_l[0];
  By_con_l_n[2] = 0.7745966692414833*By_l[7]+1.5*By_l[6]+0.4472135954999579*By_l[5]+1.118033988749895*By_l[4]+1.1618950038622249*By_l[3]+0.6708203932499369*By_l[2]+0.8660254037844386*By_l[1]+0.5*By_l[0];
  Bz_con_l_n[2] = 0.7745966692414833*Bz_l[7]+1.5*Bz_l[6]+0.4472135954999579*Bz_l[5]+1.118033988749895*Bz_l[4]+1.1618950038622249*Bz_l[3]+0.6708203932499369*Bz_l[2]+0.8660254037844386*Bz_l[1]+0.5*Bz_l[0];
  phi_l_n[2] = 0.7745966692414833*phi_l[7]+1.5*phi_l[6]+0.4472135954999579*phi_l[5]+1.118033988749895*phi_l[4]+1.1618950038622249*phi_l[3]+0.6708203932499369*phi_l[2]+0.8660254037844386*phi_l[1]+0.5*phi_l[0];
  psi_l_n[2] = 0.7745966692414833*psi_l[7]+1.5*psi_l[6]+0.4472135954999579*psi_l[5]+1.118033988749895*psi_l[4]+1.1618950038622249*psi_l[3]+0.6708203932499369*psi_l[2]+0.8660254037844386*psi_l[1]+0.5*psi_l[0];
  
  // If at a theta pole, zero out B^(theta), and D^(theta) 
  if ( theta_pole ) { 
    Dy_con_l_n[0] = 0.0;
    By_con_l_n[0] = 0.0;
    Dy_con_l_n[1] = 0.0;
    By_con_l_n[1] = 0.0;
    Dy_con_l_n[2] = 0.0;
    By_con_l_n[2] = 0.0;
  } 
  
  // Away from the theta pole use conservative variables, otherwise use primative 
  if ( theta_pole == 0 ) { 
    for (int i=0; i<3; ++i) {
      Hx_l_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*JBx_con_l_n[i] + h_xy_nodal[i]*JBy_con_l_n[i] + h_xz_nodal[i]*JBz_con_l_n[i] ) / J_c[i] - ( shift_nodal_y[i]*JDz_con_l_n[i] - shift_nodal_z[i]*JDy_con_l_n[i]); 
      Hy_l_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*JBx_con_l_n[i] + h_yy_nodal[i]*JBy_con_l_n[i] + h_yz_nodal[i]*JBz_con_l_n[i] ) / J_c[i] - ( - shift_nodal_x[i]*JDz_con_l_n[i] + shift_nodal_z[i]*JDx_con_l_n[i]); 
      Hz_l_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*JBx_con_l_n[i] + h_yz_nodal[i]*JBy_con_l_n[i] + h_zz_nodal[i]*JBz_con_l_n[i] ) / J_c[i] - ( shift_nodal_x[i]*JDy_con_l_n[i] - shift_nodal_y[i]*JDx_con_l_n[i]); 
      Ex_l_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*JDx_con_l_n[i] + h_xy_nodal[i]*JDy_con_l_n[i] + h_xz_nodal[i]*JDz_con_l_n[i] ) / J_c[i] + ( shift_nodal_y[i]*JBz_con_l_n[i] - shift_nodal_z[i]*JBy_con_l_n[i]); 
      Ey_l_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*JDx_con_l_n[i] + h_yy_nodal[i]*JDy_con_l_n[i] + h_yz_nodal[i]*JDz_con_l_n[i] ) / J_c[i] + ( - shift_nodal_x[i]*JBz_con_l_n[i] + shift_nodal_z[i]*JBx_con_l_n[i]); 
      Ez_l_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*JDx_con_l_n[i] + h_yz_nodal[i]*JDy_con_l_n[i] + h_zz_nodal[i]*JDz_con_l_n[i] ) / J_c[i] + ( shift_nodal_x[i]*JBy_con_l_n[i] - shift_nodal_y[i]*JBx_con_l_n[i]); 
      FD11_l_n[i] = chi*h_xx_inv_nodal[i]*Jphi_l_n[i]; 
      FD12_l_n[i] = chi*h_xy_inv_nodal[i]*Jphi_l_n[i]; 
      FD13_l_n[i] = chi*h_xz_inv_nodal[i]*Jphi_l_n[i]; 
      FD21_l_n[i] = chi*h_xy_inv_nodal[i]*Jphi_l_n[i]; 
      FD22_l_n[i] = chi*h_yy_inv_nodal[i]*Jphi_l_n[i]; 
      FD23_l_n[i] = chi*h_yz_inv_nodal[i]*Jphi_l_n[i]; 
      FD31_l_n[i] = chi*h_xz_inv_nodal[i]*Jphi_l_n[i]; 
      FD32_l_n[i] = chi*h_yz_inv_nodal[i]*Jphi_l_n[i]; 
      FD33_l_n[i] = chi*h_zz_inv_nodal[i]*Jphi_l_n[i]; 
      FB11_l_n[i] = gamma*h_xx_inv_nodal[i]*Jpsi_l_n[i]; 
      FB12_l_n[i] = gamma*h_xy_inv_nodal[i]*Jpsi_l_n[i]; 
      FB13_l_n[i] = gamma*h_xz_inv_nodal[i]*Jpsi_l_n[i]; 
      FB21_l_n[i] = gamma*h_xy_inv_nodal[i]*Jpsi_l_n[i]; 
      FB22_l_n[i] = gamma*h_yy_inv_nodal[i]*Jpsi_l_n[i]; 
      FB23_l_n[i] = gamma*h_yz_inv_nodal[i]*Jpsi_l_n[i]; 
      FB31_l_n[i] = gamma*h_xz_inv_nodal[i]*Jpsi_l_n[i]; 
      FB32_l_n[i] = gamma*h_yz_inv_nodal[i]*Jpsi_l_n[i]; 
      FB33_l_n[i] = gamma*h_zz_inv_nodal[i]*Jpsi_l_n[i]; 
      FPhi1_l_n[i] = chi*JDx_con_l_n[i]; 
      FPhi2_l_n[i] = chi*JDy_con_l_n[i]; 
      FPhi3_l_n[i] = chi*JDz_con_l_n[i]; 
      FPsi1_l_n[i] = gamma*JBx_con_l_n[i]; 
      FPsi2_l_n[i] = gamma*JBy_con_l_n[i]; 
      FPsi3_l_n[i] = gamma*JBz_con_l_n[i]; 
    }
  }
  else {
    for (int i=0; i<3; ++i) {
      Hx_l_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Bx_con_l_n[i] + h_xy_nodal[i]*By_con_l_n[i] + h_xz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( shift_nodal_y[i]*Dz_con_l_n[i] - shift_nodal_z[i]*Dy_con_l_n[i]); 
      Hy_l_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Bx_con_l_n[i] + h_yy_nodal[i]*By_con_l_n[i] + h_yz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( - shift_nodal_x[i]*Dz_con_l_n[i] + shift_nodal_z[i]*Dx_con_l_n[i]); 
      Hz_l_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Bx_con_l_n[i] + h_yz_nodal[i]*By_con_l_n[i] + h_zz_nodal[i]*Bz_con_l_n[i] ) - J_c[i]*( shift_nodal_x[i]*Dy_con_l_n[i] - shift_nodal_y[i]*Dx_con_l_n[i]); 
      Ex_l_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Dx_con_l_n[i] + h_xy_nodal[i]*Dy_con_l_n[i] + h_xz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( shift_nodal_y[i]*Bz_con_l_n[i] - shift_nodal_z[i]*By_con_l_n[i]); 
      Ey_l_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Dx_con_l_n[i] + h_yy_nodal[i]*Dy_con_l_n[i] + h_yz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( - shift_nodal_x[i]*Bz_con_l_n[i] + shift_nodal_z[i]*Bx_con_l_n[i]); 
      Ez_l_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Dx_con_l_n[i] + h_yz_nodal[i]*Dy_con_l_n[i] + h_zz_nodal[i]*Dz_con_l_n[i] ) + J_c[i]*( shift_nodal_x[i]*By_con_l_n[i] - shift_nodal_y[i]*Bx_con_l_n[i]); 
    }
  }
  
  const double *JDx_r = &field_con_r[0]; 
  const double *JDy_r = &field_con_r[8]; 
  const double *JDz_r = &field_con_r[16]; 
  const double *JBx_r = &field_con_r[24]; 
  const double *JBy_r = &field_con_r[32]; 
  const double *JBz_r = &field_con_r[40]; 
  const double *Jphi_r = &field_con_r[48]; 
  const double *Jpsi_r = &field_con_r[56]; 
  
  JDx_con_r_n[0] = -(0.7745966692414833*JDx_r[7])-1.5*JDx_r[6]+0.4472135954999579*JDx_r[5]+1.118033988749895*JDx_r[4]+1.1618950038622249*JDx_r[3]-0.6708203932499369*JDx_r[2]-0.8660254037844386*JDx_r[1]+0.5*JDx_r[0];
  JDy_con_r_n[0] = -(0.7745966692414833*JDy_r[7])-1.5*JDy_r[6]+0.4472135954999579*JDy_r[5]+1.118033988749895*JDy_r[4]+1.1618950038622249*JDy_r[3]-0.6708203932499369*JDy_r[2]-0.8660254037844386*JDy_r[1]+0.5*JDy_r[0];
  JDz_con_r_n[0] = -(0.7745966692414833*JDz_r[7])-1.5*JDz_r[6]+0.4472135954999579*JDz_r[5]+1.118033988749895*JDz_r[4]+1.1618950038622249*JDz_r[3]-0.6708203932499369*JDz_r[2]-0.8660254037844386*JDz_r[1]+0.5*JDz_r[0];
  JBx_con_r_n[0] = -(0.7745966692414833*JBx_r[7])-1.5*JBx_r[6]+0.4472135954999579*JBx_r[5]+1.118033988749895*JBx_r[4]+1.1618950038622249*JBx_r[3]-0.6708203932499369*JBx_r[2]-0.8660254037844386*JBx_r[1]+0.5*JBx_r[0];
  JBy_con_r_n[0] = -(0.7745966692414833*JBy_r[7])-1.5*JBy_r[6]+0.4472135954999579*JBy_r[5]+1.118033988749895*JBy_r[4]+1.1618950038622249*JBy_r[3]-0.6708203932499369*JBy_r[2]-0.8660254037844386*JBy_r[1]+0.5*JBy_r[0];
  JBz_con_r_n[0] = -(0.7745966692414833*JBz_r[7])-1.5*JBz_r[6]+0.4472135954999579*JBz_r[5]+1.118033988749895*JBz_r[4]+1.1618950038622249*JBz_r[3]-0.6708203932499369*JBz_r[2]-0.8660254037844386*JBz_r[1]+0.5*JBz_r[0];
  Jphi_r_n[0] = -(0.7745966692414833*Jphi_r[7])-1.5*Jphi_r[6]+0.4472135954999579*Jphi_r[5]+1.118033988749895*Jphi_r[4]+1.1618950038622249*Jphi_r[3]-0.6708203932499369*Jphi_r[2]-0.8660254037844386*Jphi_r[1]+0.5*Jphi_r[0];
  Jpsi_r_n[0] = -(0.7745966692414833*Jpsi_r[7])-1.5*Jpsi_r[6]+0.4472135954999579*Jpsi_r[5]+1.118033988749895*Jpsi_r[4]+1.1618950038622249*Jpsi_r[3]-0.6708203932499369*Jpsi_r[2]-0.8660254037844386*Jpsi_r[1]+0.5*Jpsi_r[0];
  
  JDx_con_r_n[1] = 0.9682458365518543*JDx_r[7]-0.5590169943749475*JDx_r[5]+1.118033988749895*JDx_r[4]-0.8660254037844386*JDx_r[1]+0.5*JDx_r[0];
  JDy_con_r_n[1] = 0.9682458365518543*JDy_r[7]-0.5590169943749475*JDy_r[5]+1.118033988749895*JDy_r[4]-0.8660254037844386*JDy_r[1]+0.5*JDy_r[0];
  JDz_con_r_n[1] = 0.9682458365518543*JDz_r[7]-0.5590169943749475*JDz_r[5]+1.118033988749895*JDz_r[4]-0.8660254037844386*JDz_r[1]+0.5*JDz_r[0];
  JBx_con_r_n[1] = 0.9682458365518543*JBx_r[7]-0.5590169943749475*JBx_r[5]+1.118033988749895*JBx_r[4]-0.8660254037844386*JBx_r[1]+0.5*JBx_r[0];
  JBy_con_r_n[1] = 0.9682458365518543*JBy_r[7]-0.5590169943749475*JBy_r[5]+1.118033988749895*JBy_r[4]-0.8660254037844386*JBy_r[1]+0.5*JBy_r[0];
  JBz_con_r_n[1] = 0.9682458365518543*JBz_r[7]-0.5590169943749475*JBz_r[5]+1.118033988749895*JBz_r[4]-0.8660254037844386*JBz_r[1]+0.5*JBz_r[0];
  Jphi_r_n[1] = 0.9682458365518543*Jphi_r[7]-0.5590169943749475*Jphi_r[5]+1.118033988749895*Jphi_r[4]-0.8660254037844386*Jphi_r[1]+0.5*Jphi_r[0];
  Jpsi_r_n[1] = 0.9682458365518543*Jpsi_r[7]-0.5590169943749475*Jpsi_r[5]+1.118033988749895*Jpsi_r[4]-0.8660254037844386*Jpsi_r[1]+0.5*Jpsi_r[0];
  
  JDx_con_r_n[2] = -(0.7745966692414833*JDx_r[7])+1.5*JDx_r[6]+0.4472135954999579*JDx_r[5]+1.118033988749895*JDx_r[4]-1.1618950038622249*JDx_r[3]+0.6708203932499369*JDx_r[2]-0.8660254037844386*JDx_r[1]+0.5*JDx_r[0];
  JDy_con_r_n[2] = -(0.7745966692414833*JDy_r[7])+1.5*JDy_r[6]+0.4472135954999579*JDy_r[5]+1.118033988749895*JDy_r[4]-1.1618950038622249*JDy_r[3]+0.6708203932499369*JDy_r[2]-0.8660254037844386*JDy_r[1]+0.5*JDy_r[0];
  JDz_con_r_n[2] = -(0.7745966692414833*JDz_r[7])+1.5*JDz_r[6]+0.4472135954999579*JDz_r[5]+1.118033988749895*JDz_r[4]-1.1618950038622249*JDz_r[3]+0.6708203932499369*JDz_r[2]-0.8660254037844386*JDz_r[1]+0.5*JDz_r[0];
  JBx_con_r_n[2] = -(0.7745966692414833*JBx_r[7])+1.5*JBx_r[6]+0.4472135954999579*JBx_r[5]+1.118033988749895*JBx_r[4]-1.1618950038622249*JBx_r[3]+0.6708203932499369*JBx_r[2]-0.8660254037844386*JBx_r[1]+0.5*JBx_r[0];
  JBy_con_r_n[2] = -(0.7745966692414833*JBy_r[7])+1.5*JBy_r[6]+0.4472135954999579*JBy_r[5]+1.118033988749895*JBy_r[4]-1.1618950038622249*JBy_r[3]+0.6708203932499369*JBy_r[2]-0.8660254037844386*JBy_r[1]+0.5*JBy_r[0];
  JBz_con_r_n[2] = -(0.7745966692414833*JBz_r[7])+1.5*JBz_r[6]+0.4472135954999579*JBz_r[5]+1.118033988749895*JBz_r[4]-1.1618950038622249*JBz_r[3]+0.6708203932499369*JBz_r[2]-0.8660254037844386*JBz_r[1]+0.5*JBz_r[0];
  Jphi_r_n[2] = -(0.7745966692414833*Jphi_r[7])+1.5*Jphi_r[6]+0.4472135954999579*Jphi_r[5]+1.118033988749895*Jphi_r[4]-1.1618950038622249*Jphi_r[3]+0.6708203932499369*Jphi_r[2]-0.8660254037844386*Jphi_r[1]+0.5*Jphi_r[0];
  Jpsi_r_n[2] = -(0.7745966692414833*Jpsi_r[7])+1.5*Jpsi_r[6]+0.4472135954999579*Jpsi_r[5]+1.118033988749895*Jpsi_r[4]-1.1618950038622249*Jpsi_r[3]+0.6708203932499369*Jpsi_r[2]-0.8660254037844386*Jpsi_r[1]+0.5*Jpsi_r[0];
  
  const double *Dx_r = &field_no_J_con_r[0]; 
  const double *Dy_r = &field_no_J_con_r[8]; 
  const double *Dz_r = &field_no_J_con_r[16]; 
  const double *Bx_r = &field_no_J_con_r[24]; 
  const double *By_r = &field_no_J_con_r[32]; 
  const double *Bz_r = &field_no_J_con_r[40]; 
  const double *phi_r = &field_no_J_con_r[48]; 
  const double *psi_r = &field_no_J_con_r[56]; 
  
  Dx_con_r_n[0] = -(0.7745966692414833*Dx_r[7])-1.5*Dx_r[6]+0.4472135954999579*Dx_r[5]+1.118033988749895*Dx_r[4]+1.1618950038622249*Dx_r[3]-0.6708203932499369*Dx_r[2]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
  Dy_con_r_n[0] = -(0.7745966692414833*Dy_r[7])-1.5*Dy_r[6]+0.4472135954999579*Dy_r[5]+1.118033988749895*Dy_r[4]+1.1618950038622249*Dy_r[3]-0.6708203932499369*Dy_r[2]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
  Dz_con_r_n[0] = -(0.7745966692414833*Dz_r[7])-1.5*Dz_r[6]+0.4472135954999579*Dz_r[5]+1.118033988749895*Dz_r[4]+1.1618950038622249*Dz_r[3]-0.6708203932499369*Dz_r[2]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
  Bx_con_r_n[0] = -(0.7745966692414833*Bx_r[7])-1.5*Bx_r[6]+0.4472135954999579*Bx_r[5]+1.118033988749895*Bx_r[4]+1.1618950038622249*Bx_r[3]-0.6708203932499369*Bx_r[2]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
  By_con_r_n[0] = -(0.7745966692414833*By_r[7])-1.5*By_r[6]+0.4472135954999579*By_r[5]+1.118033988749895*By_r[4]+1.1618950038622249*By_r[3]-0.6708203932499369*By_r[2]-0.8660254037844386*By_r[1]+0.5*By_r[0];
  Bz_con_r_n[0] = -(0.7745966692414833*Bz_r[7])-1.5*Bz_r[6]+0.4472135954999579*Bz_r[5]+1.118033988749895*Bz_r[4]+1.1618950038622249*Bz_r[3]-0.6708203932499369*Bz_r[2]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];
  phi_r_n[0] = -(0.7745966692414833*phi_r[7])-1.5*phi_r[6]+0.4472135954999579*phi_r[5]+1.118033988749895*phi_r[4]+1.1618950038622249*phi_r[3]-0.6708203932499369*phi_r[2]-0.8660254037844386*phi_r[1]+0.5*phi_r[0];
  psi_r_n[0] = -(0.7745966692414833*psi_r[7])-1.5*psi_r[6]+0.4472135954999579*psi_r[5]+1.118033988749895*psi_r[4]+1.1618950038622249*psi_r[3]-0.6708203932499369*psi_r[2]-0.8660254037844386*psi_r[1]+0.5*psi_r[0];
  
  Dx_con_r_n[1] = 0.9682458365518543*Dx_r[7]-0.5590169943749475*Dx_r[5]+1.118033988749895*Dx_r[4]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
  Dy_con_r_n[1] = 0.9682458365518543*Dy_r[7]-0.5590169943749475*Dy_r[5]+1.118033988749895*Dy_r[4]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
  Dz_con_r_n[1] = 0.9682458365518543*Dz_r[7]-0.5590169943749475*Dz_r[5]+1.118033988749895*Dz_r[4]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
  Bx_con_r_n[1] = 0.9682458365518543*Bx_r[7]-0.5590169943749475*Bx_r[5]+1.118033988749895*Bx_r[4]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
  By_con_r_n[1] = 0.9682458365518543*By_r[7]-0.5590169943749475*By_r[5]+1.118033988749895*By_r[4]-0.8660254037844386*By_r[1]+0.5*By_r[0];
  Bz_con_r_n[1] = 0.9682458365518543*Bz_r[7]-0.5590169943749475*Bz_r[5]+1.118033988749895*Bz_r[4]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];
  phi_r_n[1] = 0.9682458365518543*phi_r[7]-0.5590169943749475*phi_r[5]+1.118033988749895*phi_r[4]-0.8660254037844386*phi_r[1]+0.5*phi_r[0];
  psi_r_n[1] = 0.9682458365518543*psi_r[7]-0.5590169943749475*psi_r[5]+1.118033988749895*psi_r[4]-0.8660254037844386*psi_r[1]+0.5*psi_r[0];
  
  Dx_con_r_n[2] = -(0.7745966692414833*Dx_r[7])+1.5*Dx_r[6]+0.4472135954999579*Dx_r[5]+1.118033988749895*Dx_r[4]-1.1618950038622249*Dx_r[3]+0.6708203932499369*Dx_r[2]-0.8660254037844386*Dx_r[1]+0.5*Dx_r[0];
  Dy_con_r_n[2] = -(0.7745966692414833*Dy_r[7])+1.5*Dy_r[6]+0.4472135954999579*Dy_r[5]+1.118033988749895*Dy_r[4]-1.1618950038622249*Dy_r[3]+0.6708203932499369*Dy_r[2]-0.8660254037844386*Dy_r[1]+0.5*Dy_r[0];
  Dz_con_r_n[2] = -(0.7745966692414833*Dz_r[7])+1.5*Dz_r[6]+0.4472135954999579*Dz_r[5]+1.118033988749895*Dz_r[4]-1.1618950038622249*Dz_r[3]+0.6708203932499369*Dz_r[2]-0.8660254037844386*Dz_r[1]+0.5*Dz_r[0];
  Bx_con_r_n[2] = -(0.7745966692414833*Bx_r[7])+1.5*Bx_r[6]+0.4472135954999579*Bx_r[5]+1.118033988749895*Bx_r[4]-1.1618950038622249*Bx_r[3]+0.6708203932499369*Bx_r[2]-0.8660254037844386*Bx_r[1]+0.5*Bx_r[0];
  By_con_r_n[2] = -(0.7745966692414833*By_r[7])+1.5*By_r[6]+0.4472135954999579*By_r[5]+1.118033988749895*By_r[4]-1.1618950038622249*By_r[3]+0.6708203932499369*By_r[2]-0.8660254037844386*By_r[1]+0.5*By_r[0];
  Bz_con_r_n[2] = -(0.7745966692414833*Bz_r[7])+1.5*Bz_r[6]+0.4472135954999579*Bz_r[5]+1.118033988749895*Bz_r[4]-1.1618950038622249*Bz_r[3]+0.6708203932499369*Bz_r[2]-0.8660254037844386*Bz_r[1]+0.5*Bz_r[0];
  phi_r_n[2] = -(0.7745966692414833*phi_r[7])+1.5*phi_r[6]+0.4472135954999579*phi_r[5]+1.118033988749895*phi_r[4]-1.1618950038622249*phi_r[3]+0.6708203932499369*phi_r[2]-0.8660254037844386*phi_r[1]+0.5*phi_r[0];
  psi_r_n[2] = -(0.7745966692414833*psi_r[7])+1.5*psi_r[6]+0.4472135954999579*psi_r[5]+1.118033988749895*psi_r[4]-1.1618950038622249*psi_r[3]+0.6708203932499369*psi_r[2]-0.8660254037844386*psi_r[1]+0.5*psi_r[0];
  
  // If at a theta pole, zero out B^(theta), and D^(theta) 
  if ( theta_pole ) { 
    Dy_con_r_n[0] = 0.0;
    By_con_r_n[0] = 0.0;
    Dy_con_r_n[1] = 0.0;
    By_con_r_n[1] = 0.0;
    Dy_con_r_n[2] = 0.0;
    By_con_r_n[2] = 0.0;
  } 
  
  // Away from the theta pole use conservative variables, otherwise use primative 
  if ( theta_pole == 0 ) { 
    for (int i=0; i<3; ++i) {
      Hx_r_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*JBx_con_r_n[i] + h_xy_nodal[i]*JBy_con_r_n[i] + h_xz_nodal[i]*JBz_con_r_n[i] ) / J_c[i] - ( shift_nodal_y[i]*JDz_con_r_n[i] - shift_nodal_z[i]*JDy_con_r_n[i]); 
      Hy_r_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*JBx_con_r_n[i] + h_yy_nodal[i]*JBy_con_r_n[i] + h_yz_nodal[i]*JBz_con_r_n[i] ) / J_c[i] - ( - shift_nodal_x[i]*JDz_con_r_n[i] + shift_nodal_z[i]*JDx_con_r_n[i]); 
      Hz_r_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*JBx_con_r_n[i] + h_yz_nodal[i]*JBy_con_r_n[i] + h_zz_nodal[i]*JBz_con_r_n[i] ) / J_c[i] - ( shift_nodal_x[i]*JDy_con_r_n[i] - shift_nodal_y[i]*JDx_con_r_n[i]); 
      Ex_r_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*JDx_con_r_n[i] + h_xy_nodal[i]*JDy_con_r_n[i] + h_xz_nodal[i]*JDz_con_r_n[i] ) / J_c[i] + ( shift_nodal_y[i]*JBz_con_r_n[i] - shift_nodal_z[i]*JBy_con_r_n[i]); 
      Ey_r_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*JDx_con_r_n[i] + h_yy_nodal[i]*JDy_con_r_n[i] + h_yz_nodal[i]*JDz_con_r_n[i] ) / J_c[i] + ( - shift_nodal_x[i]*JBz_con_r_n[i] + shift_nodal_z[i]*JBx_con_r_n[i]); 
      Ez_r_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*JDx_con_r_n[i] + h_yz_nodal[i]*JDy_con_r_n[i] + h_zz_nodal[i]*JDz_con_r_n[i] ) / J_c[i] + ( shift_nodal_x[i]*JBy_con_r_n[i] - shift_nodal_y[i]*JBx_con_r_n[i]); 
      FD11_r_n[i] = chi*h_xx_inv_nodal[i]*Jphi_r_n[i]; 
      FD12_r_n[i] = chi*h_xy_inv_nodal[i]*Jphi_r_n[i]; 
      FD13_r_n[i] = chi*h_xz_inv_nodal[i]*Jphi_r_n[i]; 
      FD21_r_n[i] = chi*h_xy_inv_nodal[i]*Jphi_r_n[i]; 
      FD22_r_n[i] = chi*h_yy_inv_nodal[i]*Jphi_r_n[i]; 
      FD23_r_n[i] = chi*h_yz_inv_nodal[i]*Jphi_r_n[i]; 
      FD31_r_n[i] = chi*h_xz_inv_nodal[i]*Jphi_r_n[i]; 
      FD32_r_n[i] = chi*h_yz_inv_nodal[i]*Jphi_r_n[i]; 
      FD33_r_n[i] = chi*h_zz_inv_nodal[i]*Jphi_r_n[i]; 
      FB11_r_n[i] = gamma*h_xx_inv_nodal[i]*Jpsi_r_n[i]; 
      FB12_r_n[i] = gamma*h_xy_inv_nodal[i]*Jpsi_r_n[i]; 
      FB13_r_n[i] = gamma*h_xz_inv_nodal[i]*Jpsi_r_n[i]; 
      FB21_r_n[i] = gamma*h_xy_inv_nodal[i]*Jpsi_r_n[i]; 
      FB22_r_n[i] = gamma*h_yy_inv_nodal[i]*Jpsi_r_n[i]; 
      FB23_r_n[i] = gamma*h_yz_inv_nodal[i]*Jpsi_r_n[i]; 
      FB31_r_n[i] = gamma*h_xz_inv_nodal[i]*Jpsi_r_n[i]; 
      FB32_r_n[i] = gamma*h_yz_inv_nodal[i]*Jpsi_r_n[i]; 
      FB33_r_n[i] = gamma*h_zz_inv_nodal[i]*Jpsi_r_n[i]; 
      FPhi1_r_n[i] = chi*JDx_con_r_n[i]; 
      FPhi2_r_n[i] = chi*JDy_con_r_n[i]; 
      FPhi3_r_n[i] = chi*JDz_con_r_n[i]; 
      FPsi1_r_n[i] = gamma*JBx_con_r_n[i]; 
      FPsi2_r_n[i] = gamma*JBy_con_r_n[i]; 
      FPsi3_r_n[i] = gamma*JBz_con_r_n[i]; 
    }
  }
  else {
    for (int i=0; i<3; ++i) {
      Hx_r_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Bx_con_r_n[i] + h_xy_nodal[i]*By_con_r_n[i] + h_xz_nodal[i]*Bz_con_r_n[i] ) - J_c[i]*( shift_nodal_y[i]*Dz_con_r_n[i] - shift_nodal_z[i]*Dy_con_r_n[i]); 
      Hy_r_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Bx_con_r_n[i] + h_yy_nodal[i]*By_con_r_n[i] + h_yz_nodal[i]*Bz_con_r_n[i] ) - J_c[i]*( - shift_nodal_x[i]*Dz_con_r_n[i] + shift_nodal_z[i]*Dx_con_r_n[i]); 
      Hz_r_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Bx_con_r_n[i] + h_yz_nodal[i]*By_con_r_n[i] + h_zz_nodal[i]*Bz_con_r_n[i] ) - J_c[i]*( shift_nodal_x[i]*Dy_con_r_n[i] - shift_nodal_y[i]*Dx_con_r_n[i]); 
      Ex_r_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*Dx_con_r_n[i] + h_xy_nodal[i]*Dy_con_r_n[i] + h_xz_nodal[i]*Dz_con_r_n[i] ) + J_c[i]*( shift_nodal_y[i]*Bz_con_r_n[i] - shift_nodal_z[i]*By_con_r_n[i]); 
      Ey_r_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*Dx_con_r_n[i] + h_yy_nodal[i]*Dy_con_r_n[i] + h_yz_nodal[i]*Dz_con_r_n[i] ) + J_c[i]*( - shift_nodal_x[i]*Bz_con_r_n[i] + shift_nodal_z[i]*Bx_con_r_n[i]); 
      Ez_r_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*Dx_con_r_n[i] + h_yz_nodal[i]*Dy_con_r_n[i] + h_zz_nodal[i]*Dz_con_r_n[i] ) + J_c[i]*( shift_nodal_x[i]*By_con_r_n[i] - shift_nodal_y[i]*Bx_con_r_n[i]); 
    }
  }
  
  double *flux_l_quad; 
  double *flux_r_quad; 
  flux_l_quad = &flux_l[0]; 
  flux_r_quad = &flux_r[0]; 
  flux_l_quad[0] = jacob_pos_l_inv*(FD11_l_n[0]); 
  flux_r_quad[0] = jacob_pos_r_inv*(FD11_r_n[0]); 
  flux_l_quad[1] = jacob_pos_l_inv*(FD11_l_n[1]); 
  flux_r_quad[1] = jacob_pos_r_inv*(FD11_r_n[1]); 
  flux_l_quad[2] = jacob_pos_l_inv*(FD11_l_n[2]); 
  flux_r_quad[2] = jacob_pos_r_inv*(FD11_r_n[2]); 
  flux_l_quad = &flux_l[3]; 
  flux_r_quad = &flux_r[3]; 
  flux_l_quad[0] = jacob_pos_l_inv*(Hz_l_n[0]+FD12_l_n[0]); 
  flux_r_quad[0] = jacob_pos_r_inv*(Hz_r_n[0]+FD12_r_n[0]); 
  flux_l_quad[1] = jacob_pos_l_inv*(Hz_l_n[1]+FD12_l_n[1]); 
  flux_r_quad[1] = jacob_pos_r_inv*(Hz_r_n[1]+FD12_r_n[1]); 
  flux_l_quad[2] = jacob_pos_l_inv*(Hz_l_n[2]+FD12_l_n[2]); 
  flux_r_quad[2] = jacob_pos_r_inv*(Hz_r_n[2]+FD12_r_n[2]); 
  flux_l_quad = &flux_l[6]; 
  flux_r_quad = &flux_r[6]; 
  flux_l_quad[0] = jacob_pos_l_inv*(FD13_l_n[0]-1.0*Hy_l_n[0]); 
  flux_r_quad[0] = jacob_pos_r_inv*(FD13_r_n[0]-1.0*Hy_r_n[0]); 
  flux_l_quad[1] = jacob_pos_l_inv*(FD13_l_n[1]-1.0*Hy_l_n[1]); 
  flux_r_quad[1] = jacob_pos_r_inv*(FD13_r_n[1]-1.0*Hy_r_n[1]); 
  flux_l_quad[2] = jacob_pos_l_inv*(FD13_l_n[2]-1.0*Hy_l_n[2]); 
  flux_r_quad[2] = jacob_pos_r_inv*(FD13_r_n[2]-1.0*Hy_r_n[2]); 
  flux_l_quad = &flux_l[9]; 
  flux_r_quad = &flux_r[9]; 
  flux_l_quad[0] = jacob_pos_l_inv*(FB11_l_n[0]); 
  flux_r_quad[0] = jacob_pos_r_inv*(FB11_r_n[0]); 
  flux_l_quad[1] = jacob_pos_l_inv*(FB11_l_n[1]); 
  flux_r_quad[1] = jacob_pos_r_inv*(FB11_r_n[1]); 
  flux_l_quad[2] = jacob_pos_l_inv*(FB11_l_n[2]); 
  flux_r_quad[2] = jacob_pos_r_inv*(FB11_r_n[2]); 
  flux_l_quad = &flux_l[12]; 
  flux_r_quad = &flux_r[12]; 
  flux_l_quad[0] = jacob_pos_l_inv*(FB12_l_n[0]-1.0*Ez_l_n[0]); 
  flux_r_quad[0] = jacob_pos_r_inv*(FB12_r_n[0]-1.0*Ez_r_n[0]); 
  flux_l_quad[1] = jacob_pos_l_inv*(FB12_l_n[1]-1.0*Ez_l_n[1]); 
  flux_r_quad[1] = jacob_pos_r_inv*(FB12_r_n[1]-1.0*Ez_r_n[1]); 
  flux_l_quad[2] = jacob_pos_l_inv*(FB12_l_n[2]-1.0*Ez_l_n[2]); 
  flux_r_quad[2] = jacob_pos_r_inv*(FB12_r_n[2]-1.0*Ez_r_n[2]); 
  flux_l_quad = &flux_l[15]; 
  flux_r_quad = &flux_r[15]; 
  flux_l_quad[0] = jacob_pos_l_inv*(FB13_l_n[0]+Ey_l_n[0]); 
  flux_r_quad[0] = jacob_pos_r_inv*(FB13_r_n[0]+Ey_r_n[0]); 
  flux_l_quad[1] = jacob_pos_l_inv*(FB13_l_n[1]+Ey_l_n[1]); 
  flux_r_quad[1] = jacob_pos_r_inv*(FB13_r_n[1]+Ey_r_n[1]); 
  flux_l_quad[2] = jacob_pos_l_inv*(FB13_l_n[2]+Ey_l_n[2]); 
  flux_r_quad[2] = jacob_pos_r_inv*(FB13_r_n[2]+Ey_r_n[2]); 
  flux_l_quad = &flux_l[18]; 
  flux_r_quad = &flux_r[18]; 
  flux_l_quad[0] = jacob_pos_l_inv*(FPhi1_l_n[0]); 
  flux_r_quad[0] = jacob_pos_r_inv*(FPhi1_r_n[0]); 
  flux_l_quad[1] = jacob_pos_l_inv*(FPhi1_l_n[1]); 
  flux_r_quad[1] = jacob_pos_r_inv*(FPhi1_r_n[1]); 
  flux_l_quad[2] = jacob_pos_l_inv*(FPhi1_l_n[2]); 
  flux_r_quad[2] = jacob_pos_r_inv*(FPhi1_r_n[2]); 
  flux_l_quad = &flux_l[21]; 
  flux_r_quad = &flux_r[21]; 
  flux_l_quad[0] = jacob_pos_l_inv*(FPsi1_l_n[0]); 
  flux_r_quad[0] = jacob_pos_r_inv*(FPsi1_r_n[0]); 
  flux_l_quad[1] = jacob_pos_l_inv*(FPsi1_l_n[1]); 
  flux_r_quad[1] = jacob_pos_r_inv*(FPsi1_r_n[1]); 
  flux_l_quad[2] = jacob_pos_l_inv*(FPsi1_l_n[2]); 
  flux_r_quad[2] = jacob_pos_r_inv*(FPsi1_r_n[2]); 
  
  double lambda_1[3] = {0.0};
  double lambda_2[3] = {0.0};
  double lambda_3[3] = {0.0};
  double lambda_4[3] = {0.0};
  double lambda_5[3] = {0.0};
  double lambda_6[3] = {0.0};
  for (int i=0; i<3; ++i) {
    if (theta_pole == 0) {
    lambda_1[i] =   chi * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i];
    lambda_2[i] = - chi * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i];
    lambda_3[i] =   gamma * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i];
    lambda_4[i] = - gamma * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i];
    lambda_5[i] = -shift_nodal_x[i] + lapse_nodal[i] * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i];
    lambda_6[i] = -shift_nodal_x[i] - lapse_nodal[i] * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] )/J_c[i];
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_1[i] ));
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_2[i] ));
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_3[i] ));
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_4[i] ));
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_5[i] ));
    max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_6[i] ));
    }
    else {      lambda_1[i] =   chi * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] );
      lambda_2[i] = - chi * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] );
      lambda_3[i] =   gamma * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] );
      lambda_4[i] = - gamma * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] );
      lambda_5[i] =   lapse_nodal[i] * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] );
      lambda_6[i] = - lapse_nodal[i] * sqrt( h_yy_nodal[i] * h_zz_nodal[i] - h_yz_nodal[i] * h_yz_nodal[i] );
      max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_1[i] ));
      max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_2[i] ));
      max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_3[i] ));
      max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_4[i] ));
      max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_5[i] ));
      max_alpha_quad[i] = fmax(max_alpha_quad[i], fabs( lambda_6[i] ));
  }
  }

} 
