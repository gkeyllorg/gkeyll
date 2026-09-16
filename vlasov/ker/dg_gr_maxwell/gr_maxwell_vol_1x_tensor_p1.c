#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void gr_maxwell_vol_1x_tensor_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out) 
{ 
  const double chi = meq->chi, gamma = meq->gamma; 
  double dx0 = 2.0/dx[0]; 
  
  double JDx_con_n[2] = {0.0};
  double JDy_con_n[2] = {0.0};
  double JDz_con_n[2] = {0.0};
  double JBx_con_n[2] = {0.0};
  double JBy_con_n[2] = {0.0};
  double JBz_con_n[2] = {0.0};
  double Jphi_n[2] = {0.0};
  double Jpsi_n[2] = {0.0};
  double Ex_n[2] = {0.0};
  double Ey_n[2] = {0.0};
  double Ez_n[2] = {0.0};
  double Hx_n[2] = {0.0};
  double Hy_n[2] = {0.0};
  double Hz_n[2] = {0.0};
  double FD11_n[2] = {0.0};
  double FD12_n[2] = {0.0};
  double FD13_n[2] = {0.0};
  double FD21_n[2] = {0.0};
  double FD22_n[2] = {0.0};
  double FD23_n[2] = {0.0};
  double FD31_n[2] = {0.0};
  double FD32_n[2] = {0.0};
  double FD33_n[2] = {0.0};
  double FB11_n[2] = {0.0};
  double FB12_n[2] = {0.0};
  double FB13_n[2] = {0.0};
  double FB21_n[2] = {0.0};
  double FB22_n[2] = {0.0};
  double FB23_n[2] = {0.0};
  double FB31_n[2] = {0.0};
  double FB32_n[2] = {0.0};
  double FB33_n[2] = {0.0};
  double FPhi1_n[2] = {0.0};
  double FPhi2_n[2] = {0.0};
  double FPhi3_n[2] = {0.0};
  double FPsi1_n[2] = {0.0};
  double FPsi2_n[2] = {0.0};
  double FPsi3_n[2] = {0.0};
  
  const double *JDx = &fields[0]; 
  const double *JDy = &fields[2]; 
  const double *JDz = &fields[4]; 
  const double *JBx = &fields[6]; 
  const double *JBy = &fields[8]; 
  const double *JBz = &fields[10]; 
  const double *Jphi = &fields[12]; 
  const double *Jpsi = &fields[14]; 
  
  double *outJDx = &out[0]; 
  double *outJDy = &out[2]; 
  double *outJDz = &out[4]; 
  double *outJBx = &out[6]; 
  double *outJBy = &out[8]; 
  double *outJBz = &out[10]; 
  double *outJphi = &out[12]; 
  double *outJpsi = &out[14]; 
  
  JDx_con_n[0] = 0.7071067811865475*JDx[0]-0.7071067811865475*JDx[1];
  JDy_con_n[0] = 0.7071067811865475*JDy[0]-0.7071067811865475*JDy[1];
  JDz_con_n[0] = 0.7071067811865475*JDz[0]-0.7071067811865475*JDz[1];
  JBx_con_n[0] = 0.7071067811865475*JBx[0]-0.7071067811865475*JBx[1];
  JBy_con_n[0] = 0.7071067811865475*JBy[0]-0.7071067811865475*JBy[1];
  JBz_con_n[0] = 0.7071067811865475*JBz[0]-0.7071067811865475*JBz[1];
  Jphi_n[0] = 0.7071067811865475*Jphi[0]-0.7071067811865475*Jphi[1];
  Jpsi_n[0] = 0.7071067811865475*Jpsi[0]-0.7071067811865475*Jpsi[1];
  
  JDx_con_n[1] = 0.7071067811865475*JDx[1]+0.7071067811865475*JDx[0];
  JDy_con_n[1] = 0.7071067811865475*JDy[1]+0.7071067811865475*JDy[0];
  JDz_con_n[1] = 0.7071067811865475*JDz[1]+0.7071067811865475*JDz[0];
  JBx_con_n[1] = 0.7071067811865475*JBx[1]+0.7071067811865475*JBx[0];
  JBy_con_n[1] = 0.7071067811865475*JBy[1]+0.7071067811865475*JBy[0];
  JBz_con_n[1] = 0.7071067811865475*JBz[1]+0.7071067811865475*JBz[0];
  Jphi_n[1] = 0.7071067811865475*Jphi[1]+0.7071067811865475*Jphi[0];
  Jpsi_n[1] = 0.7071067811865475*Jpsi[1]+0.7071067811865475*Jpsi[0];
  
  const double *h_xx_nodal = &h_ij_nodal[0]; 
  const double *h_xy_nodal = &h_ij_nodal[2]; 
  const double *h_xz_nodal = &h_ij_nodal[4]; 
  const double *h_yy_nodal = &h_ij_nodal[6]; 
  const double *h_yz_nodal = &h_ij_nodal[8]; 
  const double *h_zz_nodal = &h_ij_nodal[10]; 
  const double *h_xx_inv_nodal = &h_ij_inv_nodal[0]; 
  const double *h_xy_inv_nodal = &h_ij_inv_nodal[2]; 
  const double *h_xz_inv_nodal = &h_ij_inv_nodal[4]; 
  const double *h_yy_inv_nodal = &h_ij_inv_nodal[6]; 
  const double *h_yz_inv_nodal = &h_ij_inv_nodal[8]; 
  const double *h_zz_inv_nodal = &h_ij_inv_nodal[10]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[2]; 
  const double *shift_nodal_z = &shift_nodal[4]; 
  
  for (int i=0; i<2; ++i) {
    Hx_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*JBx_con_n[i] + h_xy_nodal[i]*JBy_con_n[i] + h_xz_nodal[i]*JBz_con_n[i] ) / J_c[i] - ( shift_nodal_y[i]*JDz_con_n[i] - shift_nodal_z[i]*JDy_con_n[i]); 
    Hy_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*JBx_con_n[i] + h_yy_nodal[i]*JBy_con_n[i] + h_yz_nodal[i]*JBz_con_n[i] ) / J_c[i] - ( - shift_nodal_x[i]*JDz_con_n[i] + shift_nodal_z[i]*JDx_con_n[i]); 
    Hz_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*JBx_con_n[i] + h_yz_nodal[i]*JBy_con_n[i] + h_zz_nodal[i]*JBz_con_n[i] ) / J_c[i] - ( shift_nodal_x[i]*JDy_con_n[i] - shift_nodal_y[i]*JDx_con_n[i]); 
    Ex_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*JDx_con_n[i] + h_xy_nodal[i]*JDy_con_n[i] + h_xz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( shift_nodal_y[i]*JBz_con_n[i] - shift_nodal_z[i]*JBy_con_n[i]); 
    Ey_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*JDx_con_n[i] + h_yy_nodal[i]*JDy_con_n[i] + h_yz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( - shift_nodal_x[i]*JBz_con_n[i] + shift_nodal_z[i]*JBx_con_n[i]); 
    Ez_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*JDx_con_n[i] + h_yz_nodal[i]*JDy_con_n[i] + h_zz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( shift_nodal_x[i]*JBy_con_n[i] - shift_nodal_y[i]*JBx_con_n[i]); 
    FD11_n[i] = chi*h_xx_inv_nodal[i]*Jphi_n[i]; 
    FD12_n[i] = chi*h_xy_inv_nodal[i]*Jphi_n[i]; 
    FD13_n[i] = chi*h_xz_inv_nodal[i]*Jphi_n[i]; 
    FD21_n[i] = chi*h_xy_inv_nodal[i]*Jphi_n[i]; 
    FD22_n[i] = chi*h_yy_inv_nodal[i]*Jphi_n[i]; 
    FD23_n[i] = chi*h_yz_inv_nodal[i]*Jphi_n[i]; 
    FD31_n[i] = chi*h_xz_inv_nodal[i]*Jphi_n[i]; 
    FD32_n[i] = chi*h_yz_inv_nodal[i]*Jphi_n[i]; 
    FD33_n[i] = chi*h_zz_inv_nodal[i]*Jphi_n[i]; 
    FB11_n[i] = gamma*h_xx_inv_nodal[i]*Jpsi_n[i]; 
    FB12_n[i] = gamma*h_xy_inv_nodal[i]*Jpsi_n[i]; 
    FB13_n[i] = gamma*h_xz_inv_nodal[i]*Jpsi_n[i]; 
    FB21_n[i] = gamma*h_xy_inv_nodal[i]*Jpsi_n[i]; 
    FB22_n[i] = gamma*h_yy_inv_nodal[i]*Jpsi_n[i]; 
    FB23_n[i] = gamma*h_yz_inv_nodal[i]*Jpsi_n[i]; 
    FB31_n[i] = gamma*h_xz_inv_nodal[i]*Jpsi_n[i]; 
    FB32_n[i] = gamma*h_yz_inv_nodal[i]*Jpsi_n[i]; 
    FB33_n[i] = gamma*h_zz_inv_nodal[i]*Jpsi_n[i]; 
    FPhi1_n[i] = chi*JDx_con_n[i]; 
    FPhi2_n[i] = chi*JDy_con_n[i]; 
    FPhi3_n[i] = chi*JDz_con_n[i]; 
    FPsi1_n[i] = gamma*JBx_con_n[i]; 
    FPsi2_n[i] = gamma*JBy_con_n[i]; 
    FPsi3_n[i] = gamma*JBz_con_n[i]; 
  }
  
  const double jacob_pos_x_inv = 1.0/jacob_pos[0]; 
  outJDx[1] += 1.224744871391589*FD11_n[1]*dx0*jacob_pos_x_inv+1.224744871391589*FD11_n[0]*dx0*jacob_pos_x_inv; 

  outJDy[1] += 1.224744871391589*Hz_n[1]*dx0*jacob_pos_x_inv+1.224744871391589*FD12_n[1]*dx0*jacob_pos_x_inv+1.224744871391589*Hz_n[0]*dx0*jacob_pos_x_inv+1.224744871391589*FD12_n[0]*dx0*jacob_pos_x_inv; 

  outJDz[1] += -(1.224744871391589*Hy_n[1]*dx0*jacob_pos_x_inv)+1.224744871391589*FD13_n[1]*dx0*jacob_pos_x_inv-1.224744871391589*Hy_n[0]*dx0*jacob_pos_x_inv+1.224744871391589*FD13_n[0]*dx0*jacob_pos_x_inv; 

  outJBx[1] += 1.224744871391589*FB11_n[1]*dx0*jacob_pos_x_inv+1.224744871391589*FB11_n[0]*dx0*jacob_pos_x_inv; 

  outJBy[1] += 1.224744871391589*FB12_n[1]*dx0*jacob_pos_x_inv-1.224744871391589*Ez_n[1]*dx0*jacob_pos_x_inv+1.224744871391589*FB12_n[0]*dx0*jacob_pos_x_inv-1.224744871391589*Ez_n[0]*dx0*jacob_pos_x_inv; 

  outJBz[1] += 1.224744871391589*FB13_n[1]*dx0*jacob_pos_x_inv+1.224744871391589*Ey_n[1]*dx0*jacob_pos_x_inv+1.224744871391589*FB13_n[0]*dx0*jacob_pos_x_inv+1.224744871391589*Ey_n[0]*dx0*jacob_pos_x_inv; 

  outJphi[1] += 1.224744871391589*FPhi1_n[1]*dx0*jacob_pos_x_inv+1.224744871391589*FPhi1_n[0]*dx0*jacob_pos_x_inv; 

  outJpsi[1] += 1.224744871391589*FPsi1_n[1]*dx0*jacob_pos_x_inv+1.224744871391589*FPsi1_n[0]*dx0*jacob_pos_x_inv; 

} 
