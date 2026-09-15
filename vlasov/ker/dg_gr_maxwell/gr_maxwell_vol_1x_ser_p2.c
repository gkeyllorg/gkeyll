#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void gr_maxwell_vol_1x_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, double* GKYL_RESTRICT out) 
{ 
  const double chi = meq->chi, gamma = meq->gamma; 
  double dx0 = 2.0/dx[0]; 
  
  double JDx_con_n[3] = {0.0};
  double JDy_con_n[3] = {0.0};
  double JDz_con_n[3] = {0.0};
  double JBx_con_n[3] = {0.0};
  double JBy_con_n[3] = {0.0};
  double JBz_con_n[3] = {0.0};
  double Jphi_n[3] = {0.0};
  double Jpsi_n[3] = {0.0};
  double Ex_n[3] = {0.0};
  double Ey_n[3] = {0.0};
  double Ez_n[3] = {0.0};
  double Hx_n[3] = {0.0};
  double Hy_n[3] = {0.0};
  double Hz_n[3] = {0.0};
  double FD11_n[3] = {0.0};
  double FD12_n[3] = {0.0};
  double FD13_n[3] = {0.0};
  double FD21_n[3] = {0.0};
  double FD22_n[3] = {0.0};
  double FD23_n[3] = {0.0};
  double FD31_n[3] = {0.0};
  double FD32_n[3] = {0.0};
  double FD33_n[3] = {0.0};
  double FB11_n[3] = {0.0};
  double FB12_n[3] = {0.0};
  double FB13_n[3] = {0.0};
  double FB21_n[3] = {0.0};
  double FB22_n[3] = {0.0};
  double FB23_n[3] = {0.0};
  double FB31_n[3] = {0.0};
  double FB32_n[3] = {0.0};
  double FB33_n[3] = {0.0};
  double FPhi1_n[3] = {0.0};
  double FPhi2_n[3] = {0.0};
  double FPhi3_n[3] = {0.0};
  double FPsi1_n[3] = {0.0};
  double FPsi2_n[3] = {0.0};
  double FPsi3_n[3] = {0.0};
  
  const double *JDx = &fields[0]; 
  const double *JDy = &fields[3]; 
  const double *JDz = &fields[6]; 
  const double *JBx = &fields[9]; 
  const double *JBy = &fields[12]; 
  const double *JBz = &fields[15]; 
  const double *Jphi = &fields[18]; 
  const double *Jpsi = &fields[21]; 
  
  double *outJDx = &out[0]; 
  double *outJDy = &out[3]; 
  double *outJDz = &out[6]; 
  double *outJBx = &out[9]; 
  double *outJBy = &out[12]; 
  double *outJBz = &out[15]; 
  double *outJphi = &out[18]; 
  double *outJpsi = &out[21]; 
  
  JDx_con_n[0] = 0.6324555320336759*JDx[2]-0.9486832980505137*JDx[1]+0.7071067811865475*JDx[0];
  JDy_con_n[0] = 0.6324555320336759*JDy[2]-0.9486832980505137*JDy[1]+0.7071067811865475*JDy[0];
  JDz_con_n[0] = 0.6324555320336759*JDz[2]-0.9486832980505137*JDz[1]+0.7071067811865475*JDz[0];
  JBx_con_n[0] = 0.6324555320336759*JBx[2]-0.9486832980505137*JBx[1]+0.7071067811865475*JBx[0];
  JBy_con_n[0] = 0.6324555320336759*JBy[2]-0.9486832980505137*JBy[1]+0.7071067811865475*JBy[0];
  JBz_con_n[0] = 0.6324555320336759*JBz[2]-0.9486832980505137*JBz[1]+0.7071067811865475*JBz[0];
  Jphi_n[0] = 0.6324555320336759*Jphi[2]-0.9486832980505137*Jphi[1]+0.7071067811865475*Jphi[0];
  Jpsi_n[0] = 0.6324555320336759*Jpsi[2]-0.9486832980505137*Jpsi[1]+0.7071067811865475*Jpsi[0];
  
  JDx_con_n[1] = 0.7071067811865475*JDx[0]-0.7905694150420947*JDx[2];
  JDy_con_n[1] = 0.7071067811865475*JDy[0]-0.7905694150420947*JDy[2];
  JDz_con_n[1] = 0.7071067811865475*JDz[0]-0.7905694150420947*JDz[2];
  JBx_con_n[1] = 0.7071067811865475*JBx[0]-0.7905694150420947*JBx[2];
  JBy_con_n[1] = 0.7071067811865475*JBy[0]-0.7905694150420947*JBy[2];
  JBz_con_n[1] = 0.7071067811865475*JBz[0]-0.7905694150420947*JBz[2];
  Jphi_n[1] = 0.7071067811865475*Jphi[0]-0.7905694150420947*Jphi[2];
  Jpsi_n[1] = 0.7071067811865475*Jpsi[0]-0.7905694150420947*Jpsi[2];
  
  JDx_con_n[2] = 0.6324555320336759*JDx[2]+0.9486832980505137*JDx[1]+0.7071067811865475*JDx[0];
  JDy_con_n[2] = 0.6324555320336759*JDy[2]+0.9486832980505137*JDy[1]+0.7071067811865475*JDy[0];
  JDz_con_n[2] = 0.6324555320336759*JDz[2]+0.9486832980505137*JDz[1]+0.7071067811865475*JDz[0];
  JBx_con_n[2] = 0.6324555320336759*JBx[2]+0.9486832980505137*JBx[1]+0.7071067811865475*JBx[0];
  JBy_con_n[2] = 0.6324555320336759*JBy[2]+0.9486832980505137*JBy[1]+0.7071067811865475*JBy[0];
  JBz_con_n[2] = 0.6324555320336759*JBz[2]+0.9486832980505137*JBz[1]+0.7071067811865475*JBz[0];
  Jphi_n[2] = 0.6324555320336759*Jphi[2]+0.9486832980505137*Jphi[1]+0.7071067811865475*Jphi[0];
  Jpsi_n[2] = 0.6324555320336759*Jpsi[2]+0.9486832980505137*Jpsi[1]+0.7071067811865475*Jpsi[0];
  
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
  
  for (int i=0; i<3; ++i) {
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
  outJDx[1] += 0.6804138174397718*FD11_n[2]*dx0*jacob_pos_x_inv+1.0886621079036354*FD11_n[1]*dx0*jacob_pos_x_inv+0.6804138174397718*FD11_n[0]*dx0*jacob_pos_x_inv; 
  outJDx[2] += 2.041241452319315*FD11_n[2]*dx0*jacob_pos_x_inv-2.041241452319315*FD11_n[0]*dx0*jacob_pos_x_inv; 

  outJDy[1] += 0.6804138174397718*Hz_n[2]*dx0*jacob_pos_x_inv+0.6804138174397718*FD12_n[2]*dx0*jacob_pos_x_inv+1.0886621079036354*Hz_n[1]*dx0*jacob_pos_x_inv+1.0886621079036354*FD12_n[1]*dx0*jacob_pos_x_inv+0.6804138174397718*Hz_n[0]*dx0*jacob_pos_x_inv+0.6804138174397718*FD12_n[0]*dx0*jacob_pos_x_inv; 
  outJDy[2] += 2.041241452319315*Hz_n[2]*dx0*jacob_pos_x_inv+2.041241452319315*FD12_n[2]*dx0*jacob_pos_x_inv-2.041241452319315*Hz_n[0]*dx0*jacob_pos_x_inv-2.041241452319315*FD12_n[0]*dx0*jacob_pos_x_inv; 

  outJDz[1] += -(0.6804138174397718*Hy_n[2]*dx0*jacob_pos_x_inv)+0.6804138174397718*FD13_n[2]*dx0*jacob_pos_x_inv-1.0886621079036354*Hy_n[1]*dx0*jacob_pos_x_inv+1.0886621079036354*FD13_n[1]*dx0*jacob_pos_x_inv-0.6804138174397718*Hy_n[0]*dx0*jacob_pos_x_inv+0.6804138174397718*FD13_n[0]*dx0*jacob_pos_x_inv; 
  outJDz[2] += -(2.041241452319315*Hy_n[2]*dx0*jacob_pos_x_inv)+2.041241452319315*FD13_n[2]*dx0*jacob_pos_x_inv+2.041241452319315*Hy_n[0]*dx0*jacob_pos_x_inv-2.041241452319315*FD13_n[0]*dx0*jacob_pos_x_inv; 

  outJBx[1] += 0.6804138174397718*FB11_n[2]*dx0*jacob_pos_x_inv+1.0886621079036354*FB11_n[1]*dx0*jacob_pos_x_inv+0.6804138174397718*FB11_n[0]*dx0*jacob_pos_x_inv; 
  outJBx[2] += 2.041241452319315*FB11_n[2]*dx0*jacob_pos_x_inv-2.041241452319315*FB11_n[0]*dx0*jacob_pos_x_inv; 

  outJBy[1] += 0.6804138174397718*FB12_n[2]*dx0*jacob_pos_x_inv-0.6804138174397718*Ez_n[2]*dx0*jacob_pos_x_inv+1.0886621079036354*FB12_n[1]*dx0*jacob_pos_x_inv-1.0886621079036354*Ez_n[1]*dx0*jacob_pos_x_inv+0.6804138174397718*FB12_n[0]*dx0*jacob_pos_x_inv-0.6804138174397718*Ez_n[0]*dx0*jacob_pos_x_inv; 
  outJBy[2] += 2.041241452319315*FB12_n[2]*dx0*jacob_pos_x_inv-2.041241452319315*Ez_n[2]*dx0*jacob_pos_x_inv-2.041241452319315*FB12_n[0]*dx0*jacob_pos_x_inv+2.041241452319315*Ez_n[0]*dx0*jacob_pos_x_inv; 

  outJBz[1] += 0.6804138174397718*FB13_n[2]*dx0*jacob_pos_x_inv+0.6804138174397718*Ey_n[2]*dx0*jacob_pos_x_inv+1.0886621079036354*FB13_n[1]*dx0*jacob_pos_x_inv+1.0886621079036354*Ey_n[1]*dx0*jacob_pos_x_inv+0.6804138174397718*FB13_n[0]*dx0*jacob_pos_x_inv+0.6804138174397718*Ey_n[0]*dx0*jacob_pos_x_inv; 
  outJBz[2] += 2.041241452319315*FB13_n[2]*dx0*jacob_pos_x_inv+2.041241452319315*Ey_n[2]*dx0*jacob_pos_x_inv-2.041241452319315*FB13_n[0]*dx0*jacob_pos_x_inv-2.041241452319315*Ey_n[0]*dx0*jacob_pos_x_inv; 

  outJphi[1] += 0.6804138174397718*FPhi1_n[2]*dx0*jacob_pos_x_inv+1.0886621079036354*FPhi1_n[1]*dx0*jacob_pos_x_inv+0.6804138174397718*FPhi1_n[0]*dx0*jacob_pos_x_inv; 
  outJphi[2] += 2.041241452319315*FPhi1_n[2]*dx0*jacob_pos_x_inv-2.041241452319315*FPhi1_n[0]*dx0*jacob_pos_x_inv; 

  outJpsi[1] += 0.6804138174397718*FPsi1_n[2]*dx0*jacob_pos_x_inv+1.0886621079036354*FPsi1_n[1]*dx0*jacob_pos_x_inv+0.6804138174397718*FPsi1_n[0]*dx0*jacob_pos_x_inv; 
  outJpsi[2] += 2.041241452319315*FPsi1_n[2]*dx0*jacob_pos_x_inv-2.041241452319315*FPsi1_n[0]*dx0*jacob_pos_x_inv; 

} 
