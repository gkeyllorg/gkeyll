#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_3x3v_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal, 
    const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force) 
{ 
  const double chi = meq->chi, gamma = meq->gamma; 
  const double jacob_pos_tot_inv = 1.0/(jacob_pos[0]*jacob_pos[2]*jacob_pos[4]); 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 
  double dx2 = 2.0/dx[2]; 
  
  double JDx_con_n[8] = {0.0};
  double JDy_con_n[8] = {0.0};
  double JDz_con_n[8] = {0.0};
  double JBx_con_n[8] = {0.0};
  double JBy_con_n[8] = {0.0};
  double JBz_con_n[8] = {0.0};
  double Ex_n[8] = {0.0};
  double Ey_n[8] = {0.0};
  double Ez_n[8] = {0.0};
  double Bx_n[8] = {0.0};
  double By_n[8] = {0.0};
  double Bz_n[8] = {0.0};
  double E_conf_for_force_n[24] = {0.0};
  double B_conf_for_force_n[24] = {0.0};
  
  const double *JDx = &fields[0]; 
  const double *JDy = &fields[8]; 
  const double *JDz = &fields[16]; 
  const double *JBx = &fields[24]; 
  const double *JBy = &fields[32]; 
  const double *JBz = &fields[40]; 
  
  JDx_con_n[0] = -(0.3535533905932737*JDx[7])+0.3535533905932737*JDx[6]+0.3535533905932737*JDx[5]+0.3535533905932737*JDx[4]-0.3535533905932737*JDx[3]-0.3535533905932737*JDx[2]-0.3535533905932737*JDx[1]+0.3535533905932737*JDx[0];
  JDy_con_n[0] = -(0.3535533905932737*JDy[7])+0.3535533905932737*JDy[6]+0.3535533905932737*JDy[5]+0.3535533905932737*JDy[4]-0.3535533905932737*JDy[3]-0.3535533905932737*JDy[2]-0.3535533905932737*JDy[1]+0.3535533905932737*JDy[0];
  JDz_con_n[0] = -(0.3535533905932737*JDz[7])+0.3535533905932737*JDz[6]+0.3535533905932737*JDz[5]+0.3535533905932737*JDz[4]-0.3535533905932737*JDz[3]-0.3535533905932737*JDz[2]-0.3535533905932737*JDz[1]+0.3535533905932737*JDz[0];
  JBx_con_n[0] = -(0.3535533905932737*JBx[7])+0.3535533905932737*JBx[6]+0.3535533905932737*JBx[5]+0.3535533905932737*JBx[4]-0.3535533905932737*JBx[3]-0.3535533905932737*JBx[2]-0.3535533905932737*JBx[1]+0.3535533905932737*JBx[0];
  JBy_con_n[0] = -(0.3535533905932737*JBy[7])+0.3535533905932737*JBy[6]+0.3535533905932737*JBy[5]+0.3535533905932737*JBy[4]-0.3535533905932737*JBy[3]-0.3535533905932737*JBy[2]-0.3535533905932737*JBy[1]+0.3535533905932737*JBy[0];
  JBz_con_n[0] = -(0.3535533905932737*JBz[7])+0.3535533905932737*JBz[6]+0.3535533905932737*JBz[5]+0.3535533905932737*JBz[4]-0.3535533905932737*JBz[3]-0.3535533905932737*JBz[2]-0.3535533905932737*JBz[1]+0.3535533905932737*JBz[0];
  
  JDx_con_n[1] = 0.3535533905932737*JDx[7]-0.3535533905932737*JDx[6]-0.3535533905932737*JDx[5]+0.3535533905932737*JDx[4]+0.3535533905932737*JDx[3]-0.3535533905932737*JDx[2]-0.3535533905932737*JDx[1]+0.3535533905932737*JDx[0];
  JDy_con_n[1] = 0.3535533905932737*JDy[7]-0.3535533905932737*JDy[6]-0.3535533905932737*JDy[5]+0.3535533905932737*JDy[4]+0.3535533905932737*JDy[3]-0.3535533905932737*JDy[2]-0.3535533905932737*JDy[1]+0.3535533905932737*JDy[0];
  JDz_con_n[1] = 0.3535533905932737*JDz[7]-0.3535533905932737*JDz[6]-0.3535533905932737*JDz[5]+0.3535533905932737*JDz[4]+0.3535533905932737*JDz[3]-0.3535533905932737*JDz[2]-0.3535533905932737*JDz[1]+0.3535533905932737*JDz[0];
  JBx_con_n[1] = 0.3535533905932737*JBx[7]-0.3535533905932737*JBx[6]-0.3535533905932737*JBx[5]+0.3535533905932737*JBx[4]+0.3535533905932737*JBx[3]-0.3535533905932737*JBx[2]-0.3535533905932737*JBx[1]+0.3535533905932737*JBx[0];
  JBy_con_n[1] = 0.3535533905932737*JBy[7]-0.3535533905932737*JBy[6]-0.3535533905932737*JBy[5]+0.3535533905932737*JBy[4]+0.3535533905932737*JBy[3]-0.3535533905932737*JBy[2]-0.3535533905932737*JBy[1]+0.3535533905932737*JBy[0];
  JBz_con_n[1] = 0.3535533905932737*JBz[7]-0.3535533905932737*JBz[6]-0.3535533905932737*JBz[5]+0.3535533905932737*JBz[4]+0.3535533905932737*JBz[3]-0.3535533905932737*JBz[2]-0.3535533905932737*JBz[1]+0.3535533905932737*JBz[0];
  
  JDx_con_n[2] = 0.3535533905932737*JDx[7]-0.3535533905932737*JDx[6]+0.3535533905932737*JDx[5]-0.3535533905932737*JDx[4]-0.3535533905932737*JDx[3]+0.3535533905932737*JDx[2]-0.3535533905932737*JDx[1]+0.3535533905932737*JDx[0];
  JDy_con_n[2] = 0.3535533905932737*JDy[7]-0.3535533905932737*JDy[6]+0.3535533905932737*JDy[5]-0.3535533905932737*JDy[4]-0.3535533905932737*JDy[3]+0.3535533905932737*JDy[2]-0.3535533905932737*JDy[1]+0.3535533905932737*JDy[0];
  JDz_con_n[2] = 0.3535533905932737*JDz[7]-0.3535533905932737*JDz[6]+0.3535533905932737*JDz[5]-0.3535533905932737*JDz[4]-0.3535533905932737*JDz[3]+0.3535533905932737*JDz[2]-0.3535533905932737*JDz[1]+0.3535533905932737*JDz[0];
  JBx_con_n[2] = 0.3535533905932737*JBx[7]-0.3535533905932737*JBx[6]+0.3535533905932737*JBx[5]-0.3535533905932737*JBx[4]-0.3535533905932737*JBx[3]+0.3535533905932737*JBx[2]-0.3535533905932737*JBx[1]+0.3535533905932737*JBx[0];
  JBy_con_n[2] = 0.3535533905932737*JBy[7]-0.3535533905932737*JBy[6]+0.3535533905932737*JBy[5]-0.3535533905932737*JBy[4]-0.3535533905932737*JBy[3]+0.3535533905932737*JBy[2]-0.3535533905932737*JBy[1]+0.3535533905932737*JBy[0];
  JBz_con_n[2] = 0.3535533905932737*JBz[7]-0.3535533905932737*JBz[6]+0.3535533905932737*JBz[5]-0.3535533905932737*JBz[4]-0.3535533905932737*JBz[3]+0.3535533905932737*JBz[2]-0.3535533905932737*JBz[1]+0.3535533905932737*JBz[0];
  
  JDx_con_n[3] = -(0.3535533905932737*JDx[7])+0.3535533905932737*JDx[6]-0.3535533905932737*JDx[5]-0.3535533905932737*JDx[4]+0.3535533905932737*JDx[3]+0.3535533905932737*JDx[2]-0.3535533905932737*JDx[1]+0.3535533905932737*JDx[0];
  JDy_con_n[3] = -(0.3535533905932737*JDy[7])+0.3535533905932737*JDy[6]-0.3535533905932737*JDy[5]-0.3535533905932737*JDy[4]+0.3535533905932737*JDy[3]+0.3535533905932737*JDy[2]-0.3535533905932737*JDy[1]+0.3535533905932737*JDy[0];
  JDz_con_n[3] = -(0.3535533905932737*JDz[7])+0.3535533905932737*JDz[6]-0.3535533905932737*JDz[5]-0.3535533905932737*JDz[4]+0.3535533905932737*JDz[3]+0.3535533905932737*JDz[2]-0.3535533905932737*JDz[1]+0.3535533905932737*JDz[0];
  JBx_con_n[3] = -(0.3535533905932737*JBx[7])+0.3535533905932737*JBx[6]-0.3535533905932737*JBx[5]-0.3535533905932737*JBx[4]+0.3535533905932737*JBx[3]+0.3535533905932737*JBx[2]-0.3535533905932737*JBx[1]+0.3535533905932737*JBx[0];
  JBy_con_n[3] = -(0.3535533905932737*JBy[7])+0.3535533905932737*JBy[6]-0.3535533905932737*JBy[5]-0.3535533905932737*JBy[4]+0.3535533905932737*JBy[3]+0.3535533905932737*JBy[2]-0.3535533905932737*JBy[1]+0.3535533905932737*JBy[0];
  JBz_con_n[3] = -(0.3535533905932737*JBz[7])+0.3535533905932737*JBz[6]-0.3535533905932737*JBz[5]-0.3535533905932737*JBz[4]+0.3535533905932737*JBz[3]+0.3535533905932737*JBz[2]-0.3535533905932737*JBz[1]+0.3535533905932737*JBz[0];
  
  JDx_con_n[4] = 0.3535533905932737*JDx[7]+0.3535533905932737*JDx[6]-0.3535533905932737*JDx[5]-0.3535533905932737*JDx[4]-0.3535533905932737*JDx[3]-0.3535533905932737*JDx[2]+0.3535533905932737*JDx[1]+0.3535533905932737*JDx[0];
  JDy_con_n[4] = 0.3535533905932737*JDy[7]+0.3535533905932737*JDy[6]-0.3535533905932737*JDy[5]-0.3535533905932737*JDy[4]-0.3535533905932737*JDy[3]-0.3535533905932737*JDy[2]+0.3535533905932737*JDy[1]+0.3535533905932737*JDy[0];
  JDz_con_n[4] = 0.3535533905932737*JDz[7]+0.3535533905932737*JDz[6]-0.3535533905932737*JDz[5]-0.3535533905932737*JDz[4]-0.3535533905932737*JDz[3]-0.3535533905932737*JDz[2]+0.3535533905932737*JDz[1]+0.3535533905932737*JDz[0];
  JBx_con_n[4] = 0.3535533905932737*JBx[7]+0.3535533905932737*JBx[6]-0.3535533905932737*JBx[5]-0.3535533905932737*JBx[4]-0.3535533905932737*JBx[3]-0.3535533905932737*JBx[2]+0.3535533905932737*JBx[1]+0.3535533905932737*JBx[0];
  JBy_con_n[4] = 0.3535533905932737*JBy[7]+0.3535533905932737*JBy[6]-0.3535533905932737*JBy[5]-0.3535533905932737*JBy[4]-0.3535533905932737*JBy[3]-0.3535533905932737*JBy[2]+0.3535533905932737*JBy[1]+0.3535533905932737*JBy[0];
  JBz_con_n[4] = 0.3535533905932737*JBz[7]+0.3535533905932737*JBz[6]-0.3535533905932737*JBz[5]-0.3535533905932737*JBz[4]-0.3535533905932737*JBz[3]-0.3535533905932737*JBz[2]+0.3535533905932737*JBz[1]+0.3535533905932737*JBz[0];
  
  JDx_con_n[5] = -(0.3535533905932737*JDx[7])-0.3535533905932737*JDx[6]+0.3535533905932737*JDx[5]-0.3535533905932737*JDx[4]+0.3535533905932737*JDx[3]-0.3535533905932737*JDx[2]+0.3535533905932737*JDx[1]+0.3535533905932737*JDx[0];
  JDy_con_n[5] = -(0.3535533905932737*JDy[7])-0.3535533905932737*JDy[6]+0.3535533905932737*JDy[5]-0.3535533905932737*JDy[4]+0.3535533905932737*JDy[3]-0.3535533905932737*JDy[2]+0.3535533905932737*JDy[1]+0.3535533905932737*JDy[0];
  JDz_con_n[5] = -(0.3535533905932737*JDz[7])-0.3535533905932737*JDz[6]+0.3535533905932737*JDz[5]-0.3535533905932737*JDz[4]+0.3535533905932737*JDz[3]-0.3535533905932737*JDz[2]+0.3535533905932737*JDz[1]+0.3535533905932737*JDz[0];
  JBx_con_n[5] = -(0.3535533905932737*JBx[7])-0.3535533905932737*JBx[6]+0.3535533905932737*JBx[5]-0.3535533905932737*JBx[4]+0.3535533905932737*JBx[3]-0.3535533905932737*JBx[2]+0.3535533905932737*JBx[1]+0.3535533905932737*JBx[0];
  JBy_con_n[5] = -(0.3535533905932737*JBy[7])-0.3535533905932737*JBy[6]+0.3535533905932737*JBy[5]-0.3535533905932737*JBy[4]+0.3535533905932737*JBy[3]-0.3535533905932737*JBy[2]+0.3535533905932737*JBy[1]+0.3535533905932737*JBy[0];
  JBz_con_n[5] = -(0.3535533905932737*JBz[7])-0.3535533905932737*JBz[6]+0.3535533905932737*JBz[5]-0.3535533905932737*JBz[4]+0.3535533905932737*JBz[3]-0.3535533905932737*JBz[2]+0.3535533905932737*JBz[1]+0.3535533905932737*JBz[0];
  
  JDx_con_n[6] = -(0.3535533905932737*JDx[7])-0.3535533905932737*JDx[6]-0.3535533905932737*JDx[5]+0.3535533905932737*JDx[4]-0.3535533905932737*JDx[3]+0.3535533905932737*JDx[2]+0.3535533905932737*JDx[1]+0.3535533905932737*JDx[0];
  JDy_con_n[6] = -(0.3535533905932737*JDy[7])-0.3535533905932737*JDy[6]-0.3535533905932737*JDy[5]+0.3535533905932737*JDy[4]-0.3535533905932737*JDy[3]+0.3535533905932737*JDy[2]+0.3535533905932737*JDy[1]+0.3535533905932737*JDy[0];
  JDz_con_n[6] = -(0.3535533905932737*JDz[7])-0.3535533905932737*JDz[6]-0.3535533905932737*JDz[5]+0.3535533905932737*JDz[4]-0.3535533905932737*JDz[3]+0.3535533905932737*JDz[2]+0.3535533905932737*JDz[1]+0.3535533905932737*JDz[0];
  JBx_con_n[6] = -(0.3535533905932737*JBx[7])-0.3535533905932737*JBx[6]-0.3535533905932737*JBx[5]+0.3535533905932737*JBx[4]-0.3535533905932737*JBx[3]+0.3535533905932737*JBx[2]+0.3535533905932737*JBx[1]+0.3535533905932737*JBx[0];
  JBy_con_n[6] = -(0.3535533905932737*JBy[7])-0.3535533905932737*JBy[6]-0.3535533905932737*JBy[5]+0.3535533905932737*JBy[4]-0.3535533905932737*JBy[3]+0.3535533905932737*JBy[2]+0.3535533905932737*JBy[1]+0.3535533905932737*JBy[0];
  JBz_con_n[6] = -(0.3535533905932737*JBz[7])-0.3535533905932737*JBz[6]-0.3535533905932737*JBz[5]+0.3535533905932737*JBz[4]-0.3535533905932737*JBz[3]+0.3535533905932737*JBz[2]+0.3535533905932737*JBz[1]+0.3535533905932737*JBz[0];
  
  JDx_con_n[7] = 0.3535533905932737*JDx[7]+0.3535533905932737*JDx[6]+0.3535533905932737*JDx[5]+0.3535533905932737*JDx[4]+0.3535533905932737*JDx[3]+0.3535533905932737*JDx[2]+0.3535533905932737*JDx[1]+0.3535533905932737*JDx[0];
  JDy_con_n[7] = 0.3535533905932737*JDy[7]+0.3535533905932737*JDy[6]+0.3535533905932737*JDy[5]+0.3535533905932737*JDy[4]+0.3535533905932737*JDy[3]+0.3535533905932737*JDy[2]+0.3535533905932737*JDy[1]+0.3535533905932737*JDy[0];
  JDz_con_n[7] = 0.3535533905932737*JDz[7]+0.3535533905932737*JDz[6]+0.3535533905932737*JDz[5]+0.3535533905932737*JDz[4]+0.3535533905932737*JDz[3]+0.3535533905932737*JDz[2]+0.3535533905932737*JDz[1]+0.3535533905932737*JDz[0];
  JBx_con_n[7] = 0.3535533905932737*JBx[7]+0.3535533905932737*JBx[6]+0.3535533905932737*JBx[5]+0.3535533905932737*JBx[4]+0.3535533905932737*JBx[3]+0.3535533905932737*JBx[2]+0.3535533905932737*JBx[1]+0.3535533905932737*JBx[0];
  JBy_con_n[7] = 0.3535533905932737*JBy[7]+0.3535533905932737*JBy[6]+0.3535533905932737*JBy[5]+0.3535533905932737*JBy[4]+0.3535533905932737*JBy[3]+0.3535533905932737*JBy[2]+0.3535533905932737*JBy[1]+0.3535533905932737*JBy[0];
  JBz_con_n[7] = 0.3535533905932737*JBz[7]+0.3535533905932737*JBz[6]+0.3535533905932737*JBz[5]+0.3535533905932737*JBz[4]+0.3535533905932737*JBz[3]+0.3535533905932737*JBz[2]+0.3535533905932737*JBz[1]+0.3535533905932737*JBz[0];
  
  const double *h_xx_nodal = &h_ij_nodal[0]; 
  const double *h_xy_nodal = &h_ij_nodal[8]; 
  const double *h_xz_nodal = &h_ij_nodal[16]; 
  const double *h_yy_nodal = &h_ij_nodal[24]; 
  const double *h_yz_nodal = &h_ij_nodal[32]; 
  const double *h_zz_nodal = &h_ij_nodal[40]; 
  const double *h_xx_inv_nodal = &h_ij_inv_nodal[0]; 
  const double *h_xy_inv_nodal = &h_ij_inv_nodal[8]; 
  const double *h_xz_inv_nodal = &h_ij_inv_nodal[16]; 
  const double *h_yy_inv_nodal = &h_ij_inv_nodal[24]; 
  const double *h_yz_inv_nodal = &h_ij_inv_nodal[32]; 
  const double *h_zz_inv_nodal = &h_ij_inv_nodal[40]; 
  const double *vierb_xx_cov_nodal = &vierb_cov_nodal[0]; 
  const double *vierb_xy_cov_nodal = &vierb_cov_nodal[8]; 
  const double *vierb_xz_cov_nodal = &vierb_cov_nodal[16]; 
  const double *vierb_yx_cov_nodal = &vierb_cov_nodal[24]; 
  const double *vierb_yy_cov_nodal = &vierb_cov_nodal[32]; 
  const double *vierb_yz_cov_nodal = &vierb_cov_nodal[40]; 
  const double *vierb_zx_cov_nodal = &vierb_cov_nodal[48]; 
  const double *vierb_zy_cov_nodal = &vierb_cov_nodal[56]; 
  const double *vierb_zz_cov_nodal = &vierb_cov_nodal[64]; 
  const double *vierb_xx_con_nodal = &vierb_con_nodal[0]; 
  const double *vierb_xy_con_nodal = &vierb_con_nodal[8]; 
  const double *vierb_xz_con_nodal = &vierb_con_nodal[16]; 
  const double *vierb_yx_con_nodal = &vierb_con_nodal[24]; 
  const double *vierb_yy_con_nodal = &vierb_con_nodal[32]; 
  const double *vierb_yz_con_nodal = &vierb_con_nodal[40]; 
  const double *vierb_zx_con_nodal = &vierb_con_nodal[48]; 
  const double *vierb_zy_con_nodal = &vierb_con_nodal[56]; 
  const double *vierb_zz_con_nodal = &vierb_con_nodal[64]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[8]; 
  const double *shift_nodal_z = &shift_nodal[16]; 
  
  for (int i=0; i<8; ++i) {
    Ex_n[i] = lapse_nodal[i]*( h_xx_nodal[i]*JDx_con_n[i] + h_xy_nodal[i]*JDy_con_n[i] + h_xz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( shift_nodal_y[i]*JBz_con_n[i] - shift_nodal_z[i]*JBy_con_n[i]); 
    Ey_n[i] = lapse_nodal[i]*( h_xy_nodal[i]*JDx_con_n[i] + h_yy_nodal[i]*JDy_con_n[i] + h_yz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( - shift_nodal_x[i]*JBz_con_n[i] + shift_nodal_z[i]*JBx_con_n[i]); 
    Ez_n[i] = lapse_nodal[i]*( h_xz_nodal[i]*JDx_con_n[i] + h_yz_nodal[i]*JDy_con_n[i] + h_zz_nodal[i]*JDz_con_n[i] ) / J_c[i] + ( shift_nodal_x[i]*JBy_con_n[i] - shift_nodal_y[i]*JBx_con_n[i]); 
    Bx_n[i] = JBx_con_n[i] / J_c[i]; 
    By_n[i] = JBy_con_n[i] / J_c[i]; 
    Bz_n[i] = JBz_con_n[i] / J_c[i]; 
    Ex_n[i] *= jacob_pos_tot_inv; 
    Ey_n[i] *= jacob_pos_tot_inv; 
    Ez_n[i] *= jacob_pos_tot_inv; 
    Bx_n[i] *= jacob_pos_tot_inv; 
    By_n[i] *= jacob_pos_tot_inv; 
    Bz_n[i] *= jacob_pos_tot_inv; 
  }
  
  for (int i=0; i<8; ++i) {
    E_conf_for_force_n[i + 0*8] = vierb_xx_con_nodal[i]*Ex_n[i] + vierb_yx_con_nodal[i]*Ey_n[i]  + vierb_zx_con_nodal[i]*Ez_n[i];
    E_conf_for_force_n[i + 1*8] = vierb_xy_con_nodal[i]*Ex_n[i] + vierb_yy_con_nodal[i]*Ey_n[i]  + vierb_zy_con_nodal[i]*Ez_n[i];
    E_conf_for_force_n[i + 2*8] = vierb_xz_con_nodal[i]*Ex_n[i] + vierb_yz_con_nodal[i]*Ey_n[i]  + vierb_zz_con_nodal[i]*Ez_n[i];
    B_conf_for_force_n[i + 0*8] = vierb_xx_cov_nodal[i]*Bx_n[i] + vierb_yx_cov_nodal[i]*By_n[i] + vierb_zx_cov_nodal[i]*Bz_n[i];
    B_conf_for_force_n[i + 1*8] = vierb_xy_cov_nodal[i]*Bx_n[i] + vierb_yy_cov_nodal[i]*By_n[i] + vierb_zy_cov_nodal[i]*Bz_n[i];
    B_conf_for_force_n[i + 2*8] = vierb_xz_cov_nodal[i]*Bx_n[i] + vierb_yz_cov_nodal[i]*By_n[i] + vierb_zz_cov_nodal[i]*Bz_n[i];
  }
  
  double *E_conf_for_force_comp;
  double *B_conf_for_force_comp;
  E_conf_for_force_comp = &E_conf_for_force[0];
  E_conf_for_force_comp[0] = 0.3535533905932737*E_conf_for_force_n[7]+0.3535533905932737*E_conf_for_force_n[6]+0.3535533905932737*E_conf_for_force_n[5]+0.3535533905932737*E_conf_for_force_n[4]+0.3535533905932737*E_conf_for_force_n[3]+0.3535533905932737*E_conf_for_force_n[2]+0.3535533905932737*E_conf_for_force_n[1]+0.3535533905932737*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[1] = 0.3535533905932737*E_conf_for_force_n[7]+0.3535533905932737*E_conf_for_force_n[6]+0.3535533905932737*E_conf_for_force_n[5]+0.3535533905932737*E_conf_for_force_n[4]-0.3535533905932737*E_conf_for_force_n[3]-0.3535533905932737*E_conf_for_force_n[2]-0.3535533905932737*E_conf_for_force_n[1]-0.3535533905932737*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[2] = 0.3535533905932737*E_conf_for_force_n[7]+0.3535533905932737*E_conf_for_force_n[6]-0.3535533905932737*E_conf_for_force_n[5]-0.3535533905932737*E_conf_for_force_n[4]+0.3535533905932737*E_conf_for_force_n[3]+0.3535533905932737*E_conf_for_force_n[2]-0.3535533905932737*E_conf_for_force_n[1]-0.3535533905932737*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[3] = 0.3535533905932737*E_conf_for_force_n[7]-0.3535533905932737*E_conf_for_force_n[6]+0.3535533905932737*E_conf_for_force_n[5]-0.3535533905932737*E_conf_for_force_n[4]+0.3535533905932737*E_conf_for_force_n[3]-0.3535533905932737*E_conf_for_force_n[2]+0.3535533905932737*E_conf_for_force_n[1]-0.3535533905932737*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[4] = 0.3535533905932737*E_conf_for_force_n[7]+0.3535533905932737*E_conf_for_force_n[6]-0.3535533905932737*E_conf_for_force_n[5]-0.3535533905932737*E_conf_for_force_n[4]-0.3535533905932737*E_conf_for_force_n[3]-0.3535533905932737*E_conf_for_force_n[2]+0.3535533905932737*E_conf_for_force_n[1]+0.3535533905932737*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[5] = 0.3535533905932737*E_conf_for_force_n[7]-0.3535533905932737*E_conf_for_force_n[6]+0.3535533905932737*E_conf_for_force_n[5]-0.3535533905932737*E_conf_for_force_n[4]-0.3535533905932737*E_conf_for_force_n[3]+0.3535533905932737*E_conf_for_force_n[2]-0.3535533905932737*E_conf_for_force_n[1]+0.3535533905932737*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[6] = 0.3535533905932737*E_conf_for_force_n[7]-0.3535533905932737*E_conf_for_force_n[6]-0.3535533905932737*E_conf_for_force_n[5]+0.3535533905932737*E_conf_for_force_n[4]+0.3535533905932737*E_conf_for_force_n[3]-0.3535533905932737*E_conf_for_force_n[2]-0.3535533905932737*E_conf_for_force_n[1]+0.3535533905932737*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[7] = 0.3535533905932737*E_conf_for_force_n[7]-0.3535533905932737*E_conf_for_force_n[6]-0.3535533905932737*E_conf_for_force_n[5]+0.3535533905932737*E_conf_for_force_n[4]-0.3535533905932737*E_conf_for_force_n[3]+0.3535533905932737*E_conf_for_force_n[2]+0.3535533905932737*E_conf_for_force_n[1]-0.3535533905932737*E_conf_for_force_n[0]; 

  E_conf_for_force_comp = &E_conf_for_force[8];
  E_conf_for_force_comp[0] = 0.3535533905932737*E_conf_for_force_n[15]+0.3535533905932737*E_conf_for_force_n[14]+0.3535533905932737*E_conf_for_force_n[13]+0.3535533905932737*E_conf_for_force_n[12]+0.3535533905932737*E_conf_for_force_n[11]+0.3535533905932737*E_conf_for_force_n[10]+0.3535533905932737*E_conf_for_force_n[9]+0.3535533905932737*E_conf_for_force_n[8]; 
  E_conf_for_force_comp[1] = 0.3535533905932737*E_conf_for_force_n[15]+0.3535533905932737*E_conf_for_force_n[14]+0.3535533905932737*E_conf_for_force_n[13]+0.3535533905932737*E_conf_for_force_n[12]-0.3535533905932737*E_conf_for_force_n[11]-0.3535533905932737*E_conf_for_force_n[10]-0.3535533905932737*E_conf_for_force_n[9]-0.3535533905932737*E_conf_for_force_n[8]; 
  E_conf_for_force_comp[2] = 0.3535533905932737*E_conf_for_force_n[15]+0.3535533905932737*E_conf_for_force_n[14]-0.3535533905932737*E_conf_for_force_n[13]-0.3535533905932737*E_conf_for_force_n[12]+0.3535533905932737*E_conf_for_force_n[11]+0.3535533905932737*E_conf_for_force_n[10]-0.3535533905932737*E_conf_for_force_n[9]-0.3535533905932737*E_conf_for_force_n[8]; 
  E_conf_for_force_comp[3] = 0.3535533905932737*E_conf_for_force_n[15]-0.3535533905932737*E_conf_for_force_n[14]+0.3535533905932737*E_conf_for_force_n[13]-0.3535533905932737*E_conf_for_force_n[12]+0.3535533905932737*E_conf_for_force_n[11]-0.3535533905932737*E_conf_for_force_n[10]+0.3535533905932737*E_conf_for_force_n[9]-0.3535533905932737*E_conf_for_force_n[8]; 
  E_conf_for_force_comp[4] = 0.3535533905932737*E_conf_for_force_n[15]+0.3535533905932737*E_conf_for_force_n[14]-0.3535533905932737*E_conf_for_force_n[13]-0.3535533905932737*E_conf_for_force_n[12]-0.3535533905932737*E_conf_for_force_n[11]-0.3535533905932737*E_conf_for_force_n[10]+0.3535533905932737*E_conf_for_force_n[9]+0.3535533905932737*E_conf_for_force_n[8]; 
  E_conf_for_force_comp[5] = 0.3535533905932737*E_conf_for_force_n[15]-0.3535533905932737*E_conf_for_force_n[14]+0.3535533905932737*E_conf_for_force_n[13]-0.3535533905932737*E_conf_for_force_n[12]-0.3535533905932737*E_conf_for_force_n[11]+0.3535533905932737*E_conf_for_force_n[10]-0.3535533905932737*E_conf_for_force_n[9]+0.3535533905932737*E_conf_for_force_n[8]; 
  E_conf_for_force_comp[6] = 0.3535533905932737*E_conf_for_force_n[15]-0.3535533905932737*E_conf_for_force_n[14]-0.3535533905932737*E_conf_for_force_n[13]+0.3535533905932737*E_conf_for_force_n[12]+0.3535533905932737*E_conf_for_force_n[11]-0.3535533905932737*E_conf_for_force_n[10]-0.3535533905932737*E_conf_for_force_n[9]+0.3535533905932737*E_conf_for_force_n[8]; 
  E_conf_for_force_comp[7] = 0.3535533905932737*E_conf_for_force_n[15]-0.3535533905932737*E_conf_for_force_n[14]-0.3535533905932737*E_conf_for_force_n[13]+0.3535533905932737*E_conf_for_force_n[12]-0.3535533905932737*E_conf_for_force_n[11]+0.3535533905932737*E_conf_for_force_n[10]+0.3535533905932737*E_conf_for_force_n[9]-0.3535533905932737*E_conf_for_force_n[8]; 

  E_conf_for_force_comp = &E_conf_for_force[16];
  E_conf_for_force_comp[0] = 0.3535533905932737*E_conf_for_force_n[23]+0.3535533905932737*E_conf_for_force_n[22]+0.3535533905932737*E_conf_for_force_n[21]+0.3535533905932737*E_conf_for_force_n[20]+0.3535533905932737*E_conf_for_force_n[19]+0.3535533905932737*E_conf_for_force_n[18]+0.3535533905932737*E_conf_for_force_n[17]+0.3535533905932737*E_conf_for_force_n[16]; 
  E_conf_for_force_comp[1] = 0.3535533905932737*E_conf_for_force_n[23]+0.3535533905932737*E_conf_for_force_n[22]+0.3535533905932737*E_conf_for_force_n[21]+0.3535533905932737*E_conf_for_force_n[20]-0.3535533905932737*E_conf_for_force_n[19]-0.3535533905932737*E_conf_for_force_n[18]-0.3535533905932737*E_conf_for_force_n[17]-0.3535533905932737*E_conf_for_force_n[16]; 
  E_conf_for_force_comp[2] = 0.3535533905932737*E_conf_for_force_n[23]+0.3535533905932737*E_conf_for_force_n[22]-0.3535533905932737*E_conf_for_force_n[21]-0.3535533905932737*E_conf_for_force_n[20]+0.3535533905932737*E_conf_for_force_n[19]+0.3535533905932737*E_conf_for_force_n[18]-0.3535533905932737*E_conf_for_force_n[17]-0.3535533905932737*E_conf_for_force_n[16]; 
  E_conf_for_force_comp[3] = 0.3535533905932737*E_conf_for_force_n[23]-0.3535533905932737*E_conf_for_force_n[22]+0.3535533905932737*E_conf_for_force_n[21]-0.3535533905932737*E_conf_for_force_n[20]+0.3535533905932737*E_conf_for_force_n[19]-0.3535533905932737*E_conf_for_force_n[18]+0.3535533905932737*E_conf_for_force_n[17]-0.3535533905932737*E_conf_for_force_n[16]; 
  E_conf_for_force_comp[4] = 0.3535533905932737*E_conf_for_force_n[23]+0.3535533905932737*E_conf_for_force_n[22]-0.3535533905932737*E_conf_for_force_n[21]-0.3535533905932737*E_conf_for_force_n[20]-0.3535533905932737*E_conf_for_force_n[19]-0.3535533905932737*E_conf_for_force_n[18]+0.3535533905932737*E_conf_for_force_n[17]+0.3535533905932737*E_conf_for_force_n[16]; 
  E_conf_for_force_comp[5] = 0.3535533905932737*E_conf_for_force_n[23]-0.3535533905932737*E_conf_for_force_n[22]+0.3535533905932737*E_conf_for_force_n[21]-0.3535533905932737*E_conf_for_force_n[20]-0.3535533905932737*E_conf_for_force_n[19]+0.3535533905932737*E_conf_for_force_n[18]-0.3535533905932737*E_conf_for_force_n[17]+0.3535533905932737*E_conf_for_force_n[16]; 
  E_conf_for_force_comp[6] = 0.3535533905932737*E_conf_for_force_n[23]-0.3535533905932737*E_conf_for_force_n[22]-0.3535533905932737*E_conf_for_force_n[21]+0.3535533905932737*E_conf_for_force_n[20]+0.3535533905932737*E_conf_for_force_n[19]-0.3535533905932737*E_conf_for_force_n[18]-0.3535533905932737*E_conf_for_force_n[17]+0.3535533905932737*E_conf_for_force_n[16]; 
  E_conf_for_force_comp[7] = 0.3535533905932737*E_conf_for_force_n[23]-0.3535533905932737*E_conf_for_force_n[22]-0.3535533905932737*E_conf_for_force_n[21]+0.3535533905932737*E_conf_for_force_n[20]-0.3535533905932737*E_conf_for_force_n[19]+0.3535533905932737*E_conf_for_force_n[18]+0.3535533905932737*E_conf_for_force_n[17]-0.3535533905932737*E_conf_for_force_n[16]; 

  B_conf_for_force_comp = &B_conf_for_force[0];
  B_conf_for_force_comp[0] = 0.3535533905932737*B_conf_for_force_n[7]+0.3535533905932737*B_conf_for_force_n[6]+0.3535533905932737*B_conf_for_force_n[5]+0.3535533905932737*B_conf_for_force_n[4]+0.3535533905932737*B_conf_for_force_n[3]+0.3535533905932737*B_conf_for_force_n[2]+0.3535533905932737*B_conf_for_force_n[1]+0.3535533905932737*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[1] = 0.3535533905932737*B_conf_for_force_n[7]+0.3535533905932737*B_conf_for_force_n[6]+0.3535533905932737*B_conf_for_force_n[5]+0.3535533905932737*B_conf_for_force_n[4]-0.3535533905932737*B_conf_for_force_n[3]-0.3535533905932737*B_conf_for_force_n[2]-0.3535533905932737*B_conf_for_force_n[1]-0.3535533905932737*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[2] = 0.3535533905932737*B_conf_for_force_n[7]+0.3535533905932737*B_conf_for_force_n[6]-0.3535533905932737*B_conf_for_force_n[5]-0.3535533905932737*B_conf_for_force_n[4]+0.3535533905932737*B_conf_for_force_n[3]+0.3535533905932737*B_conf_for_force_n[2]-0.3535533905932737*B_conf_for_force_n[1]-0.3535533905932737*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[3] = 0.3535533905932737*B_conf_for_force_n[7]-0.3535533905932737*B_conf_for_force_n[6]+0.3535533905932737*B_conf_for_force_n[5]-0.3535533905932737*B_conf_for_force_n[4]+0.3535533905932737*B_conf_for_force_n[3]-0.3535533905932737*B_conf_for_force_n[2]+0.3535533905932737*B_conf_for_force_n[1]-0.3535533905932737*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[4] = 0.3535533905932737*B_conf_for_force_n[7]+0.3535533905932737*B_conf_for_force_n[6]-0.3535533905932737*B_conf_for_force_n[5]-0.3535533905932737*B_conf_for_force_n[4]-0.3535533905932737*B_conf_for_force_n[3]-0.3535533905932737*B_conf_for_force_n[2]+0.3535533905932737*B_conf_for_force_n[1]+0.3535533905932737*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[5] = 0.3535533905932737*B_conf_for_force_n[7]-0.3535533905932737*B_conf_for_force_n[6]+0.3535533905932737*B_conf_for_force_n[5]-0.3535533905932737*B_conf_for_force_n[4]-0.3535533905932737*B_conf_for_force_n[3]+0.3535533905932737*B_conf_for_force_n[2]-0.3535533905932737*B_conf_for_force_n[1]+0.3535533905932737*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[6] = 0.3535533905932737*B_conf_for_force_n[7]-0.3535533905932737*B_conf_for_force_n[6]-0.3535533905932737*B_conf_for_force_n[5]+0.3535533905932737*B_conf_for_force_n[4]+0.3535533905932737*B_conf_for_force_n[3]-0.3535533905932737*B_conf_for_force_n[2]-0.3535533905932737*B_conf_for_force_n[1]+0.3535533905932737*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[7] = 0.3535533905932737*B_conf_for_force_n[7]-0.3535533905932737*B_conf_for_force_n[6]-0.3535533905932737*B_conf_for_force_n[5]+0.3535533905932737*B_conf_for_force_n[4]-0.3535533905932737*B_conf_for_force_n[3]+0.3535533905932737*B_conf_for_force_n[2]+0.3535533905932737*B_conf_for_force_n[1]-0.3535533905932737*B_conf_for_force_n[0]; 

  B_conf_for_force_comp = &B_conf_for_force[8];
  B_conf_for_force_comp[0] = 0.3535533905932737*B_conf_for_force_n[15]+0.3535533905932737*B_conf_for_force_n[14]+0.3535533905932737*B_conf_for_force_n[13]+0.3535533905932737*B_conf_for_force_n[12]+0.3535533905932737*B_conf_for_force_n[11]+0.3535533905932737*B_conf_for_force_n[10]+0.3535533905932737*B_conf_for_force_n[9]+0.3535533905932737*B_conf_for_force_n[8]; 
  B_conf_for_force_comp[1] = 0.3535533905932737*B_conf_for_force_n[15]+0.3535533905932737*B_conf_for_force_n[14]+0.3535533905932737*B_conf_for_force_n[13]+0.3535533905932737*B_conf_for_force_n[12]-0.3535533905932737*B_conf_for_force_n[11]-0.3535533905932737*B_conf_for_force_n[10]-0.3535533905932737*B_conf_for_force_n[9]-0.3535533905932737*B_conf_for_force_n[8]; 
  B_conf_for_force_comp[2] = 0.3535533905932737*B_conf_for_force_n[15]+0.3535533905932737*B_conf_for_force_n[14]-0.3535533905932737*B_conf_for_force_n[13]-0.3535533905932737*B_conf_for_force_n[12]+0.3535533905932737*B_conf_for_force_n[11]+0.3535533905932737*B_conf_for_force_n[10]-0.3535533905932737*B_conf_for_force_n[9]-0.3535533905932737*B_conf_for_force_n[8]; 
  B_conf_for_force_comp[3] = 0.3535533905932737*B_conf_for_force_n[15]-0.3535533905932737*B_conf_for_force_n[14]+0.3535533905932737*B_conf_for_force_n[13]-0.3535533905932737*B_conf_for_force_n[12]+0.3535533905932737*B_conf_for_force_n[11]-0.3535533905932737*B_conf_for_force_n[10]+0.3535533905932737*B_conf_for_force_n[9]-0.3535533905932737*B_conf_for_force_n[8]; 
  B_conf_for_force_comp[4] = 0.3535533905932737*B_conf_for_force_n[15]+0.3535533905932737*B_conf_for_force_n[14]-0.3535533905932737*B_conf_for_force_n[13]-0.3535533905932737*B_conf_for_force_n[12]-0.3535533905932737*B_conf_for_force_n[11]-0.3535533905932737*B_conf_for_force_n[10]+0.3535533905932737*B_conf_for_force_n[9]+0.3535533905932737*B_conf_for_force_n[8]; 
  B_conf_for_force_comp[5] = 0.3535533905932737*B_conf_for_force_n[15]-0.3535533905932737*B_conf_for_force_n[14]+0.3535533905932737*B_conf_for_force_n[13]-0.3535533905932737*B_conf_for_force_n[12]-0.3535533905932737*B_conf_for_force_n[11]+0.3535533905932737*B_conf_for_force_n[10]-0.3535533905932737*B_conf_for_force_n[9]+0.3535533905932737*B_conf_for_force_n[8]; 
  B_conf_for_force_comp[6] = 0.3535533905932737*B_conf_for_force_n[15]-0.3535533905932737*B_conf_for_force_n[14]-0.3535533905932737*B_conf_for_force_n[13]+0.3535533905932737*B_conf_for_force_n[12]+0.3535533905932737*B_conf_for_force_n[11]-0.3535533905932737*B_conf_for_force_n[10]-0.3535533905932737*B_conf_for_force_n[9]+0.3535533905932737*B_conf_for_force_n[8]; 
  B_conf_for_force_comp[7] = 0.3535533905932737*B_conf_for_force_n[15]-0.3535533905932737*B_conf_for_force_n[14]-0.3535533905932737*B_conf_for_force_n[13]+0.3535533905932737*B_conf_for_force_n[12]-0.3535533905932737*B_conf_for_force_n[11]+0.3535533905932737*B_conf_for_force_n[10]+0.3535533905932737*B_conf_for_force_n[9]-0.3535533905932737*B_conf_for_force_n[8]; 

  B_conf_for_force_comp = &B_conf_for_force[16];
  B_conf_for_force_comp[0] = 0.3535533905932737*B_conf_for_force_n[23]+0.3535533905932737*B_conf_for_force_n[22]+0.3535533905932737*B_conf_for_force_n[21]+0.3535533905932737*B_conf_for_force_n[20]+0.3535533905932737*B_conf_for_force_n[19]+0.3535533905932737*B_conf_for_force_n[18]+0.3535533905932737*B_conf_for_force_n[17]+0.3535533905932737*B_conf_for_force_n[16]; 
  B_conf_for_force_comp[1] = 0.3535533905932737*B_conf_for_force_n[23]+0.3535533905932737*B_conf_for_force_n[22]+0.3535533905932737*B_conf_for_force_n[21]+0.3535533905932737*B_conf_for_force_n[20]-0.3535533905932737*B_conf_for_force_n[19]-0.3535533905932737*B_conf_for_force_n[18]-0.3535533905932737*B_conf_for_force_n[17]-0.3535533905932737*B_conf_for_force_n[16]; 
  B_conf_for_force_comp[2] = 0.3535533905932737*B_conf_for_force_n[23]+0.3535533905932737*B_conf_for_force_n[22]-0.3535533905932737*B_conf_for_force_n[21]-0.3535533905932737*B_conf_for_force_n[20]+0.3535533905932737*B_conf_for_force_n[19]+0.3535533905932737*B_conf_for_force_n[18]-0.3535533905932737*B_conf_for_force_n[17]-0.3535533905932737*B_conf_for_force_n[16]; 
  B_conf_for_force_comp[3] = 0.3535533905932737*B_conf_for_force_n[23]-0.3535533905932737*B_conf_for_force_n[22]+0.3535533905932737*B_conf_for_force_n[21]-0.3535533905932737*B_conf_for_force_n[20]+0.3535533905932737*B_conf_for_force_n[19]-0.3535533905932737*B_conf_for_force_n[18]+0.3535533905932737*B_conf_for_force_n[17]-0.3535533905932737*B_conf_for_force_n[16]; 
  B_conf_for_force_comp[4] = 0.3535533905932737*B_conf_for_force_n[23]+0.3535533905932737*B_conf_for_force_n[22]-0.3535533905932737*B_conf_for_force_n[21]-0.3535533905932737*B_conf_for_force_n[20]-0.3535533905932737*B_conf_for_force_n[19]-0.3535533905932737*B_conf_for_force_n[18]+0.3535533905932737*B_conf_for_force_n[17]+0.3535533905932737*B_conf_for_force_n[16]; 
  B_conf_for_force_comp[5] = 0.3535533905932737*B_conf_for_force_n[23]-0.3535533905932737*B_conf_for_force_n[22]+0.3535533905932737*B_conf_for_force_n[21]-0.3535533905932737*B_conf_for_force_n[20]-0.3535533905932737*B_conf_for_force_n[19]+0.3535533905932737*B_conf_for_force_n[18]-0.3535533905932737*B_conf_for_force_n[17]+0.3535533905932737*B_conf_for_force_n[16]; 
  B_conf_for_force_comp[6] = 0.3535533905932737*B_conf_for_force_n[23]-0.3535533905932737*B_conf_for_force_n[22]-0.3535533905932737*B_conf_for_force_n[21]+0.3535533905932737*B_conf_for_force_n[20]+0.3535533905932737*B_conf_for_force_n[19]-0.3535533905932737*B_conf_for_force_n[18]-0.3535533905932737*B_conf_for_force_n[17]+0.3535533905932737*B_conf_for_force_n[16]; 
  B_conf_for_force_comp[7] = 0.3535533905932737*B_conf_for_force_n[23]-0.3535533905932737*B_conf_for_force_n[22]-0.3535533905932737*B_conf_for_force_n[21]+0.3535533905932737*B_conf_for_force_n[20]-0.3535533905932737*B_conf_for_force_n[19]+0.3535533905932737*B_conf_for_force_n[18]+0.3535533905932737*B_conf_for_force_n[17]-0.3535533905932737*B_conf_for_force_n[16]; 

} 
