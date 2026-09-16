#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_1x3v_ser_p1(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal, 
    const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force) 
{ 
  const double chi = meq->chi, gamma = meq->gamma; 
  const double jacob_pos_tot_inv = 1.0/(jacob_pos[0]); 
  double dx0 = 2.0/dx[0]; 
  
  double JDx_con_n[2] = {0.0};
  double JDy_con_n[2] = {0.0};
  double JDz_con_n[2] = {0.0};
  double JBx_con_n[2] = {0.0};
  double JBy_con_n[2] = {0.0};
  double JBz_con_n[2] = {0.0};
  double Ex_n[2] = {0.0};
  double Ey_n[2] = {0.0};
  double Ez_n[2] = {0.0};
  double Bx_n[2] = {0.0};
  double By_n[2] = {0.0};
  double Bz_n[2] = {0.0};
  double E_conf_for_force_n[6] = {0.0};
  double B_conf_for_force_n[6] = {0.0};
  
  const double *JDx = &fields[0]; 
  const double *JDy = &fields[2]; 
  const double *JDz = &fields[4]; 
  const double *JBx = &fields[6]; 
  const double *JBy = &fields[8]; 
  const double *JBz = &fields[10]; 
  
  JDx_con_n[0] = 0.7071067811865475*JDx[0]-0.7071067811865475*JDx[1];
  JDy_con_n[0] = 0.7071067811865475*JDy[0]-0.7071067811865475*JDy[1];
  JDz_con_n[0] = 0.7071067811865475*JDz[0]-0.7071067811865475*JDz[1];
  JBx_con_n[0] = 0.7071067811865475*JBx[0]-0.7071067811865475*JBx[1];
  JBy_con_n[0] = 0.7071067811865475*JBy[0]-0.7071067811865475*JBy[1];
  JBz_con_n[0] = 0.7071067811865475*JBz[0]-0.7071067811865475*JBz[1];
  
  JDx_con_n[1] = 0.7071067811865475*JDx[1]+0.7071067811865475*JDx[0];
  JDy_con_n[1] = 0.7071067811865475*JDy[1]+0.7071067811865475*JDy[0];
  JDz_con_n[1] = 0.7071067811865475*JDz[1]+0.7071067811865475*JDz[0];
  JBx_con_n[1] = 0.7071067811865475*JBx[1]+0.7071067811865475*JBx[0];
  JBy_con_n[1] = 0.7071067811865475*JBy[1]+0.7071067811865475*JBy[0];
  JBz_con_n[1] = 0.7071067811865475*JBz[1]+0.7071067811865475*JBz[0];
  
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
  const double *vierb_xx_cov_nodal = &vierb_cov_nodal[0]; 
  const double *vierb_xy_cov_nodal = &vierb_cov_nodal[2]; 
  const double *vierb_xz_cov_nodal = &vierb_cov_nodal[4]; 
  const double *vierb_yx_cov_nodal = &vierb_cov_nodal[6]; 
  const double *vierb_yy_cov_nodal = &vierb_cov_nodal[8]; 
  const double *vierb_yz_cov_nodal = &vierb_cov_nodal[10]; 
  const double *vierb_zx_cov_nodal = &vierb_cov_nodal[12]; 
  const double *vierb_zy_cov_nodal = &vierb_cov_nodal[14]; 
  const double *vierb_zz_cov_nodal = &vierb_cov_nodal[16]; 
  const double *vierb_xx_con_nodal = &vierb_con_nodal[0]; 
  const double *vierb_xy_con_nodal = &vierb_con_nodal[2]; 
  const double *vierb_xz_con_nodal = &vierb_con_nodal[4]; 
  const double *vierb_yx_con_nodal = &vierb_con_nodal[6]; 
  const double *vierb_yy_con_nodal = &vierb_con_nodal[8]; 
  const double *vierb_yz_con_nodal = &vierb_con_nodal[10]; 
  const double *vierb_zx_con_nodal = &vierb_con_nodal[12]; 
  const double *vierb_zy_con_nodal = &vierb_con_nodal[14]; 
  const double *vierb_zz_con_nodal = &vierb_con_nodal[16]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[2]; 
  const double *shift_nodal_z = &shift_nodal[4]; 
  
  for (int i=0; i<2; ++i) {
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
  
  for (int i=0; i<2; ++i) {
    E_conf_for_force_n[i + 0*2] = vierb_xx_con_nodal[i]*Ex_n[i] + vierb_yx_con_nodal[i]*Ey_n[i]  + vierb_zx_con_nodal[i]*Ez_n[i];
    E_conf_for_force_n[i + 1*2] = vierb_xy_con_nodal[i]*Ex_n[i] + vierb_yy_con_nodal[i]*Ey_n[i]  + vierb_zy_con_nodal[i]*Ez_n[i];
    E_conf_for_force_n[i + 2*2] = vierb_xz_con_nodal[i]*Ex_n[i] + vierb_yz_con_nodal[i]*Ey_n[i]  + vierb_zz_con_nodal[i]*Ez_n[i];
    B_conf_for_force_n[i + 0*2] = vierb_xx_cov_nodal[i]*Bx_n[i] + vierb_yx_cov_nodal[i]*By_n[i] + vierb_zx_cov_nodal[i]*Bz_n[i];
    B_conf_for_force_n[i + 1*2] = vierb_xy_cov_nodal[i]*Bx_n[i] + vierb_yy_cov_nodal[i]*By_n[i] + vierb_zy_cov_nodal[i]*Bz_n[i];
    B_conf_for_force_n[i + 2*2] = vierb_xz_cov_nodal[i]*Bx_n[i] + vierb_yz_cov_nodal[i]*By_n[i] + vierb_zz_cov_nodal[i]*Bz_n[i];
  }
  
  double *E_conf_for_force_comp;
  double *B_conf_for_force_comp;
  E_conf_for_force_comp = &E_conf_for_force[0];
  E_conf_for_force_comp[0] = 0.7071067811865475*E_conf_for_force_n[1]+0.7071067811865475*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[1] = 0.7071067811865475*E_conf_for_force_n[1]-0.7071067811865475*E_conf_for_force_n[0]; 

  E_conf_for_force_comp = &E_conf_for_force[2];
  E_conf_for_force_comp[0] = 0.7071067811865475*E_conf_for_force_n[3]+0.7071067811865475*E_conf_for_force_n[2]; 
  E_conf_for_force_comp[1] = 0.7071067811865475*E_conf_for_force_n[3]-0.7071067811865475*E_conf_for_force_n[2]; 

  E_conf_for_force_comp = &E_conf_for_force[4];
  E_conf_for_force_comp[0] = 0.7071067811865475*E_conf_for_force_n[5]+0.7071067811865475*E_conf_for_force_n[4]; 
  E_conf_for_force_comp[1] = 0.7071067811865475*E_conf_for_force_n[5]-0.7071067811865475*E_conf_for_force_n[4]; 

  B_conf_for_force_comp = &B_conf_for_force[0];
  B_conf_for_force_comp[0] = 0.7071067811865475*B_conf_for_force_n[1]+0.7071067811865475*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[1] = 0.7071067811865475*B_conf_for_force_n[1]-0.7071067811865475*B_conf_for_force_n[0]; 

  B_conf_for_force_comp = &B_conf_for_force[2];
  B_conf_for_force_comp[0] = 0.7071067811865475*B_conf_for_force_n[3]+0.7071067811865475*B_conf_for_force_n[2]; 
  B_conf_for_force_comp[1] = 0.7071067811865475*B_conf_for_force_n[3]-0.7071067811865475*B_conf_for_force_n[2]; 

  B_conf_for_force_comp = &B_conf_for_force[4];
  B_conf_for_force_comp[0] = 0.7071067811865475*B_conf_for_force_n[5]+0.7071067811865475*B_conf_for_force_n[4]; 
  B_conf_for_force_comp[1] = 0.7071067811865475*B_conf_for_force_n[5]-0.7071067811865475*B_conf_for_force_n[4]; 

} 
