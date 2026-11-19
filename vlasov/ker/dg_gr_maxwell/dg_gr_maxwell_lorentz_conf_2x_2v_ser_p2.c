#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_lorentz_conf_2x2v_ser_p2(const gkyl_dg_gr_maxwell_inp *meq, const double *dx, const double *jacob_pos, const double *lapse_nodal, const double *shift_nodal,
  const double *h_ij_nodal, const double *h_ij_inv_nodal, const double *J_c, const double *fields, const double *vierb_cov_nodal, 
    const double *vierb_con_nodal, double* GKYL_RESTRICT E_conf_for_force, double* GKYL_RESTRICT B_conf_for_force) 
{ 
  const double chi = meq->chi, gamma = meq->gamma; 
  const double jacob_pos_tot_inv = 1.0/(jacob_pos[0]*jacob_pos[3]); 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 
  
  double JDx_con_n[9] = {0.0};
  double JDy_con_n[9] = {0.0};
  double JDz_con_n[9] = {0.0};
  double JBx_con_n[9] = {0.0};
  double JBy_con_n[9] = {0.0};
  double JBz_con_n[9] = {0.0};
  double Ex_n[9] = {0.0};
  double Ey_n[9] = {0.0};
  double Ez_n[9] = {0.0};
  double Bx_n[9] = {0.0};
  double By_n[9] = {0.0};
  double Bz_n[9] = {0.0};
  double E_conf_for_force_n[27] = {0.0};
  double B_conf_for_force_n[27] = {0.0};
  
  const double *JDx = &fields[0]; 
  const double *JDy = &fields[8]; 
  const double *JDz = &fields[16]; 
  const double *JBx = &fields[24]; 
  const double *JBy = &fields[32]; 
  const double *JBz = &fields[40]; 
  
  JDx_con_n[0] = -(0.5999999999999995*JDx[7])-0.5999999999999999*JDx[6]+0.4472135954999579*JDx[5]+0.4472135954999579*JDx[4]+0.9*JDx[3]-0.6708203932499369*JDx[2]-0.6708203932499369*JDx[1]+0.5*JDx[0];
  JDy_con_n[0] = -(0.5999999999999995*JDy[7])-0.5999999999999999*JDy[6]+0.4472135954999579*JDy[5]+0.4472135954999579*JDy[4]+0.9*JDy[3]-0.6708203932499369*JDy[2]-0.6708203932499369*JDy[1]+0.5*JDy[0];
  JDz_con_n[0] = -(0.5999999999999995*JDz[7])-0.5999999999999999*JDz[6]+0.4472135954999579*JDz[5]+0.4472135954999579*JDz[4]+0.9*JDz[3]-0.6708203932499369*JDz[2]-0.6708203932499369*JDz[1]+0.5*JDz[0];
  JBx_con_n[0] = -(0.5999999999999995*JBx[7])-0.5999999999999999*JBx[6]+0.4472135954999579*JBx[5]+0.4472135954999579*JBx[4]+0.9*JBx[3]-0.6708203932499369*JBx[2]-0.6708203932499369*JBx[1]+0.5*JBx[0];
  JBy_con_n[0] = -(0.5999999999999995*JBy[7])-0.5999999999999999*JBy[6]+0.4472135954999579*JBy[5]+0.4472135954999579*JBy[4]+0.9*JBy[3]-0.6708203932499369*JBy[2]-0.6708203932499369*JBy[1]+0.5*JBy[0];
  JBz_con_n[0] = -(0.5999999999999995*JBz[7])-0.5999999999999999*JBz[6]+0.4472135954999579*JBz[5]+0.4472135954999579*JBz[4]+0.9*JBz[3]-0.6708203932499369*JBz[2]-0.6708203932499369*JBz[1]+0.5*JBz[0];
  
  JDx_con_n[1] = 0.75*JDx[7]-0.5590169943749475*JDx[5]+0.4472135954999579*JDx[4]-0.6708203932499369*JDx[1]+0.5*JDx[0];
  JDy_con_n[1] = 0.75*JDy[7]-0.5590169943749475*JDy[5]+0.4472135954999579*JDy[4]-0.6708203932499369*JDy[1]+0.5*JDy[0];
  JDz_con_n[1] = 0.75*JDz[7]-0.5590169943749475*JDz[5]+0.4472135954999579*JDz[4]-0.6708203932499369*JDz[1]+0.5*JDz[0];
  JBx_con_n[1] = 0.75*JBx[7]-0.5590169943749475*JBx[5]+0.4472135954999579*JBx[4]-0.6708203932499369*JBx[1]+0.5*JBx[0];
  JBy_con_n[1] = 0.75*JBy[7]-0.5590169943749475*JBy[5]+0.4472135954999579*JBy[4]-0.6708203932499369*JBy[1]+0.5*JBy[0];
  JBz_con_n[1] = 0.75*JBz[7]-0.5590169943749475*JBz[5]+0.4472135954999579*JBz[4]-0.6708203932499369*JBz[1]+0.5*JBz[0];
  
  JDx_con_n[2] = -(0.5999999999999995*JDx[7])+0.5999999999999999*JDx[6]+0.4472135954999579*JDx[5]+0.4472135954999579*JDx[4]-0.9*JDx[3]+0.6708203932499369*JDx[2]-0.6708203932499369*JDx[1]+0.5*JDx[0];
  JDy_con_n[2] = -(0.5999999999999995*JDy[7])+0.5999999999999999*JDy[6]+0.4472135954999579*JDy[5]+0.4472135954999579*JDy[4]-0.9*JDy[3]+0.6708203932499369*JDy[2]-0.6708203932499369*JDy[1]+0.5*JDy[0];
  JDz_con_n[2] = -(0.5999999999999995*JDz[7])+0.5999999999999999*JDz[6]+0.4472135954999579*JDz[5]+0.4472135954999579*JDz[4]-0.9*JDz[3]+0.6708203932499369*JDz[2]-0.6708203932499369*JDz[1]+0.5*JDz[0];
  JBx_con_n[2] = -(0.5999999999999995*JBx[7])+0.5999999999999999*JBx[6]+0.4472135954999579*JBx[5]+0.4472135954999579*JBx[4]-0.9*JBx[3]+0.6708203932499369*JBx[2]-0.6708203932499369*JBx[1]+0.5*JBx[0];
  JBy_con_n[2] = -(0.5999999999999995*JBy[7])+0.5999999999999999*JBy[6]+0.4472135954999579*JBy[5]+0.4472135954999579*JBy[4]-0.9*JBy[3]+0.6708203932499369*JBy[2]-0.6708203932499369*JBy[1]+0.5*JBy[0];
  JBz_con_n[2] = -(0.5999999999999995*JBz[7])+0.5999999999999999*JBz[6]+0.4472135954999579*JBz[5]+0.4472135954999579*JBz[4]-0.9*JBz[3]+0.6708203932499369*JBz[2]-0.6708203932499369*JBz[1]+0.5*JBz[0];
  
  JDx_con_n[3] = 0.75*JDx[6]+0.4472135954999579*JDx[5]-0.5590169943749475*JDx[4]-0.6708203932499369*JDx[2]+0.5*JDx[0];
  JDy_con_n[3] = 0.75*JDy[6]+0.4472135954999579*JDy[5]-0.5590169943749475*JDy[4]-0.6708203932499369*JDy[2]+0.5*JDy[0];
  JDz_con_n[3] = 0.75*JDz[6]+0.4472135954999579*JDz[5]-0.5590169943749475*JDz[4]-0.6708203932499369*JDz[2]+0.5*JDz[0];
  JBx_con_n[3] = 0.75*JBx[6]+0.4472135954999579*JBx[5]-0.5590169943749475*JBx[4]-0.6708203932499369*JBx[2]+0.5*JBx[0];
  JBy_con_n[3] = 0.75*JBy[6]+0.4472135954999579*JBy[5]-0.5590169943749475*JBy[4]-0.6708203932499369*JBy[2]+0.5*JBy[0];
  JBz_con_n[3] = 0.75*JBz[6]+0.4472135954999579*JBz[5]-0.5590169943749475*JBz[4]-0.6708203932499369*JBz[2]+0.5*JBz[0];
  
  JDx_con_n[4] = -(0.5590169943749475*JDx[5])-0.5590169943749475*JDx[4]+0.5*JDx[0];
  JDy_con_n[4] = -(0.5590169943749475*JDy[5])-0.5590169943749475*JDy[4]+0.5*JDy[0];
  JDz_con_n[4] = -(0.5590169943749475*JDz[5])-0.5590169943749475*JDz[4]+0.5*JDz[0];
  JBx_con_n[4] = -(0.5590169943749475*JBx[5])-0.5590169943749475*JBx[4]+0.5*JBx[0];
  JBy_con_n[4] = -(0.5590169943749475*JBy[5])-0.5590169943749475*JBy[4]+0.5*JBy[0];
  JBz_con_n[4] = -(0.5590169943749475*JBz[5])-0.5590169943749475*JBz[4]+0.5*JBz[0];
  
  JDx_con_n[5] = -(0.75*JDx[6])+0.4472135954999579*JDx[5]-0.5590169943749475*JDx[4]+0.6708203932499369*JDx[2]+0.5*JDx[0];
  JDy_con_n[5] = -(0.75*JDy[6])+0.4472135954999579*JDy[5]-0.5590169943749475*JDy[4]+0.6708203932499369*JDy[2]+0.5*JDy[0];
  JDz_con_n[5] = -(0.75*JDz[6])+0.4472135954999579*JDz[5]-0.5590169943749475*JDz[4]+0.6708203932499369*JDz[2]+0.5*JDz[0];
  JBx_con_n[5] = -(0.75*JBx[6])+0.4472135954999579*JBx[5]-0.5590169943749475*JBx[4]+0.6708203932499369*JBx[2]+0.5*JBx[0];
  JBy_con_n[5] = -(0.75*JBy[6])+0.4472135954999579*JBy[5]-0.5590169943749475*JBy[4]+0.6708203932499369*JBy[2]+0.5*JBy[0];
  JBz_con_n[5] = -(0.75*JBz[6])+0.4472135954999579*JBz[5]-0.5590169943749475*JBz[4]+0.6708203932499369*JBz[2]+0.5*JBz[0];
  
  JDx_con_n[6] = 0.5999999999999995*JDx[7]-0.5999999999999999*JDx[6]+0.4472135954999579*JDx[5]+0.4472135954999579*JDx[4]-0.9*JDx[3]-0.6708203932499369*JDx[2]+0.6708203932499369*JDx[1]+0.5*JDx[0];
  JDy_con_n[6] = 0.5999999999999995*JDy[7]-0.5999999999999999*JDy[6]+0.4472135954999579*JDy[5]+0.4472135954999579*JDy[4]-0.9*JDy[3]-0.6708203932499369*JDy[2]+0.6708203932499369*JDy[1]+0.5*JDy[0];
  JDz_con_n[6] = 0.5999999999999995*JDz[7]-0.5999999999999999*JDz[6]+0.4472135954999579*JDz[5]+0.4472135954999579*JDz[4]-0.9*JDz[3]-0.6708203932499369*JDz[2]+0.6708203932499369*JDz[1]+0.5*JDz[0];
  JBx_con_n[6] = 0.5999999999999995*JBx[7]-0.5999999999999999*JBx[6]+0.4472135954999579*JBx[5]+0.4472135954999579*JBx[4]-0.9*JBx[3]-0.6708203932499369*JBx[2]+0.6708203932499369*JBx[1]+0.5*JBx[0];
  JBy_con_n[6] = 0.5999999999999995*JBy[7]-0.5999999999999999*JBy[6]+0.4472135954999579*JBy[5]+0.4472135954999579*JBy[4]-0.9*JBy[3]-0.6708203932499369*JBy[2]+0.6708203932499369*JBy[1]+0.5*JBy[0];
  JBz_con_n[6] = 0.5999999999999995*JBz[7]-0.5999999999999999*JBz[6]+0.4472135954999579*JBz[5]+0.4472135954999579*JBz[4]-0.9*JBz[3]-0.6708203932499369*JBz[2]+0.6708203932499369*JBz[1]+0.5*JBz[0];
  
  JDx_con_n[7] = -(0.75*JDx[7])-0.5590169943749475*JDx[5]+0.4472135954999579*JDx[4]+0.6708203932499369*JDx[1]+0.5*JDx[0];
  JDy_con_n[7] = -(0.75*JDy[7])-0.5590169943749475*JDy[5]+0.4472135954999579*JDy[4]+0.6708203932499369*JDy[1]+0.5*JDy[0];
  JDz_con_n[7] = -(0.75*JDz[7])-0.5590169943749475*JDz[5]+0.4472135954999579*JDz[4]+0.6708203932499369*JDz[1]+0.5*JDz[0];
  JBx_con_n[7] = -(0.75*JBx[7])-0.5590169943749475*JBx[5]+0.4472135954999579*JBx[4]+0.6708203932499369*JBx[1]+0.5*JBx[0];
  JBy_con_n[7] = -(0.75*JBy[7])-0.5590169943749475*JBy[5]+0.4472135954999579*JBy[4]+0.6708203932499369*JBy[1]+0.5*JBy[0];
  JBz_con_n[7] = -(0.75*JBz[7])-0.5590169943749475*JBz[5]+0.4472135954999579*JBz[4]+0.6708203932499369*JBz[1]+0.5*JBz[0];
  
  JDx_con_n[8] = 0.5999999999999995*JDx[7]+0.5999999999999999*JDx[6]+0.4472135954999579*JDx[5]+0.4472135954999579*JDx[4]+0.9*JDx[3]+0.6708203932499369*JDx[2]+0.6708203932499369*JDx[1]+0.5*JDx[0];
  JDy_con_n[8] = 0.5999999999999995*JDy[7]+0.5999999999999999*JDy[6]+0.4472135954999579*JDy[5]+0.4472135954999579*JDy[4]+0.9*JDy[3]+0.6708203932499369*JDy[2]+0.6708203932499369*JDy[1]+0.5*JDy[0];
  JDz_con_n[8] = 0.5999999999999995*JDz[7]+0.5999999999999999*JDz[6]+0.4472135954999579*JDz[5]+0.4472135954999579*JDz[4]+0.9*JDz[3]+0.6708203932499369*JDz[2]+0.6708203932499369*JDz[1]+0.5*JDz[0];
  JBx_con_n[8] = 0.5999999999999995*JBx[7]+0.5999999999999999*JBx[6]+0.4472135954999579*JBx[5]+0.4472135954999579*JBx[4]+0.9*JBx[3]+0.6708203932499369*JBx[2]+0.6708203932499369*JBx[1]+0.5*JBx[0];
  JBy_con_n[8] = 0.5999999999999995*JBy[7]+0.5999999999999999*JBy[6]+0.4472135954999579*JBy[5]+0.4472135954999579*JBy[4]+0.9*JBy[3]+0.6708203932499369*JBy[2]+0.6708203932499369*JBy[1]+0.5*JBy[0];
  JBz_con_n[8] = 0.5999999999999995*JBz[7]+0.5999999999999999*JBz[6]+0.4472135954999579*JBz[5]+0.4472135954999579*JBz[4]+0.9*JBz[3]+0.6708203932499369*JBz[2]+0.6708203932499369*JBz[1]+0.5*JBz[0];
  
  const double *h_xx_nodal = &h_ij_nodal[0]; 
  const double *h_xy_nodal = &h_ij_nodal[9]; 
  const double *h_xz_nodal = &h_ij_nodal[18]; 
  const double *h_yy_nodal = &h_ij_nodal[27]; 
  const double *h_yz_nodal = &h_ij_nodal[36]; 
  const double *h_zz_nodal = &h_ij_nodal[45]; 
  const double *h_xx_inv_nodal = &h_ij_inv_nodal[0]; 
  const double *h_xy_inv_nodal = &h_ij_inv_nodal[9]; 
  const double *h_xz_inv_nodal = &h_ij_inv_nodal[18]; 
  const double *h_yy_inv_nodal = &h_ij_inv_nodal[27]; 
  const double *h_yz_inv_nodal = &h_ij_inv_nodal[36]; 
  const double *h_zz_inv_nodal = &h_ij_inv_nodal[45]; 
  const double *vierb_xx_cov_nodal = &vierb_cov_nodal[0]; 
  const double *vierb_xy_cov_nodal = &vierb_cov_nodal[9]; 
  const double vierb_xz_cov_nodal[9] = {0.0}; 
  const double *vierb_yx_cov_nodal = &vierb_cov_nodal[18]; 
  const double *vierb_yy_cov_nodal = &vierb_cov_nodal[27]; 
  const double vierb_yz_cov_nodal[9] = {0.0}; 
  const double vierb_zx_cov_nodal[9] = {0.0}; 
  const double vierb_zy_cov_nodal[9] = {0.0}; 
  double vierb_zz_cov_nodal[9] = {0.0}; 
  const double *vierb_xx_con_nodal = &vierb_con_nodal[0]; 
  const double *vierb_xy_con_nodal = &vierb_con_nodal[9]; 
  const double vierb_xz_con_nodal[9] = {0.0}; 
  const double *vierb_yx_con_nodal = &vierb_con_nodal[18]; 
  const double *vierb_yy_con_nodal = &vierb_con_nodal[27]; 
  const double vierb_yz_con_nodal[9] = {0.0}; 
  const double vierb_zx_con_nodal[9] = {0.0}; 
  const double vierb_zy_con_nodal[9] = {0.0}; 
  double vierb_zz_con_nodal[9] = {0.0}; 
  for (int i=0; i<9; ++i) {
    vierb_zz_cov_nodal[i] = 1.0;
    vierb_zz_con_nodal[i] = 1.0;
  }
  
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[9]; 
  const double *shift_nodal_z = &shift_nodal[18]; 
  
  for (int i=0; i<9; ++i) {
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
  
  for (int i=0; i<9; ++i) {
    E_conf_for_force_n[i + 0*9] = vierb_xx_con_nodal[i]*Ex_n[i] + vierb_yx_con_nodal[i]*Ey_n[i]  + vierb_zx_con_nodal[i]*Ez_n[i];
    E_conf_for_force_n[i + 1*9] = vierb_xy_con_nodal[i]*Ex_n[i] + vierb_yy_con_nodal[i]*Ey_n[i]  + vierb_zy_con_nodal[i]*Ez_n[i];
    E_conf_for_force_n[i + 2*9] = vierb_xz_con_nodal[i]*Ex_n[i] + vierb_yz_con_nodal[i]*Ey_n[i]  + vierb_zz_con_nodal[i]*Ez_n[i];
    B_conf_for_force_n[i + 0*9] = vierb_xx_cov_nodal[i]*Bx_n[i] + vierb_yx_cov_nodal[i]*By_n[i] + vierb_zx_cov_nodal[i]*Bz_n[i];
    B_conf_for_force_n[i + 1*9] = vierb_xy_cov_nodal[i]*Bx_n[i] + vierb_yy_cov_nodal[i]*By_n[i] + vierb_zy_cov_nodal[i]*Bz_n[i];
    B_conf_for_force_n[i + 2*9] = vierb_xz_cov_nodal[i]*Bx_n[i] + vierb_yz_cov_nodal[i]*By_n[i] + vierb_zz_cov_nodal[i]*Bz_n[i];
  }
  
  double *E_conf_for_force_comp;
  double *B_conf_for_force_comp;
  E_conf_for_force_comp = &E_conf_for_force[0];
  E_conf_for_force_comp[0] = 0.15432098765432098*E_conf_for_force_n[8]+0.24691358024691357*E_conf_for_force_n[7]+0.15432098765432098*E_conf_for_force_n[6]+0.24691358024691357*E_conf_for_force_n[5]+0.3950617283950617*E_conf_for_force_n[4]+0.24691358024691357*E_conf_for_force_n[3]+0.15432098765432098*E_conf_for_force_n[2]+0.24691358024691357*E_conf_for_force_n[1]+0.15432098765432098*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[1] = 0.20704333124998056*E_conf_for_force_n[8]+0.33126932999996883*E_conf_for_force_n[7]+0.20704333124998056*E_conf_for_force_n[6]-0.20704333124998056*E_conf_for_force_n[2]-0.33126932999996883*E_conf_for_force_n[1]-0.20704333124998056*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[2] = 0.20704333124998056*E_conf_for_force_n[8]-0.20704333124998056*E_conf_for_force_n[6]+0.33126932999996883*E_conf_for_force_n[5]-0.33126932999996883*E_conf_for_force_n[3]+0.20704333124998056*E_conf_for_force_n[2]-0.20704333124998056*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[3] = 0.2777777777777778*E_conf_for_force_n[8]-0.2777777777777778*E_conf_for_force_n[6]-0.2777777777777778*E_conf_for_force_n[2]+0.2777777777777778*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[4] = 0.13802888749998704*E_conf_for_force_n[8]+0.22084621999997922*E_conf_for_force_n[7]+0.13802888749998704*E_conf_for_force_n[6]-0.2760577749999741*E_conf_for_force_n[5]-0.44169243999995844*E_conf_for_force_n[4]-0.2760577749999741*E_conf_for_force_n[3]+0.13802888749998704*E_conf_for_force_n[2]+0.22084621999997922*E_conf_for_force_n[1]+0.13802888749998704*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[5] = 0.13802888749998704*E_conf_for_force_n[8]-0.2760577749999741*E_conf_for_force_n[7]+0.13802888749998704*E_conf_for_force_n[6]+0.22084621999997922*E_conf_for_force_n[5]-0.44169243999995844*E_conf_for_force_n[4]+0.22084621999997922*E_conf_for_force_n[3]+0.13802888749998704*E_conf_for_force_n[2]-0.2760577749999741*E_conf_for_force_n[1]+0.13802888749998704*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[6] = 0.18518518518518526*E_conf_for_force_n[8]-0.18518518518518526*E_conf_for_force_n[6]-0.3703703703703705*E_conf_for_force_n[5]+0.3703703703703705*E_conf_for_force_n[3]+0.18518518518518526*E_conf_for_force_n[2]-0.18518518518518526*E_conf_for_force_n[0]; 
  E_conf_for_force_comp[7] = 0.18518518518518526*E_conf_for_force_n[8]-0.3703703703703705*E_conf_for_force_n[7]+0.18518518518518526*E_conf_for_force_n[6]-0.18518518518518526*E_conf_for_force_n[2]+0.3703703703703705*E_conf_for_force_n[1]-0.18518518518518526*E_conf_for_force_n[0]; 

  E_conf_for_force_comp = &E_conf_for_force[8];
  E_conf_for_force_comp[0] = 0.15432098765432098*E_conf_for_force_n[17]+0.24691358024691357*E_conf_for_force_n[16]+0.15432098765432098*E_conf_for_force_n[15]+0.24691358024691357*E_conf_for_force_n[14]+0.3950617283950617*E_conf_for_force_n[13]+0.24691358024691357*E_conf_for_force_n[12]+0.15432098765432098*E_conf_for_force_n[11]+0.24691358024691357*E_conf_for_force_n[10]+0.15432098765432098*E_conf_for_force_n[9]; 
  E_conf_for_force_comp[1] = 0.20704333124998056*E_conf_for_force_n[17]+0.33126932999996883*E_conf_for_force_n[16]+0.20704333124998056*E_conf_for_force_n[15]-0.20704333124998056*E_conf_for_force_n[11]-0.33126932999996883*E_conf_for_force_n[10]-0.20704333124998056*E_conf_for_force_n[9]; 
  E_conf_for_force_comp[2] = 0.20704333124998056*E_conf_for_force_n[17]-0.20704333124998056*E_conf_for_force_n[15]+0.33126932999996883*E_conf_for_force_n[14]-0.33126932999996883*E_conf_for_force_n[12]+0.20704333124998056*E_conf_for_force_n[11]-0.20704333124998056*E_conf_for_force_n[9]; 
  E_conf_for_force_comp[3] = 0.2777777777777778*E_conf_for_force_n[17]-0.2777777777777778*E_conf_for_force_n[15]-0.2777777777777778*E_conf_for_force_n[11]+0.2777777777777778*E_conf_for_force_n[9]; 
  E_conf_for_force_comp[4] = 0.13802888749998704*E_conf_for_force_n[17]+0.22084621999997922*E_conf_for_force_n[16]+0.13802888749998704*E_conf_for_force_n[15]-0.2760577749999741*E_conf_for_force_n[14]-0.44169243999995844*E_conf_for_force_n[13]-0.2760577749999741*E_conf_for_force_n[12]+0.13802888749998704*E_conf_for_force_n[11]+0.22084621999997922*E_conf_for_force_n[10]+0.13802888749998704*E_conf_for_force_n[9]; 
  E_conf_for_force_comp[5] = 0.13802888749998704*E_conf_for_force_n[17]-0.2760577749999741*E_conf_for_force_n[16]+0.13802888749998704*E_conf_for_force_n[15]+0.22084621999997922*E_conf_for_force_n[14]-0.44169243999995844*E_conf_for_force_n[13]+0.22084621999997922*E_conf_for_force_n[12]+0.13802888749998704*E_conf_for_force_n[11]-0.2760577749999741*E_conf_for_force_n[10]+0.13802888749998704*E_conf_for_force_n[9]; 
  E_conf_for_force_comp[6] = 0.18518518518518526*E_conf_for_force_n[17]-0.18518518518518526*E_conf_for_force_n[15]-0.3703703703703705*E_conf_for_force_n[14]+0.3703703703703705*E_conf_for_force_n[12]+0.18518518518518526*E_conf_for_force_n[11]-0.18518518518518526*E_conf_for_force_n[9]; 
  E_conf_for_force_comp[7] = 0.18518518518518526*E_conf_for_force_n[17]-0.3703703703703705*E_conf_for_force_n[16]+0.18518518518518526*E_conf_for_force_n[15]-0.18518518518518526*E_conf_for_force_n[11]+0.3703703703703705*E_conf_for_force_n[10]-0.18518518518518526*E_conf_for_force_n[9]; 

  E_conf_for_force_comp = &E_conf_for_force[16];
  E_conf_for_force_comp[0] = 0.15432098765432098*E_conf_for_force_n[26]+0.24691358024691357*E_conf_for_force_n[25]+0.15432098765432098*E_conf_for_force_n[24]+0.24691358024691357*E_conf_for_force_n[23]+0.3950617283950617*E_conf_for_force_n[22]+0.24691358024691357*E_conf_for_force_n[21]+0.15432098765432098*E_conf_for_force_n[20]+0.24691358024691357*E_conf_for_force_n[19]+0.15432098765432098*E_conf_for_force_n[18]; 
  E_conf_for_force_comp[1] = 0.20704333124998056*E_conf_for_force_n[26]+0.33126932999996883*E_conf_for_force_n[25]+0.20704333124998056*E_conf_for_force_n[24]-0.20704333124998056*E_conf_for_force_n[20]-0.33126932999996883*E_conf_for_force_n[19]-0.20704333124998056*E_conf_for_force_n[18]; 
  E_conf_for_force_comp[2] = 0.20704333124998056*E_conf_for_force_n[26]-0.20704333124998056*E_conf_for_force_n[24]+0.33126932999996883*E_conf_for_force_n[23]-0.33126932999996883*E_conf_for_force_n[21]+0.20704333124998056*E_conf_for_force_n[20]-0.20704333124998056*E_conf_for_force_n[18]; 
  E_conf_for_force_comp[3] = 0.2777777777777778*E_conf_for_force_n[26]-0.2777777777777778*E_conf_for_force_n[24]-0.2777777777777778*E_conf_for_force_n[20]+0.2777777777777778*E_conf_for_force_n[18]; 
  E_conf_for_force_comp[4] = 0.13802888749998704*E_conf_for_force_n[26]+0.22084621999997922*E_conf_for_force_n[25]+0.13802888749998704*E_conf_for_force_n[24]-0.2760577749999741*E_conf_for_force_n[23]-0.44169243999995844*E_conf_for_force_n[22]-0.2760577749999741*E_conf_for_force_n[21]+0.13802888749998704*E_conf_for_force_n[20]+0.22084621999997922*E_conf_for_force_n[19]+0.13802888749998704*E_conf_for_force_n[18]; 
  E_conf_for_force_comp[5] = 0.13802888749998704*E_conf_for_force_n[26]-0.2760577749999741*E_conf_for_force_n[25]+0.13802888749998704*E_conf_for_force_n[24]+0.22084621999997922*E_conf_for_force_n[23]-0.44169243999995844*E_conf_for_force_n[22]+0.22084621999997922*E_conf_for_force_n[21]+0.13802888749998704*E_conf_for_force_n[20]-0.2760577749999741*E_conf_for_force_n[19]+0.13802888749998704*E_conf_for_force_n[18]; 
  E_conf_for_force_comp[6] = 0.18518518518518526*E_conf_for_force_n[26]-0.18518518518518526*E_conf_for_force_n[24]-0.3703703703703705*E_conf_for_force_n[23]+0.3703703703703705*E_conf_for_force_n[21]+0.18518518518518526*E_conf_for_force_n[20]-0.18518518518518526*E_conf_for_force_n[18]; 
  E_conf_for_force_comp[7] = 0.18518518518518526*E_conf_for_force_n[26]-0.3703703703703705*E_conf_for_force_n[25]+0.18518518518518526*E_conf_for_force_n[24]-0.18518518518518526*E_conf_for_force_n[20]+0.3703703703703705*E_conf_for_force_n[19]-0.18518518518518526*E_conf_for_force_n[18]; 

  B_conf_for_force_comp = &B_conf_for_force[0];
  B_conf_for_force_comp[0] = 0.15432098765432098*B_conf_for_force_n[8]+0.24691358024691357*B_conf_for_force_n[7]+0.15432098765432098*B_conf_for_force_n[6]+0.24691358024691357*B_conf_for_force_n[5]+0.3950617283950617*B_conf_for_force_n[4]+0.24691358024691357*B_conf_for_force_n[3]+0.15432098765432098*B_conf_for_force_n[2]+0.24691358024691357*B_conf_for_force_n[1]+0.15432098765432098*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[1] = 0.20704333124998056*B_conf_for_force_n[8]+0.33126932999996883*B_conf_for_force_n[7]+0.20704333124998056*B_conf_for_force_n[6]-0.20704333124998056*B_conf_for_force_n[2]-0.33126932999996883*B_conf_for_force_n[1]-0.20704333124998056*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[2] = 0.20704333124998056*B_conf_for_force_n[8]-0.20704333124998056*B_conf_for_force_n[6]+0.33126932999996883*B_conf_for_force_n[5]-0.33126932999996883*B_conf_for_force_n[3]+0.20704333124998056*B_conf_for_force_n[2]-0.20704333124998056*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[3] = 0.2777777777777778*B_conf_for_force_n[8]-0.2777777777777778*B_conf_for_force_n[6]-0.2777777777777778*B_conf_for_force_n[2]+0.2777777777777778*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[4] = 0.13802888749998704*B_conf_for_force_n[8]+0.22084621999997922*B_conf_for_force_n[7]+0.13802888749998704*B_conf_for_force_n[6]-0.2760577749999741*B_conf_for_force_n[5]-0.44169243999995844*B_conf_for_force_n[4]-0.2760577749999741*B_conf_for_force_n[3]+0.13802888749998704*B_conf_for_force_n[2]+0.22084621999997922*B_conf_for_force_n[1]+0.13802888749998704*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[5] = 0.13802888749998704*B_conf_for_force_n[8]-0.2760577749999741*B_conf_for_force_n[7]+0.13802888749998704*B_conf_for_force_n[6]+0.22084621999997922*B_conf_for_force_n[5]-0.44169243999995844*B_conf_for_force_n[4]+0.22084621999997922*B_conf_for_force_n[3]+0.13802888749998704*B_conf_for_force_n[2]-0.2760577749999741*B_conf_for_force_n[1]+0.13802888749998704*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[6] = 0.18518518518518526*B_conf_for_force_n[8]-0.18518518518518526*B_conf_for_force_n[6]-0.3703703703703705*B_conf_for_force_n[5]+0.3703703703703705*B_conf_for_force_n[3]+0.18518518518518526*B_conf_for_force_n[2]-0.18518518518518526*B_conf_for_force_n[0]; 
  B_conf_for_force_comp[7] = 0.18518518518518526*B_conf_for_force_n[8]-0.3703703703703705*B_conf_for_force_n[7]+0.18518518518518526*B_conf_for_force_n[6]-0.18518518518518526*B_conf_for_force_n[2]+0.3703703703703705*B_conf_for_force_n[1]-0.18518518518518526*B_conf_for_force_n[0]; 

  B_conf_for_force_comp = &B_conf_for_force[8];
  B_conf_for_force_comp[0] = 0.15432098765432098*B_conf_for_force_n[17]+0.24691358024691357*B_conf_for_force_n[16]+0.15432098765432098*B_conf_for_force_n[15]+0.24691358024691357*B_conf_for_force_n[14]+0.3950617283950617*B_conf_for_force_n[13]+0.24691358024691357*B_conf_for_force_n[12]+0.15432098765432098*B_conf_for_force_n[11]+0.24691358024691357*B_conf_for_force_n[10]+0.15432098765432098*B_conf_for_force_n[9]; 
  B_conf_for_force_comp[1] = 0.20704333124998056*B_conf_for_force_n[17]+0.33126932999996883*B_conf_for_force_n[16]+0.20704333124998056*B_conf_for_force_n[15]-0.20704333124998056*B_conf_for_force_n[11]-0.33126932999996883*B_conf_for_force_n[10]-0.20704333124998056*B_conf_for_force_n[9]; 
  B_conf_for_force_comp[2] = 0.20704333124998056*B_conf_for_force_n[17]-0.20704333124998056*B_conf_for_force_n[15]+0.33126932999996883*B_conf_for_force_n[14]-0.33126932999996883*B_conf_for_force_n[12]+0.20704333124998056*B_conf_for_force_n[11]-0.20704333124998056*B_conf_for_force_n[9]; 
  B_conf_for_force_comp[3] = 0.2777777777777778*B_conf_for_force_n[17]-0.2777777777777778*B_conf_for_force_n[15]-0.2777777777777778*B_conf_for_force_n[11]+0.2777777777777778*B_conf_for_force_n[9]; 
  B_conf_for_force_comp[4] = 0.13802888749998704*B_conf_for_force_n[17]+0.22084621999997922*B_conf_for_force_n[16]+0.13802888749998704*B_conf_for_force_n[15]-0.2760577749999741*B_conf_for_force_n[14]-0.44169243999995844*B_conf_for_force_n[13]-0.2760577749999741*B_conf_for_force_n[12]+0.13802888749998704*B_conf_for_force_n[11]+0.22084621999997922*B_conf_for_force_n[10]+0.13802888749998704*B_conf_for_force_n[9]; 
  B_conf_for_force_comp[5] = 0.13802888749998704*B_conf_for_force_n[17]-0.2760577749999741*B_conf_for_force_n[16]+0.13802888749998704*B_conf_for_force_n[15]+0.22084621999997922*B_conf_for_force_n[14]-0.44169243999995844*B_conf_for_force_n[13]+0.22084621999997922*B_conf_for_force_n[12]+0.13802888749998704*B_conf_for_force_n[11]-0.2760577749999741*B_conf_for_force_n[10]+0.13802888749998704*B_conf_for_force_n[9]; 
  B_conf_for_force_comp[6] = 0.18518518518518526*B_conf_for_force_n[17]-0.18518518518518526*B_conf_for_force_n[15]-0.3703703703703705*B_conf_for_force_n[14]+0.3703703703703705*B_conf_for_force_n[12]+0.18518518518518526*B_conf_for_force_n[11]-0.18518518518518526*B_conf_for_force_n[9]; 
  B_conf_for_force_comp[7] = 0.18518518518518526*B_conf_for_force_n[17]-0.3703703703703705*B_conf_for_force_n[16]+0.18518518518518526*B_conf_for_force_n[15]-0.18518518518518526*B_conf_for_force_n[11]+0.3703703703703705*B_conf_for_force_n[10]-0.18518518518518526*B_conf_for_force_n[9]; 

  B_conf_for_force_comp = &B_conf_for_force[16];
  B_conf_for_force_comp[0] = 0.15432098765432098*B_conf_for_force_n[26]+0.24691358024691357*B_conf_for_force_n[25]+0.15432098765432098*B_conf_for_force_n[24]+0.24691358024691357*B_conf_for_force_n[23]+0.3950617283950617*B_conf_for_force_n[22]+0.24691358024691357*B_conf_for_force_n[21]+0.15432098765432098*B_conf_for_force_n[20]+0.24691358024691357*B_conf_for_force_n[19]+0.15432098765432098*B_conf_for_force_n[18]; 
  B_conf_for_force_comp[1] = 0.20704333124998056*B_conf_for_force_n[26]+0.33126932999996883*B_conf_for_force_n[25]+0.20704333124998056*B_conf_for_force_n[24]-0.20704333124998056*B_conf_for_force_n[20]-0.33126932999996883*B_conf_for_force_n[19]-0.20704333124998056*B_conf_for_force_n[18]; 
  B_conf_for_force_comp[2] = 0.20704333124998056*B_conf_for_force_n[26]-0.20704333124998056*B_conf_for_force_n[24]+0.33126932999996883*B_conf_for_force_n[23]-0.33126932999996883*B_conf_for_force_n[21]+0.20704333124998056*B_conf_for_force_n[20]-0.20704333124998056*B_conf_for_force_n[18]; 
  B_conf_for_force_comp[3] = 0.2777777777777778*B_conf_for_force_n[26]-0.2777777777777778*B_conf_for_force_n[24]-0.2777777777777778*B_conf_for_force_n[20]+0.2777777777777778*B_conf_for_force_n[18]; 
  B_conf_for_force_comp[4] = 0.13802888749998704*B_conf_for_force_n[26]+0.22084621999997922*B_conf_for_force_n[25]+0.13802888749998704*B_conf_for_force_n[24]-0.2760577749999741*B_conf_for_force_n[23]-0.44169243999995844*B_conf_for_force_n[22]-0.2760577749999741*B_conf_for_force_n[21]+0.13802888749998704*B_conf_for_force_n[20]+0.22084621999997922*B_conf_for_force_n[19]+0.13802888749998704*B_conf_for_force_n[18]; 
  B_conf_for_force_comp[5] = 0.13802888749998704*B_conf_for_force_n[26]-0.2760577749999741*B_conf_for_force_n[25]+0.13802888749998704*B_conf_for_force_n[24]+0.22084621999997922*B_conf_for_force_n[23]-0.44169243999995844*B_conf_for_force_n[22]+0.22084621999997922*B_conf_for_force_n[21]+0.13802888749998704*B_conf_for_force_n[20]-0.2760577749999741*B_conf_for_force_n[19]+0.13802888749998704*B_conf_for_force_n[18]; 
  B_conf_for_force_comp[6] = 0.18518518518518526*B_conf_for_force_n[26]-0.18518518518518526*B_conf_for_force_n[24]-0.3703703703703705*B_conf_for_force_n[23]+0.3703703703703705*B_conf_for_force_n[21]+0.18518518518518526*B_conf_for_force_n[20]-0.18518518518518526*B_conf_for_force_n[18]; 
  B_conf_for_force_comp[7] = 0.18518518518518526*B_conf_for_force_n[26]-0.3703703703703705*B_conf_for_force_n[25]+0.18518518518518526*B_conf_for_force_n[24]-0.18518518518518526*B_conf_for_force_n[20]+0.3703703703703705*B_conf_for_force_n[19]-0.18518518518518526*B_conf_for_force_n[18]; 

} 
