#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void gr_maxwell_current_dep_2x_3v_ser_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
  const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out) 
{ 
  
  double Jc_rho_n[4] = {0.0};
  double Jc_Jx_flat_n[4] = {0.0};
  double Jc_Jy_flat_n[4] = {0.0};
  double Jc_Jz_flat_n[4] = {0.0};
  
  const double *Jc_rho = &m0[0]; 
  const double *Jc_Jx_flat = &m1i[0]; 
  const double *Jc_Jy_flat = &m1i[4]; 
  const double *Jc_Jz_flat = &m1i[8]; 
  
  double Current_n[12] = {0.0};
  
  Jc_rho_n[0] = 0.5*Jc_rho[3]-0.5*Jc_rho[2]-0.5*Jc_rho[1]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[0] = 0.5*Jc_Jx_flat[3]-0.5*Jc_Jx_flat[2]-0.5*Jc_Jx_flat[1]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[0] = 0.5*Jc_Jy_flat[3]-0.5*Jc_Jy_flat[2]-0.5*Jc_Jy_flat[1]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[0] = 0.5*Jc_Jz_flat[3]-0.5*Jc_Jz_flat[2]-0.5*Jc_Jz_flat[1]+0.5*Jc_Jz_flat[0];
  
  Jc_rho_n[1] = -(0.5*Jc_rho[3])+0.5*Jc_rho[2]-0.5*Jc_rho[1]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[1] = -(0.5*Jc_Jx_flat[3])+0.5*Jc_Jx_flat[2]-0.5*Jc_Jx_flat[1]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[1] = -(0.5*Jc_Jy_flat[3])+0.5*Jc_Jy_flat[2]-0.5*Jc_Jy_flat[1]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[1] = -(0.5*Jc_Jz_flat[3])+0.5*Jc_Jz_flat[2]-0.5*Jc_Jz_flat[1]+0.5*Jc_Jz_flat[0];
  
  Jc_rho_n[2] = -(0.5*Jc_rho[3])-0.5*Jc_rho[2]+0.5*Jc_rho[1]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[2] = -(0.5*Jc_Jx_flat[3])-0.5*Jc_Jx_flat[2]+0.5*Jc_Jx_flat[1]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[2] = -(0.5*Jc_Jy_flat[3])-0.5*Jc_Jy_flat[2]+0.5*Jc_Jy_flat[1]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[2] = -(0.5*Jc_Jz_flat[3])-0.5*Jc_Jz_flat[2]+0.5*Jc_Jz_flat[1]+0.5*Jc_Jz_flat[0];
  
  Jc_rho_n[3] = 0.5*Jc_rho[3]+0.5*Jc_rho[2]+0.5*Jc_rho[1]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[3] = 0.5*Jc_Jx_flat[3]+0.5*Jc_Jx_flat[2]+0.5*Jc_Jx_flat[1]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[3] = 0.5*Jc_Jy_flat[3]+0.5*Jc_Jy_flat[2]+0.5*Jc_Jy_flat[1]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[3] = 0.5*Jc_Jz_flat[3]+0.5*Jc_Jz_flat[2]+0.5*Jc_Jz_flat[1]+0.5*Jc_Jz_flat[0];
  
  const double *vierb_xx_con_nodal = &vierb_con_nodal[0]; 
  const double *vierb_xy_con_nodal = &vierb_con_nodal[4]; 
  const double *vierb_xz_con_nodal = &vierb_con_nodal[8]; 
  const double *vierb_yx_con_nodal = &vierb_con_nodal[12]; 
  const double *vierb_yy_con_nodal = &vierb_con_nodal[16]; 
  const double *vierb_yz_con_nodal = &vierb_con_nodal[20]; 
  const double *vierb_zx_con_nodal = &vierb_con_nodal[24]; 
  const double *vierb_zy_con_nodal = &vierb_con_nodal[28]; 
  const double *vierb_zz_con_nodal = &vierb_con_nodal[32]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[4]; 
  const double *shift_nodal_z = &shift_nodal[8]; 
  
  for (int i=0; i<4; ++i) {
    Current_n[i + 0*4] = q_over_eps0 * (Jc_rho_n[i] * shift_nodal_x[i] - lapse_nodal[i] * ( vierb_xx_con_nodal[i] * Jc_Jx_flat_n[i] + vierb_xy_con_nodal[i] * Jc_Jy_flat_n[i] + vierb_xz_con_nodal[i] * Jc_Jz_flat_n[i] )); 
    Current_n[i + 1*4] = q_over_eps0 * (Jc_rho_n[i] * shift_nodal_y[i] - lapse_nodal[i] * ( vierb_yx_con_nodal[i] * Jc_Jx_flat_n[i] + vierb_yy_con_nodal[i] * Jc_Jy_flat_n[i] + vierb_yz_con_nodal[i] * Jc_Jz_flat_n[i] )); 
    Current_n[i + 2*4] = q_over_eps0 * (Jc_rho_n[i] * shift_nodal_z[i] - lapse_nodal[i] * ( vierb_zx_con_nodal[i] * Jc_Jx_flat_n[i] + vierb_zy_con_nodal[i] * Jc_Jy_flat_n[i] + vierb_zz_con_nodal[i] * Jc_Jz_flat_n[i] )); 
  }
  
  double *Current_comp;
  Current_comp = &out[0];
  Current_comp[0] += 0.5*Current_n[3]+0.5*Current_n[2]+0.5*Current_n[1]+0.5*Current_n[0]; 
  Current_comp[1] += 0.5*Current_n[3]+0.5*Current_n[2]-0.5*Current_n[1]-0.5*Current_n[0]; 
  Current_comp[2] += 0.5*Current_n[3]-0.5*Current_n[2]+0.5*Current_n[1]-0.5*Current_n[0]; 
  Current_comp[3] += 0.5*Current_n[3]-0.5*Current_n[2]-0.5*Current_n[1]+0.5*Current_n[0]; 

  Current_comp = &out[4];
  Current_comp[0] += 0.5*Current_n[7]+0.5*Current_n[6]+0.5*Current_n[5]+0.5*Current_n[4]; 
  Current_comp[1] += 0.5*Current_n[7]+0.5*Current_n[6]-0.5*Current_n[5]-0.5*Current_n[4]; 
  Current_comp[2] += 0.5*Current_n[7]-0.5*Current_n[6]+0.5*Current_n[5]-0.5*Current_n[4]; 
  Current_comp[3] += 0.5*Current_n[7]-0.5*Current_n[6]-0.5*Current_n[5]+0.5*Current_n[4]; 

  Current_comp = &out[8];
  Current_comp[0] += 0.5*Current_n[11]+0.5*Current_n[10]+0.5*Current_n[9]+0.5*Current_n[8]; 
  Current_comp[1] += 0.5*Current_n[11]+0.5*Current_n[10]-0.5*Current_n[9]-0.5*Current_n[8]; 
  Current_comp[2] += 0.5*Current_n[11]-0.5*Current_n[10]+0.5*Current_n[9]-0.5*Current_n[8]; 
  Current_comp[3] += 0.5*Current_n[11]-0.5*Current_n[10]-0.5*Current_n[9]+0.5*Current_n[8]; 

}
