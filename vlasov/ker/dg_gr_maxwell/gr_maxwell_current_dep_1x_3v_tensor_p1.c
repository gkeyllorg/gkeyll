#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void gr_maxwell_current_dep_1x_3v_tensor_p1(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
  const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out) 
{ 
  
  double Jc_rho_n[2] = {0.0};
  double Jc_Jx_flat_n[2] = {0.0};
  double Jc_Jy_flat_n[2] = {0.0};
  double Jc_Jz_flat_n[2] = {0.0};
  
  const double *Jc_rho = &m0[0]; 
  const double *Jc_Jx_flat = &m1i[0]; 
  const double *Jc_Jy_flat = &m1i[2]; 
  const double *Jc_Jz_flat = &m1i[4]; 
  
  double Current_n[6] = {0.0};
  
  Jc_rho_n[0] = 0.7071067811865475*Jc_rho[0]-0.7071067811865475*Jc_rho[1];
  Jc_Jx_flat_n[0] = 0.7071067811865475*Jc_Jx_flat[0]-0.7071067811865475*Jc_Jx_flat[1];
  Jc_Jy_flat_n[0] = 0.7071067811865475*Jc_Jy_flat[0]-0.7071067811865475*Jc_Jy_flat[1];
  Jc_Jz_flat_n[0] = 0.7071067811865475*Jc_Jz_flat[0]-0.7071067811865475*Jc_Jz_flat[1];
  
  Jc_rho_n[1] = 0.7071067811865475*Jc_rho[1]+0.7071067811865475*Jc_rho[0];
  Jc_Jx_flat_n[1] = 0.7071067811865475*Jc_Jx_flat[1]+0.7071067811865475*Jc_Jx_flat[0];
  Jc_Jy_flat_n[1] = 0.7071067811865475*Jc_Jy_flat[1]+0.7071067811865475*Jc_Jy_flat[0];
  Jc_Jz_flat_n[1] = 0.7071067811865475*Jc_Jz_flat[1]+0.7071067811865475*Jc_Jz_flat[0];
  
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
    Current_n[i + 0*2] = q_over_eps0 * (Jc_rho_n[i] * shift_nodal_x[i] - lapse_nodal[i] * ( vierb_xx_con_nodal[i] * Jc_Jx_flat_n[i] + vierb_xy_con_nodal[i] * Jc_Jy_flat_n[i] + vierb_xz_con_nodal[i] * Jc_Jz_flat_n[i] )); 
    Current_n[i + 1*2] = q_over_eps0 * (Jc_rho_n[i] * shift_nodal_y[i] - lapse_nodal[i] * ( vierb_yx_con_nodal[i] * Jc_Jx_flat_n[i] + vierb_yy_con_nodal[i] * Jc_Jy_flat_n[i] + vierb_yz_con_nodal[i] * Jc_Jz_flat_n[i] )); 
    Current_n[i + 2*2] = q_over_eps0 * (Jc_rho_n[i] * shift_nodal_z[i] - lapse_nodal[i] * ( vierb_zx_con_nodal[i] * Jc_Jx_flat_n[i] + vierb_zy_con_nodal[i] * Jc_Jy_flat_n[i] + vierb_zz_con_nodal[i] * Jc_Jz_flat_n[i] )); 
  }
  
  double *Current_comp;
  Current_comp = &out[0];
  Current_comp[0] += 0.7071067811865475*Current_n[1]+0.7071067811865475*Current_n[0]; 
  Current_comp[1] += 0.7071067811865475*Current_n[1]-0.7071067811865475*Current_n[0]; 

  Current_comp = &out[2];
  Current_comp[0] += 0.7071067811865475*Current_n[3]+0.7071067811865475*Current_n[2]; 
  Current_comp[1] += 0.7071067811865475*Current_n[3]-0.7071067811865475*Current_n[2]; 

  Current_comp = &out[4];
  Current_comp[0] += 0.7071067811865475*Current_n[5]+0.7071067811865475*Current_n[4]; 
  Current_comp[1] += 0.7071067811865475*Current_n[5]-0.7071067811865475*Current_n[4]; 

}
