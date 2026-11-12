#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void gr_maxwell_current_dep_2x_3v_ser_p2(const double q_over_eps0, const double *lapse_nodal, const double *shift_nodal,
  const double *vierb_con_nodal, const double *m0, const double *m1i, double* GKYL_RESTRICT out) 
{ 
  
  double Jc_rho_n[9] = {0.0};
  double Jc_Jx_flat_n[9] = {0.0};
  double Jc_Jy_flat_n[9] = {0.0};
  double Jc_Jz_flat_n[9] = {0.0};
  
  const double *Jc_rho = &m0[0]; 
  const double *Jc_Jx_flat = &m1i[0]; 
  const double *Jc_Jy_flat = &m1i[8]; 
  const double *Jc_Jz_flat = &m1i[16]; 
  
  double Current_n[27] = {0.0};
  
  Jc_rho_n[0] = -(0.5999999999999995*Jc_rho[7])-0.5999999999999999*Jc_rho[6]+0.4472135954999579*Jc_rho[5]+0.4472135954999579*Jc_rho[4]+0.9*Jc_rho[3]-0.6708203932499369*Jc_rho[2]-0.6708203932499369*Jc_rho[1]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[0] = -(0.5999999999999995*Jc_Jx_flat[7])-0.5999999999999999*Jc_Jx_flat[6]+0.4472135954999579*Jc_Jx_flat[5]+0.4472135954999579*Jc_Jx_flat[4]+0.9*Jc_Jx_flat[3]-0.6708203932499369*Jc_Jx_flat[2]-0.6708203932499369*Jc_Jx_flat[1]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[0] = -(0.5999999999999995*Jc_Jy_flat[7])-0.5999999999999999*Jc_Jy_flat[6]+0.4472135954999579*Jc_Jy_flat[5]+0.4472135954999579*Jc_Jy_flat[4]+0.9*Jc_Jy_flat[3]-0.6708203932499369*Jc_Jy_flat[2]-0.6708203932499369*Jc_Jy_flat[1]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[0] = -(0.5999999999999995*Jc_Jz_flat[7])-0.5999999999999999*Jc_Jz_flat[6]+0.4472135954999579*Jc_Jz_flat[5]+0.4472135954999579*Jc_Jz_flat[4]+0.9*Jc_Jz_flat[3]-0.6708203932499369*Jc_Jz_flat[2]-0.6708203932499369*Jc_Jz_flat[1]+0.5*Jc_Jz_flat[0];
  
  Jc_rho_n[1] = 0.75*Jc_rho[7]-0.5590169943749475*Jc_rho[5]+0.4472135954999579*Jc_rho[4]-0.6708203932499369*Jc_rho[1]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[1] = 0.75*Jc_Jx_flat[7]-0.5590169943749475*Jc_Jx_flat[5]+0.4472135954999579*Jc_Jx_flat[4]-0.6708203932499369*Jc_Jx_flat[1]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[1] = 0.75*Jc_Jy_flat[7]-0.5590169943749475*Jc_Jy_flat[5]+0.4472135954999579*Jc_Jy_flat[4]-0.6708203932499369*Jc_Jy_flat[1]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[1] = 0.75*Jc_Jz_flat[7]-0.5590169943749475*Jc_Jz_flat[5]+0.4472135954999579*Jc_Jz_flat[4]-0.6708203932499369*Jc_Jz_flat[1]+0.5*Jc_Jz_flat[0];
  
  Jc_rho_n[2] = -(0.5999999999999995*Jc_rho[7])+0.5999999999999999*Jc_rho[6]+0.4472135954999579*Jc_rho[5]+0.4472135954999579*Jc_rho[4]-0.9*Jc_rho[3]+0.6708203932499369*Jc_rho[2]-0.6708203932499369*Jc_rho[1]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[2] = -(0.5999999999999995*Jc_Jx_flat[7])+0.5999999999999999*Jc_Jx_flat[6]+0.4472135954999579*Jc_Jx_flat[5]+0.4472135954999579*Jc_Jx_flat[4]-0.9*Jc_Jx_flat[3]+0.6708203932499369*Jc_Jx_flat[2]-0.6708203932499369*Jc_Jx_flat[1]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[2] = -(0.5999999999999995*Jc_Jy_flat[7])+0.5999999999999999*Jc_Jy_flat[6]+0.4472135954999579*Jc_Jy_flat[5]+0.4472135954999579*Jc_Jy_flat[4]-0.9*Jc_Jy_flat[3]+0.6708203932499369*Jc_Jy_flat[2]-0.6708203932499369*Jc_Jy_flat[1]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[2] = -(0.5999999999999995*Jc_Jz_flat[7])+0.5999999999999999*Jc_Jz_flat[6]+0.4472135954999579*Jc_Jz_flat[5]+0.4472135954999579*Jc_Jz_flat[4]-0.9*Jc_Jz_flat[3]+0.6708203932499369*Jc_Jz_flat[2]-0.6708203932499369*Jc_Jz_flat[1]+0.5*Jc_Jz_flat[0];
  
  Jc_rho_n[3] = 0.75*Jc_rho[6]+0.4472135954999579*Jc_rho[5]-0.5590169943749475*Jc_rho[4]-0.6708203932499369*Jc_rho[2]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[3] = 0.75*Jc_Jx_flat[6]+0.4472135954999579*Jc_Jx_flat[5]-0.5590169943749475*Jc_Jx_flat[4]-0.6708203932499369*Jc_Jx_flat[2]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[3] = 0.75*Jc_Jy_flat[6]+0.4472135954999579*Jc_Jy_flat[5]-0.5590169943749475*Jc_Jy_flat[4]-0.6708203932499369*Jc_Jy_flat[2]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[3] = 0.75*Jc_Jz_flat[6]+0.4472135954999579*Jc_Jz_flat[5]-0.5590169943749475*Jc_Jz_flat[4]-0.6708203932499369*Jc_Jz_flat[2]+0.5*Jc_Jz_flat[0];
  
  Jc_rho_n[4] = -(0.5590169943749475*Jc_rho[5])-0.5590169943749475*Jc_rho[4]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[4] = -(0.5590169943749475*Jc_Jx_flat[5])-0.5590169943749475*Jc_Jx_flat[4]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[4] = -(0.5590169943749475*Jc_Jy_flat[5])-0.5590169943749475*Jc_Jy_flat[4]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[4] = -(0.5590169943749475*Jc_Jz_flat[5])-0.5590169943749475*Jc_Jz_flat[4]+0.5*Jc_Jz_flat[0];
  
  Jc_rho_n[5] = -(0.75*Jc_rho[6])+0.4472135954999579*Jc_rho[5]-0.5590169943749475*Jc_rho[4]+0.6708203932499369*Jc_rho[2]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[5] = -(0.75*Jc_Jx_flat[6])+0.4472135954999579*Jc_Jx_flat[5]-0.5590169943749475*Jc_Jx_flat[4]+0.6708203932499369*Jc_Jx_flat[2]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[5] = -(0.75*Jc_Jy_flat[6])+0.4472135954999579*Jc_Jy_flat[5]-0.5590169943749475*Jc_Jy_flat[4]+0.6708203932499369*Jc_Jy_flat[2]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[5] = -(0.75*Jc_Jz_flat[6])+0.4472135954999579*Jc_Jz_flat[5]-0.5590169943749475*Jc_Jz_flat[4]+0.6708203932499369*Jc_Jz_flat[2]+0.5*Jc_Jz_flat[0];
  
  Jc_rho_n[6] = 0.5999999999999995*Jc_rho[7]-0.5999999999999999*Jc_rho[6]+0.4472135954999579*Jc_rho[5]+0.4472135954999579*Jc_rho[4]-0.9*Jc_rho[3]-0.6708203932499369*Jc_rho[2]+0.6708203932499369*Jc_rho[1]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[6] = 0.5999999999999995*Jc_Jx_flat[7]-0.5999999999999999*Jc_Jx_flat[6]+0.4472135954999579*Jc_Jx_flat[5]+0.4472135954999579*Jc_Jx_flat[4]-0.9*Jc_Jx_flat[3]-0.6708203932499369*Jc_Jx_flat[2]+0.6708203932499369*Jc_Jx_flat[1]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[6] = 0.5999999999999995*Jc_Jy_flat[7]-0.5999999999999999*Jc_Jy_flat[6]+0.4472135954999579*Jc_Jy_flat[5]+0.4472135954999579*Jc_Jy_flat[4]-0.9*Jc_Jy_flat[3]-0.6708203932499369*Jc_Jy_flat[2]+0.6708203932499369*Jc_Jy_flat[1]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[6] = 0.5999999999999995*Jc_Jz_flat[7]-0.5999999999999999*Jc_Jz_flat[6]+0.4472135954999579*Jc_Jz_flat[5]+0.4472135954999579*Jc_Jz_flat[4]-0.9*Jc_Jz_flat[3]-0.6708203932499369*Jc_Jz_flat[2]+0.6708203932499369*Jc_Jz_flat[1]+0.5*Jc_Jz_flat[0];
  
  Jc_rho_n[7] = -(0.75*Jc_rho[7])-0.5590169943749475*Jc_rho[5]+0.4472135954999579*Jc_rho[4]+0.6708203932499369*Jc_rho[1]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[7] = -(0.75*Jc_Jx_flat[7])-0.5590169943749475*Jc_Jx_flat[5]+0.4472135954999579*Jc_Jx_flat[4]+0.6708203932499369*Jc_Jx_flat[1]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[7] = -(0.75*Jc_Jy_flat[7])-0.5590169943749475*Jc_Jy_flat[5]+0.4472135954999579*Jc_Jy_flat[4]+0.6708203932499369*Jc_Jy_flat[1]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[7] = -(0.75*Jc_Jz_flat[7])-0.5590169943749475*Jc_Jz_flat[5]+0.4472135954999579*Jc_Jz_flat[4]+0.6708203932499369*Jc_Jz_flat[1]+0.5*Jc_Jz_flat[0];
  
  Jc_rho_n[8] = 0.5999999999999995*Jc_rho[7]+0.5999999999999999*Jc_rho[6]+0.4472135954999579*Jc_rho[5]+0.4472135954999579*Jc_rho[4]+0.9*Jc_rho[3]+0.6708203932499369*Jc_rho[2]+0.6708203932499369*Jc_rho[1]+0.5*Jc_rho[0];
  Jc_Jx_flat_n[8] = 0.5999999999999995*Jc_Jx_flat[7]+0.5999999999999999*Jc_Jx_flat[6]+0.4472135954999579*Jc_Jx_flat[5]+0.4472135954999579*Jc_Jx_flat[4]+0.9*Jc_Jx_flat[3]+0.6708203932499369*Jc_Jx_flat[2]+0.6708203932499369*Jc_Jx_flat[1]+0.5*Jc_Jx_flat[0];
  Jc_Jy_flat_n[8] = 0.5999999999999995*Jc_Jy_flat[7]+0.5999999999999999*Jc_Jy_flat[6]+0.4472135954999579*Jc_Jy_flat[5]+0.4472135954999579*Jc_Jy_flat[4]+0.9*Jc_Jy_flat[3]+0.6708203932499369*Jc_Jy_flat[2]+0.6708203932499369*Jc_Jy_flat[1]+0.5*Jc_Jy_flat[0];
  Jc_Jz_flat_n[8] = 0.5999999999999995*Jc_Jz_flat[7]+0.5999999999999999*Jc_Jz_flat[6]+0.4472135954999579*Jc_Jz_flat[5]+0.4472135954999579*Jc_Jz_flat[4]+0.9*Jc_Jz_flat[3]+0.6708203932499369*Jc_Jz_flat[2]+0.6708203932499369*Jc_Jz_flat[1]+0.5*Jc_Jz_flat[0];
  
  const double *vierb_xx_con_nodal = &vierb_con_nodal[0]; 
  const double *vierb_xy_con_nodal = &vierb_con_nodal[9]; 
  const double *vierb_xz_con_nodal = &vierb_con_nodal[18]; 
  const double *vierb_yx_con_nodal = &vierb_con_nodal[27]; 
  const double *vierb_yy_con_nodal = &vierb_con_nodal[36]; 
  const double *vierb_yz_con_nodal = &vierb_con_nodal[45]; 
  const double *vierb_zx_con_nodal = &vierb_con_nodal[54]; 
  const double *vierb_zy_con_nodal = &vierb_con_nodal[63]; 
  const double *vierb_zz_con_nodal = &vierb_con_nodal[72]; 
  const double *shift_nodal_x = &shift_nodal[0]; 
  const double *shift_nodal_y = &shift_nodal[9]; 
  const double *shift_nodal_z = &shift_nodal[18]; 
  
  for (int i=0; i<9; ++i) {
    Current_n[i + 0*9] = q_over_eps0 * (Jc_rho_n[i] * shift_nodal_x[i] - lapse_nodal[i] * ( vierb_xx_con_nodal[i] * Jc_Jx_flat_n[i] + vierb_xy_con_nodal[i] * Jc_Jy_flat_n[i] + vierb_xz_con_nodal[i] * Jc_Jz_flat_n[i] )); 
    Current_n[i + 1*9] = q_over_eps0 * (Jc_rho_n[i] * shift_nodal_y[i] - lapse_nodal[i] * ( vierb_yx_con_nodal[i] * Jc_Jx_flat_n[i] + vierb_yy_con_nodal[i] * Jc_Jy_flat_n[i] + vierb_yz_con_nodal[i] * Jc_Jz_flat_n[i] )); 
    Current_n[i + 2*9] = q_over_eps0 * (Jc_rho_n[i] * shift_nodal_z[i] - lapse_nodal[i] * ( vierb_zx_con_nodal[i] * Jc_Jx_flat_n[i] + vierb_zy_con_nodal[i] * Jc_Jy_flat_n[i] + vierb_zz_con_nodal[i] * Jc_Jz_flat_n[i] )); 
  }
  
  double *Current_comp;
  Current_comp = &out[0];
  Current_comp[0] += 0.15432098765432098*Current_n[8]+0.24691358024691357*Current_n[7]+0.15432098765432098*Current_n[6]+0.24691358024691357*Current_n[5]+0.3950617283950617*Current_n[4]+0.24691358024691357*Current_n[3]+0.15432098765432098*Current_n[2]+0.24691358024691357*Current_n[1]+0.15432098765432098*Current_n[0]; 
  Current_comp[1] += 0.20704333124998056*Current_n[8]+0.33126932999996883*Current_n[7]+0.20704333124998056*Current_n[6]-0.20704333124998056*Current_n[2]-0.33126932999996883*Current_n[1]-0.20704333124998056*Current_n[0]; 
  Current_comp[2] += 0.20704333124998056*Current_n[8]-0.20704333124998056*Current_n[6]+0.33126932999996883*Current_n[5]-0.33126932999996883*Current_n[3]+0.20704333124998056*Current_n[2]-0.20704333124998056*Current_n[0]; 
  Current_comp[3] += 0.2777777777777778*Current_n[8]-0.2777777777777778*Current_n[6]-0.2777777777777778*Current_n[2]+0.2777777777777778*Current_n[0]; 
  Current_comp[4] += 0.13802888749998704*Current_n[8]+0.22084621999997922*Current_n[7]+0.13802888749998704*Current_n[6]-0.2760577749999741*Current_n[5]-0.44169243999995844*Current_n[4]-0.2760577749999741*Current_n[3]+0.13802888749998704*Current_n[2]+0.22084621999997922*Current_n[1]+0.13802888749998704*Current_n[0]; 
  Current_comp[5] += 0.13802888749998704*Current_n[8]-0.2760577749999741*Current_n[7]+0.13802888749998704*Current_n[6]+0.22084621999997922*Current_n[5]-0.44169243999995844*Current_n[4]+0.22084621999997922*Current_n[3]+0.13802888749998704*Current_n[2]-0.2760577749999741*Current_n[1]+0.13802888749998704*Current_n[0]; 
  Current_comp[6] += 0.18518518518518526*Current_n[8]-0.18518518518518526*Current_n[6]-0.3703703703703705*Current_n[5]+0.3703703703703705*Current_n[3]+0.18518518518518526*Current_n[2]-0.18518518518518526*Current_n[0]; 
  Current_comp[7] += 0.18518518518518526*Current_n[8]-0.3703703703703705*Current_n[7]+0.18518518518518526*Current_n[6]-0.18518518518518526*Current_n[2]+0.3703703703703705*Current_n[1]-0.18518518518518526*Current_n[0]; 

  Current_comp = &out[8];
  Current_comp[0] += 0.15432098765432098*Current_n[17]+0.24691358024691357*Current_n[16]+0.15432098765432098*Current_n[15]+0.24691358024691357*Current_n[14]+0.3950617283950617*Current_n[13]+0.24691358024691357*Current_n[12]+0.15432098765432098*Current_n[11]+0.24691358024691357*Current_n[10]+0.15432098765432098*Current_n[9]; 
  Current_comp[1] += 0.20704333124998056*Current_n[17]+0.33126932999996883*Current_n[16]+0.20704333124998056*Current_n[15]-0.20704333124998056*Current_n[11]-0.33126932999996883*Current_n[10]-0.20704333124998056*Current_n[9]; 
  Current_comp[2] += 0.20704333124998056*Current_n[17]-0.20704333124998056*Current_n[15]+0.33126932999996883*Current_n[14]-0.33126932999996883*Current_n[12]+0.20704333124998056*Current_n[11]-0.20704333124998056*Current_n[9]; 
  Current_comp[3] += 0.2777777777777778*Current_n[17]-0.2777777777777778*Current_n[15]-0.2777777777777778*Current_n[11]+0.2777777777777778*Current_n[9]; 
  Current_comp[4] += 0.13802888749998704*Current_n[17]+0.22084621999997922*Current_n[16]+0.13802888749998704*Current_n[15]-0.2760577749999741*Current_n[14]-0.44169243999995844*Current_n[13]-0.2760577749999741*Current_n[12]+0.13802888749998704*Current_n[11]+0.22084621999997922*Current_n[10]+0.13802888749998704*Current_n[9]; 
  Current_comp[5] += 0.13802888749998704*Current_n[17]-0.2760577749999741*Current_n[16]+0.13802888749998704*Current_n[15]+0.22084621999997922*Current_n[14]-0.44169243999995844*Current_n[13]+0.22084621999997922*Current_n[12]+0.13802888749998704*Current_n[11]-0.2760577749999741*Current_n[10]+0.13802888749998704*Current_n[9]; 
  Current_comp[6] += 0.18518518518518526*Current_n[17]-0.18518518518518526*Current_n[15]-0.3703703703703705*Current_n[14]+0.3703703703703705*Current_n[12]+0.18518518518518526*Current_n[11]-0.18518518518518526*Current_n[9]; 
  Current_comp[7] += 0.18518518518518526*Current_n[17]-0.3703703703703705*Current_n[16]+0.18518518518518526*Current_n[15]-0.18518518518518526*Current_n[11]+0.3703703703703705*Current_n[10]-0.18518518518518526*Current_n[9]; 

  Current_comp = &out[16];
  Current_comp[0] += 0.15432098765432098*Current_n[26]+0.24691358024691357*Current_n[25]+0.15432098765432098*Current_n[24]+0.24691358024691357*Current_n[23]+0.3950617283950617*Current_n[22]+0.24691358024691357*Current_n[21]+0.15432098765432098*Current_n[20]+0.24691358024691357*Current_n[19]+0.15432098765432098*Current_n[18]; 
  Current_comp[1] += 0.20704333124998056*Current_n[26]+0.33126932999996883*Current_n[25]+0.20704333124998056*Current_n[24]-0.20704333124998056*Current_n[20]-0.33126932999996883*Current_n[19]-0.20704333124998056*Current_n[18]; 
  Current_comp[2] += 0.20704333124998056*Current_n[26]-0.20704333124998056*Current_n[24]+0.33126932999996883*Current_n[23]-0.33126932999996883*Current_n[21]+0.20704333124998056*Current_n[20]-0.20704333124998056*Current_n[18]; 
  Current_comp[3] += 0.2777777777777778*Current_n[26]-0.2777777777777778*Current_n[24]-0.2777777777777778*Current_n[20]+0.2777777777777778*Current_n[18]; 
  Current_comp[4] += 0.13802888749998704*Current_n[26]+0.22084621999997922*Current_n[25]+0.13802888749998704*Current_n[24]-0.2760577749999741*Current_n[23]-0.44169243999995844*Current_n[22]-0.2760577749999741*Current_n[21]+0.13802888749998704*Current_n[20]+0.22084621999997922*Current_n[19]+0.13802888749998704*Current_n[18]; 
  Current_comp[5] += 0.13802888749998704*Current_n[26]-0.2760577749999741*Current_n[25]+0.13802888749998704*Current_n[24]+0.22084621999997922*Current_n[23]-0.44169243999995844*Current_n[22]+0.22084621999997922*Current_n[21]+0.13802888749998704*Current_n[20]-0.2760577749999741*Current_n[19]+0.13802888749998704*Current_n[18]; 
  Current_comp[6] += 0.18518518518518526*Current_n[26]-0.18518518518518526*Current_n[24]-0.3703703703703705*Current_n[23]+0.3703703703703705*Current_n[21]+0.18518518518518526*Current_n[20]-0.18518518518518526*Current_n[18]; 
  Current_comp[7] += 0.18518518518518526*Current_n[26]-0.3703703703703705*Current_n[25]+0.18518518518518526*Current_n[24]-0.18518518518518526*Current_n[20]+0.3703703703703705*Current_n[19]-0.18518518518518526*Current_n[18]; 

}
