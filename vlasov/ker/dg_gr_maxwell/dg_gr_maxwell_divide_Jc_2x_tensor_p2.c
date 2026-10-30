#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_2x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J) 
{ 
  const double *Jf_comp;
  double *f_comp_no_J;
  for (int i=0; i<8; ++i) {
  Jf_comp = &Jfields[i*9];
  f_comp_no_J = &fields_no_J[i*9];
    double fields_no_J_nodal[9] = {0.0};
    
    fields_no_J_nodal[0] = 1.0/J[0]*(0.4*Jf_comp[8]-0.5999999999999994*Jf_comp[7]-0.5999999999999999*Jf_comp[6]+0.4472135954999579*Jf_comp[5]+0.4472135954999579*Jf_comp[4]+0.9*Jf_comp[3]-0.6708203932499369*Jf_comp[2]-0.6708203932499369*Jf_comp[1]+0.5*Jf_comp[0]);
    fields_no_J_nodal[1] = 1.0/J[1]*(-(0.5*Jf_comp[8])+0.75*Jf_comp[7]-0.5590169943749475*Jf_comp[5]+0.4472135954999579*Jf_comp[4]-0.6708203932499369*Jf_comp[1]+0.5*Jf_comp[0]);
    fields_no_J_nodal[2] = 1.0/J[2]*(0.4*Jf_comp[8]-0.5999999999999994*Jf_comp[7]+0.5999999999999999*Jf_comp[6]+0.4472135954999579*Jf_comp[5]+0.4472135954999579*Jf_comp[4]-0.9*Jf_comp[3]+0.6708203932499369*Jf_comp[2]-0.6708203932499369*Jf_comp[1]+0.5*Jf_comp[0]);
    fields_no_J_nodal[3] = 1.0/J[3]*(-(0.5*Jf_comp[8])+0.75*Jf_comp[6]+0.4472135954999579*Jf_comp[5]-0.5590169943749475*Jf_comp[4]-0.6708203932499369*Jf_comp[2]+0.5*Jf_comp[0]);
    fields_no_J_nodal[4] = 1.0/J[4]*(0.625*Jf_comp[8]-0.5590169943749475*Jf_comp[5]-0.5590169943749475*Jf_comp[4]+0.5*Jf_comp[0]);
    fields_no_J_nodal[5] = 1.0/J[5]*(-(0.5*Jf_comp[8])-0.75*Jf_comp[6]+0.4472135954999579*Jf_comp[5]-0.5590169943749475*Jf_comp[4]+0.6708203932499369*Jf_comp[2]+0.5*Jf_comp[0]);
    fields_no_J_nodal[6] = 1.0/J[6]*(0.4*Jf_comp[8]+0.5999999999999994*Jf_comp[7]-0.5999999999999999*Jf_comp[6]+0.4472135954999579*Jf_comp[5]+0.4472135954999579*Jf_comp[4]-0.9*Jf_comp[3]-0.6708203932499369*Jf_comp[2]+0.6708203932499369*Jf_comp[1]+0.5*Jf_comp[0]);
    fields_no_J_nodal[7] = 1.0/J[7]*(-(0.5*Jf_comp[8])-0.75*Jf_comp[7]-0.5590169943749475*Jf_comp[5]+0.4472135954999579*Jf_comp[4]+0.6708203932499369*Jf_comp[1]+0.5*Jf_comp[0]);
    fields_no_J_nodal[8] = 1.0/J[8]*(0.4*Jf_comp[8]+0.5999999999999994*Jf_comp[7]+0.5999999999999999*Jf_comp[6]+0.4472135954999579*Jf_comp[5]+0.4472135954999579*Jf_comp[4]+0.9*Jf_comp[3]+0.6708203932499369*Jf_comp[2]+0.6708203932499369*Jf_comp[1]+0.5*Jf_comp[0]);
  
  f_comp_no_J[0] = 0.15432098765432098*fields_no_J_nodal[8]+0.24691358024691357*fields_no_J_nodal[7]+0.15432098765432098*fields_no_J_nodal[6]+0.24691358024691357*fields_no_J_nodal[5]+0.3950617283950617*fields_no_J_nodal[4]+0.24691358024691357*fields_no_J_nodal[3]+0.15432098765432098*fields_no_J_nodal[2]+0.24691358024691357*fields_no_J_nodal[1]+0.15432098765432098*fields_no_J_nodal[0]; 
  f_comp_no_J[1] = 0.20704333124998056*fields_no_J_nodal[8]+0.33126932999996883*fields_no_J_nodal[7]+0.20704333124998056*fields_no_J_nodal[6]-0.20704333124998056*fields_no_J_nodal[2]-0.33126932999996883*fields_no_J_nodal[1]-0.20704333124998056*fields_no_J_nodal[0]; 
  f_comp_no_J[2] = 0.20704333124998056*fields_no_J_nodal[8]-0.20704333124998056*fields_no_J_nodal[6]+0.33126932999996883*fields_no_J_nodal[5]-0.33126932999996883*fields_no_J_nodal[3]+0.20704333124998056*fields_no_J_nodal[2]-0.20704333124998056*fields_no_J_nodal[0]; 
  f_comp_no_J[3] = 0.2777777777777778*fields_no_J_nodal[8]-0.2777777777777778*fields_no_J_nodal[6]-0.2777777777777778*fields_no_J_nodal[2]+0.2777777777777778*fields_no_J_nodal[0]; 
  f_comp_no_J[4] = 0.13802888749998704*fields_no_J_nodal[8]+0.22084621999997922*fields_no_J_nodal[7]+0.13802888749998704*fields_no_J_nodal[6]-0.2760577749999741*fields_no_J_nodal[5]-0.44169243999995844*fields_no_J_nodal[4]-0.2760577749999741*fields_no_J_nodal[3]+0.13802888749998704*fields_no_J_nodal[2]+0.22084621999997922*fields_no_J_nodal[1]+0.13802888749998704*fields_no_J_nodal[0]; 
  f_comp_no_J[5] = 0.13802888749998704*fields_no_J_nodal[8]-0.2760577749999741*fields_no_J_nodal[7]+0.13802888749998704*fields_no_J_nodal[6]+0.22084621999997922*fields_no_J_nodal[5]-0.44169243999995844*fields_no_J_nodal[4]+0.22084621999997922*fields_no_J_nodal[3]+0.13802888749998704*fields_no_J_nodal[2]-0.2760577749999741*fields_no_J_nodal[1]+0.13802888749998704*fields_no_J_nodal[0]; 
  f_comp_no_J[6] = 0.18518518518518526*fields_no_J_nodal[8]-0.18518518518518526*fields_no_J_nodal[6]-0.3703703703703705*fields_no_J_nodal[5]+0.3703703703703705*fields_no_J_nodal[3]+0.18518518518518526*fields_no_J_nodal[2]-0.18518518518518526*fields_no_J_nodal[0]; 
  f_comp_no_J[7] = 0.18518518518518526*fields_no_J_nodal[8]-0.3703703703703705*fields_no_J_nodal[7]+0.18518518518518526*fields_no_J_nodal[6]-0.18518518518518526*fields_no_J_nodal[2]+0.3703703703703705*fields_no_J_nodal[1]-0.18518518518518526*fields_no_J_nodal[0]; 
  f_comp_no_J[8] = 0.12345679012345678*fields_no_J_nodal[8]-0.24691358024691357*fields_no_J_nodal[7]+0.12345679012345678*fields_no_J_nodal[6]-0.24691358024691357*fields_no_J_nodal[5]+0.49382716049382713*fields_no_J_nodal[4]-0.24691358024691357*fields_no_J_nodal[3]+0.12345679012345678*fields_no_J_nodal[2]-0.24691358024691357*fields_no_J_nodal[1]+0.12345679012345678*fields_no_J_nodal[0]; 
  }
} 
