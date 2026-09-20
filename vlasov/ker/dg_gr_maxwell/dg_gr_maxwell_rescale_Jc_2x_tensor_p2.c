#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_tensor_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields) 
{ 
  const double *f_comp_no_J;
  double *Jf_comp;
  for (int i=0; i<8; ++i) {
    f_comp_no_J = &fields_no_J[i*9];
    Jf_comp = &Jfields[i*9];
    double Jfields_nodal[9] = {0.0};
    
    Jfields_nodal[0] = J[0]*(0.4*f_comp_no_J[8]-0.5999999999999994*f_comp_no_J[7]-0.5999999999999999*f_comp_no_J[6]+0.4472135954999579*f_comp_no_J[5]+0.4472135954999579*f_comp_no_J[4]+0.9*f_comp_no_J[3]-0.6708203932499369*f_comp_no_J[2]-0.6708203932499369*f_comp_no_J[1]+0.5*f_comp_no_J[0]);
    Jfields_nodal[1] = J[1]*(-(0.5*f_comp_no_J[8])+0.75*f_comp_no_J[7]-0.5590169943749475*f_comp_no_J[5]+0.4472135954999579*f_comp_no_J[4]-0.6708203932499369*f_comp_no_J[1]+0.5*f_comp_no_J[0]);
    Jfields_nodal[2] = J[2]*(0.4*f_comp_no_J[8]-0.5999999999999994*f_comp_no_J[7]+0.5999999999999999*f_comp_no_J[6]+0.4472135954999579*f_comp_no_J[5]+0.4472135954999579*f_comp_no_J[4]-0.9*f_comp_no_J[3]+0.6708203932499369*f_comp_no_J[2]-0.6708203932499369*f_comp_no_J[1]+0.5*f_comp_no_J[0]);
    Jfields_nodal[3] = J[3]*(-(0.5*f_comp_no_J[8])+0.75*f_comp_no_J[6]+0.4472135954999579*f_comp_no_J[5]-0.5590169943749475*f_comp_no_J[4]-0.6708203932499369*f_comp_no_J[2]+0.5*f_comp_no_J[0]);
    Jfields_nodal[4] = J[4]*(0.625*f_comp_no_J[8]-0.5590169943749475*f_comp_no_J[5]-0.5590169943749475*f_comp_no_J[4]+0.5*f_comp_no_J[0]);
    Jfields_nodal[5] = J[5]*(-(0.5*f_comp_no_J[8])-0.75*f_comp_no_J[6]+0.4472135954999579*f_comp_no_J[5]-0.5590169943749475*f_comp_no_J[4]+0.6708203932499369*f_comp_no_J[2]+0.5*f_comp_no_J[0]);
    Jfields_nodal[6] = J[6]*(0.4*f_comp_no_J[8]+0.5999999999999994*f_comp_no_J[7]-0.5999999999999999*f_comp_no_J[6]+0.4472135954999579*f_comp_no_J[5]+0.4472135954999579*f_comp_no_J[4]-0.9*f_comp_no_J[3]-0.6708203932499369*f_comp_no_J[2]+0.6708203932499369*f_comp_no_J[1]+0.5*f_comp_no_J[0]);
    Jfields_nodal[7] = J[7]*(-(0.5*f_comp_no_J[8])-0.75*f_comp_no_J[7]-0.5590169943749475*f_comp_no_J[5]+0.4472135954999579*f_comp_no_J[4]+0.6708203932499369*f_comp_no_J[1]+0.5*f_comp_no_J[0]);
    Jfields_nodal[8] = J[8]*(0.4*f_comp_no_J[8]+0.5999999999999994*f_comp_no_J[7]+0.5999999999999999*f_comp_no_J[6]+0.4472135954999579*f_comp_no_J[5]+0.4472135954999579*f_comp_no_J[4]+0.9*f_comp_no_J[3]+0.6708203932499369*f_comp_no_J[2]+0.6708203932499369*f_comp_no_J[1]+0.5*f_comp_no_J[0]);
    
  Jf_comp[0] = 0.15432098765432098*Jfields_nodal[8]+0.24691358024691357*Jfields_nodal[7]+0.15432098765432098*Jfields_nodal[6]+0.24691358024691357*Jfields_nodal[5]+0.3950617283950617*Jfields_nodal[4]+0.24691358024691357*Jfields_nodal[3]+0.15432098765432098*Jfields_nodal[2]+0.24691358024691357*Jfields_nodal[1]+0.15432098765432098*Jfields_nodal[0]; 
  Jf_comp[1] = 0.20704333124998056*Jfields_nodal[8]+0.33126932999996883*Jfields_nodal[7]+0.20704333124998056*Jfields_nodal[6]-0.20704333124998056*Jfields_nodal[2]-0.33126932999996883*Jfields_nodal[1]-0.20704333124998056*Jfields_nodal[0]; 
  Jf_comp[2] = 0.20704333124998056*Jfields_nodal[8]-0.20704333124998056*Jfields_nodal[6]+0.33126932999996883*Jfields_nodal[5]-0.33126932999996883*Jfields_nodal[3]+0.20704333124998056*Jfields_nodal[2]-0.20704333124998056*Jfields_nodal[0]; 
  Jf_comp[3] = 0.2777777777777778*Jfields_nodal[8]-0.2777777777777778*Jfields_nodal[6]-0.2777777777777778*Jfields_nodal[2]+0.2777777777777778*Jfields_nodal[0]; 
  Jf_comp[4] = 0.13802888749998704*Jfields_nodal[8]+0.22084621999997922*Jfields_nodal[7]+0.13802888749998704*Jfields_nodal[6]-0.2760577749999741*Jfields_nodal[5]-0.44169243999995844*Jfields_nodal[4]-0.2760577749999741*Jfields_nodal[3]+0.13802888749998704*Jfields_nodal[2]+0.22084621999997922*Jfields_nodal[1]+0.13802888749998704*Jfields_nodal[0]; 
  Jf_comp[5] = 0.13802888749998704*Jfields_nodal[8]-0.2760577749999741*Jfields_nodal[7]+0.13802888749998704*Jfields_nodal[6]+0.22084621999997922*Jfields_nodal[5]-0.44169243999995844*Jfields_nodal[4]+0.22084621999997922*Jfields_nodal[3]+0.13802888749998704*Jfields_nodal[2]-0.2760577749999741*Jfields_nodal[1]+0.13802888749998704*Jfields_nodal[0]; 
  Jf_comp[6] = 0.18518518518518526*Jfields_nodal[8]-0.18518518518518526*Jfields_nodal[6]-0.3703703703703705*Jfields_nodal[5]+0.3703703703703705*Jfields_nodal[3]+0.18518518518518526*Jfields_nodal[2]-0.18518518518518526*Jfields_nodal[0]; 
  Jf_comp[7] = 0.18518518518518526*Jfields_nodal[8]-0.3703703703703705*Jfields_nodal[7]+0.18518518518518526*Jfields_nodal[6]-0.18518518518518526*Jfields_nodal[2]+0.3703703703703705*Jfields_nodal[1]-0.18518518518518526*Jfields_nodal[0]; 
  Jf_comp[8] = 0.12345679012345678*Jfields_nodal[8]-0.24691358024691357*Jfields_nodal[7]+0.12345679012345678*Jfields_nodal[6]-0.24691358024691357*Jfields_nodal[5]+0.49382716049382713*Jfields_nodal[4]-0.24691358024691357*Jfields_nodal[3]+0.12345679012345678*Jfields_nodal[2]-0.24691358024691357*Jfields_nodal[1]+0.12345679012345678*Jfields_nodal[0]; 
  }
} 
