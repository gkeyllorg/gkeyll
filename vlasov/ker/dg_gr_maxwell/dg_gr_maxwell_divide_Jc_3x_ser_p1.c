#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_3x_ser_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J) 
{ 
  const double *Jf_comp;
  double *f_comp_no_J;
  for (int i=0; i<8; ++i) {
  Jf_comp = &Jfields[i*8];
  f_comp_no_J = &fields_no_J[i*8];
    double fields_no_J_nodal[8] = {0.0};
    
    fields_no_J_nodal[0] = 1.0/J[0]*(-(0.3535533905932737*Jf_comp[7])+0.3535533905932737*Jf_comp[6]+0.3535533905932737*Jf_comp[5]+0.3535533905932737*Jf_comp[4]-0.3535533905932737*Jf_comp[3]-0.3535533905932737*Jf_comp[2]-0.3535533905932737*Jf_comp[1]+0.3535533905932737*Jf_comp[0]);
    fields_no_J_nodal[1] = 1.0/J[1]*(0.3535533905932737*Jf_comp[7]-0.3535533905932737*Jf_comp[6]-0.3535533905932737*Jf_comp[5]+0.3535533905932737*Jf_comp[4]+0.3535533905932737*Jf_comp[3]-0.3535533905932737*Jf_comp[2]-0.3535533905932737*Jf_comp[1]+0.3535533905932737*Jf_comp[0]);
    fields_no_J_nodal[2] = 1.0/J[2]*(0.3535533905932737*Jf_comp[7]-0.3535533905932737*Jf_comp[6]+0.3535533905932737*Jf_comp[5]-0.3535533905932737*Jf_comp[4]-0.3535533905932737*Jf_comp[3]+0.3535533905932737*Jf_comp[2]-0.3535533905932737*Jf_comp[1]+0.3535533905932737*Jf_comp[0]);
    fields_no_J_nodal[3] = 1.0/J[3]*(-(0.3535533905932737*Jf_comp[7])+0.3535533905932737*Jf_comp[6]-0.3535533905932737*Jf_comp[5]-0.3535533905932737*Jf_comp[4]+0.3535533905932737*Jf_comp[3]+0.3535533905932737*Jf_comp[2]-0.3535533905932737*Jf_comp[1]+0.3535533905932737*Jf_comp[0]);
    fields_no_J_nodal[4] = 1.0/J[4]*(0.3535533905932737*Jf_comp[7]+0.3535533905932737*Jf_comp[6]-0.3535533905932737*Jf_comp[5]-0.3535533905932737*Jf_comp[4]-0.3535533905932737*Jf_comp[3]-0.3535533905932737*Jf_comp[2]+0.3535533905932737*Jf_comp[1]+0.3535533905932737*Jf_comp[0]);
    fields_no_J_nodal[5] = 1.0/J[5]*(-(0.3535533905932737*Jf_comp[7])-0.3535533905932737*Jf_comp[6]+0.3535533905932737*Jf_comp[5]-0.3535533905932737*Jf_comp[4]+0.3535533905932737*Jf_comp[3]-0.3535533905932737*Jf_comp[2]+0.3535533905932737*Jf_comp[1]+0.3535533905932737*Jf_comp[0]);
    fields_no_J_nodal[6] = 1.0/J[6]*(-(0.3535533905932737*Jf_comp[7])-0.3535533905932737*Jf_comp[6]-0.3535533905932737*Jf_comp[5]+0.3535533905932737*Jf_comp[4]-0.3535533905932737*Jf_comp[3]+0.3535533905932737*Jf_comp[2]+0.3535533905932737*Jf_comp[1]+0.3535533905932737*Jf_comp[0]);
    fields_no_J_nodal[7] = 1.0/J[7]*(0.3535533905932737*Jf_comp[7]+0.3535533905932737*Jf_comp[6]+0.3535533905932737*Jf_comp[5]+0.3535533905932737*Jf_comp[4]+0.3535533905932737*Jf_comp[3]+0.3535533905932737*Jf_comp[2]+0.3535533905932737*Jf_comp[1]+0.3535533905932737*Jf_comp[0]);
  
  f_comp_no_J[0] = 0.3535533905932737*fields_no_J_nodal[7]+0.3535533905932737*fields_no_J_nodal[6]+0.3535533905932737*fields_no_J_nodal[5]+0.3535533905932737*fields_no_J_nodal[4]+0.3535533905932737*fields_no_J_nodal[3]+0.3535533905932737*fields_no_J_nodal[2]+0.3535533905932737*fields_no_J_nodal[1]+0.3535533905932737*fields_no_J_nodal[0]; 
  f_comp_no_J[1] = 0.3535533905932737*fields_no_J_nodal[7]+0.3535533905932737*fields_no_J_nodal[6]+0.3535533905932737*fields_no_J_nodal[5]+0.3535533905932737*fields_no_J_nodal[4]-0.3535533905932737*fields_no_J_nodal[3]-0.3535533905932737*fields_no_J_nodal[2]-0.3535533905932737*fields_no_J_nodal[1]-0.3535533905932737*fields_no_J_nodal[0]; 
  f_comp_no_J[2] = 0.3535533905932737*fields_no_J_nodal[7]+0.3535533905932737*fields_no_J_nodal[6]-0.3535533905932737*fields_no_J_nodal[5]-0.3535533905932737*fields_no_J_nodal[4]+0.3535533905932737*fields_no_J_nodal[3]+0.3535533905932737*fields_no_J_nodal[2]-0.3535533905932737*fields_no_J_nodal[1]-0.3535533905932737*fields_no_J_nodal[0]; 
  f_comp_no_J[3] = 0.3535533905932737*fields_no_J_nodal[7]-0.3535533905932737*fields_no_J_nodal[6]+0.3535533905932737*fields_no_J_nodal[5]-0.3535533905932737*fields_no_J_nodal[4]+0.3535533905932737*fields_no_J_nodal[3]-0.3535533905932737*fields_no_J_nodal[2]+0.3535533905932737*fields_no_J_nodal[1]-0.3535533905932737*fields_no_J_nodal[0]; 
  f_comp_no_J[4] = 0.3535533905932737*fields_no_J_nodal[7]+0.3535533905932737*fields_no_J_nodal[6]-0.3535533905932737*fields_no_J_nodal[5]-0.3535533905932737*fields_no_J_nodal[4]-0.3535533905932737*fields_no_J_nodal[3]-0.3535533905932737*fields_no_J_nodal[2]+0.3535533905932737*fields_no_J_nodal[1]+0.3535533905932737*fields_no_J_nodal[0]; 
  f_comp_no_J[5] = 0.3535533905932737*fields_no_J_nodal[7]-0.3535533905932737*fields_no_J_nodal[6]+0.3535533905932737*fields_no_J_nodal[5]-0.3535533905932737*fields_no_J_nodal[4]-0.3535533905932737*fields_no_J_nodal[3]+0.3535533905932737*fields_no_J_nodal[2]-0.3535533905932737*fields_no_J_nodal[1]+0.3535533905932737*fields_no_J_nodal[0]; 
  f_comp_no_J[6] = 0.3535533905932737*fields_no_J_nodal[7]-0.3535533905932737*fields_no_J_nodal[6]-0.3535533905932737*fields_no_J_nodal[5]+0.3535533905932737*fields_no_J_nodal[4]+0.3535533905932737*fields_no_J_nodal[3]-0.3535533905932737*fields_no_J_nodal[2]-0.3535533905932737*fields_no_J_nodal[1]+0.3535533905932737*fields_no_J_nodal[0]; 
  f_comp_no_J[7] = 0.3535533905932737*fields_no_J_nodal[7]-0.3535533905932737*fields_no_J_nodal[6]-0.3535533905932737*fields_no_J_nodal[5]+0.3535533905932737*fields_no_J_nodal[4]-0.3535533905932737*fields_no_J_nodal[3]+0.3535533905932737*fields_no_J_nodal[2]+0.3535533905932737*fields_no_J_nodal[1]-0.3535533905932737*fields_no_J_nodal[0]; 
  }
} 
