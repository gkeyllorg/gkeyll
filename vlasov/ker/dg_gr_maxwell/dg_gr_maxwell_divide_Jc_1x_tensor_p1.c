#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_tensor_p1(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J) 
{ 
  const double *Jf_comp;
  double *f_comp_no_J;
  for (int i=0; i<8; ++i) {
  Jf_comp = &Jfields[i*2];
  f_comp_no_J = &fields_no_J[i*2];
    double fields_no_J_nodal[2] = {0.0};
    
    fields_no_J_nodal[0] = 1.0/J[0]*(0.7071067811865475*Jf_comp[0]-0.7071067811865475*Jf_comp[1]);
    fields_no_J_nodal[1] = 1.0/J[1]*(0.7071067811865475*Jf_comp[1]+0.7071067811865475*Jf_comp[0]);
  
  f_comp_no_J[0] = 0.7071067811865475*fields_no_J_nodal[1]+0.7071067811865475*fields_no_J_nodal[0]; 
  f_comp_no_J[1] = 0.7071067811865475*fields_no_J_nodal[1]-0.7071067811865475*fields_no_J_nodal[0]; 
  }
} 
