#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_divide_Jc_1x_tensor_p2(const double *J, const double *Jfields, double* GKYL_RESTRICT fields_no_J) 
{ 
  const double *Jf_comp;
  double *f_comp_no_J;
  for (int i=0; i<8; ++i) {
  Jf_comp = &Jfields[i*3];
  f_comp_no_J = &fields_no_J[i*3];
    double fields_no_J_nodal[3] = {0.0};
    
    fields_no_J_nodal[0] = 1.0/J[0]*(0.6324555320336759*Jf_comp[2]-0.9486832980505137*Jf_comp[1]+0.7071067811865475*Jf_comp[0]);
    fields_no_J_nodal[1] = 1.0/J[1]*(0.7071067811865475*Jf_comp[0]-0.7905694150420947*Jf_comp[2]);
    fields_no_J_nodal[2] = 1.0/J[2]*(0.6324555320336759*Jf_comp[2]+0.9486832980505137*Jf_comp[1]+0.7071067811865475*Jf_comp[0]);
  
  f_comp_no_J[0] = 0.39283710065919303*fields_no_J_nodal[2]+0.6285393610547091*fields_no_J_nodal[1]+0.39283710065919303*fields_no_J_nodal[0]; 
  f_comp_no_J[1] = 0.5270462766947298*fields_no_J_nodal[2]-0.5270462766947298*fields_no_J_nodal[0]; 
  f_comp_no_J[2] = 0.35136418446315326*fields_no_J_nodal[2]-0.7027283689263066*fields_no_J_nodal[1]+0.35136418446315326*fields_no_J_nodal[0]; 
  }
} 
