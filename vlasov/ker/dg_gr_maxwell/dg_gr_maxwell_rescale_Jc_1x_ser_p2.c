#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_1x_ser_p2(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields) 
{ 
  const double *f_comp_no_J;
  double *Jf_comp;
  for (int i=0; i<8; ++i) {
    f_comp_no_J = &fields_no_J[i*3];
    Jf_comp = &Jfields[i*3];
    double Jfields_nodal[3] = {0.0};
    
    Jfields_nodal[0] = J[0]*(0.6324555320336759*f_comp_no_J[2]-0.9486832980505137*f_comp_no_J[1]+0.7071067811865475*f_comp_no_J[0]);
    Jfields_nodal[1] = J[1]*(0.7071067811865475*f_comp_no_J[0]-0.7905694150420947*f_comp_no_J[2]);
    Jfields_nodal[2] = J[2]*(0.6324555320336759*f_comp_no_J[2]+0.9486832980505137*f_comp_no_J[1]+0.7071067811865475*f_comp_no_J[0]);
    
  Jf_comp[0] = 0.39283710065919303*Jfields_nodal[2]+0.6285393610547091*Jfields_nodal[1]+0.39283710065919303*Jfields_nodal[0]; 
  Jf_comp[1] = 0.5270462766947298*Jfields_nodal[2]-0.5270462766947298*Jfields_nodal[0]; 
  Jf_comp[2] = 0.35136418446315326*Jfields_nodal[2]-0.7027283689263066*Jfields_nodal[1]+0.35136418446315326*Jfields_nodal[0]; 
  }
} 
