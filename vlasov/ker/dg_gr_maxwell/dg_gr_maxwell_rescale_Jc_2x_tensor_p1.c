#include <gkyl_dg_gr_maxwell_kernels.h> 
GKYL_CU_DH void dg_gr_maxwell_rescale_Jc_2x_tensor_p1(const double *J, const double *fields_no_J, double* GKYL_RESTRICT Jfields) 
{ 
  const double *f_comp_no_J;
  double *Jf_comp;
  for (int i=0; i<8; ++i) {
    f_comp_no_J = &fields_no_J[i*4];
    Jf_comp = &Jfields[i*4];
    double Jfields_nodal[4] = {0.0};
    
    Jfields_nodal[0] = J[0]*(0.5*f_comp_no_J[3]-0.5*f_comp_no_J[2]-0.5*f_comp_no_J[1]+0.5*f_comp_no_J[0]);
    Jfields_nodal[1] = J[1]*(-(0.5*f_comp_no_J[3])+0.5*f_comp_no_J[2]-0.5*f_comp_no_J[1]+0.5*f_comp_no_J[0]);
    Jfields_nodal[2] = J[2]*(-(0.5*f_comp_no_J[3])-0.5*f_comp_no_J[2]+0.5*f_comp_no_J[1]+0.5*f_comp_no_J[0]);
    Jfields_nodal[3] = J[3]*(0.5*f_comp_no_J[3]+0.5*f_comp_no_J[2]+0.5*f_comp_no_J[1]+0.5*f_comp_no_J[0]);
    
  Jf_comp[0] = 0.5*Jfields_nodal[3]+0.5*Jfields_nodal[2]+0.5*Jfields_nodal[1]+0.5*Jfields_nodal[0]; 
  Jf_comp[1] = 0.5*Jfields_nodal[3]+0.5*Jfields_nodal[2]-0.5*Jfields_nodal[1]-0.5*Jfields_nodal[0]; 
  Jf_comp[2] = 0.5*Jfields_nodal[3]-0.5*Jfields_nodal[2]+0.5*Jfields_nodal[1]-0.5*Jfields_nodal[0]; 
  Jf_comp[3] = 0.5*Jfields_nodal[3]-0.5*Jfields_nodal[2]-0.5*Jfields_nodal[1]+0.5*Jfields_nodal[0]; 
  }
} 
