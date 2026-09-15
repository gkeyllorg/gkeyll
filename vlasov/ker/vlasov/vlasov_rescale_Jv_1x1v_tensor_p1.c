#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_rescale_Jv_1x1v_tensor_p1(const double *jacob_vel_gauss, const double *f_no_J, double* GKYL_RESTRICT Jf) 
{ 
  double Jf_nodal[6] = {0.0};

  Jf_nodal[0] = jacob_vel_gauss[0]*(0.5*f_no_J[3]-0.5*f_no_J[2]-0.5*f_no_J[1]+0.5*f_no_J[0]);
  Jf_nodal[1] = jacob_vel_gauss[1]*(-(0.5*f_no_J[3])+0.5*f_no_J[2]-0.5*f_no_J[1]+0.5*f_no_J[0]);
  Jf_nodal[2] = jacob_vel_gauss[0]*(-(0.5*f_no_J[3])-0.5*f_no_J[2]+0.5*f_no_J[1]+0.5*f_no_J[0]);
  Jf_nodal[3] = jacob_vel_gauss[1]*(0.5*f_no_J[3]+0.5*f_no_J[2]+0.5*f_no_J[1]+0.5*f_no_J[0]);

  Jf[0] = 0.5*Jf_nodal[3]+0.5*Jf_nodal[2]+0.5*Jf_nodal[1]+0.5*Jf_nodal[0]; 
  Jf[1] = 0.5*Jf_nodal[3]+0.5*Jf_nodal[2]-0.5*Jf_nodal[1]-0.5*Jf_nodal[0]; 
  Jf[2] = 0.5*Jf_nodal[3]-0.5*Jf_nodal[2]+0.5*Jf_nodal[1]-0.5*Jf_nodal[0]; 
  Jf[3] = 0.5*Jf_nodal[3]-0.5*Jf_nodal[2]-0.5*Jf_nodal[1]+0.5*Jf_nodal[0]; 
} 
