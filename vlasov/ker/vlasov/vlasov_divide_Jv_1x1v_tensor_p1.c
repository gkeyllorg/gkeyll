#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_divide_Jv_1x1v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J) 
{ 
  double f_no_J_nodal[6] = {0.0};

  f_no_J_nodal[0] = 1.0/jacob_vel_gauss[0]*(0.5*Jf[3]-0.5*Jf[2]-0.5*Jf[1]+0.5*Jf[0]);
  f_no_J_nodal[1] = 1.0/jacob_vel_gauss[1]*(-(0.5*Jf[3])+0.5*Jf[2]-0.5*Jf[1]+0.5*Jf[0]);
  f_no_J_nodal[2] = 1.0/jacob_vel_gauss[0]*(-(0.5*Jf[3])-0.5*Jf[2]+0.5*Jf[1]+0.5*Jf[0]);
  f_no_J_nodal[3] = 1.0/jacob_vel_gauss[1]*(0.5*Jf[3]+0.5*Jf[2]+0.5*Jf[1]+0.5*Jf[0]);

  f_no_J[0] = 0.5*f_no_J_nodal[3]+0.5*f_no_J_nodal[2]+0.5*f_no_J_nodal[1]+0.5*f_no_J_nodal[0]; 
  f_no_J[1] = 0.5*f_no_J_nodal[3]+0.5*f_no_J_nodal[2]-0.5*f_no_J_nodal[1]-0.5*f_no_J_nodal[0]; 
  f_no_J[2] = 0.5*f_no_J_nodal[3]-0.5*f_no_J_nodal[2]+0.5*f_no_J_nodal[1]-0.5*f_no_J_nodal[0]; 
  f_no_J[3] = 0.5*f_no_J_nodal[3]-0.5*f_no_J_nodal[2]-0.5*f_no_J_nodal[1]+0.5*f_no_J_nodal[0]; 
} 
