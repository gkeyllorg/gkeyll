#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_divide_Jv_1x1v_tensor_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J) 
{ 
  // Velocity-space Jacobian at the velocity Gauss-Legendre nodes; each configuration-space mode
  // is divided by it, evaluated at the nodes and projected back one velocity direction at a time.
  double jinv[3]; 
  jinv[0] = 1.0/jacob_vel_gauss[0]; 
  jinv[1] = 1.0/jacob_vel_gauss[1]; 
  jinv[2] = 1.0/jacob_vel_gauss[2]; 
  double ta[3]; 

  // Configuration-space mode 0. 
  ta[0] = 0.6324555320336759*Jf[4]-0.9486832980505137*Jf[2]+0.7071067811865475*Jf[0]; 
  ta[1] = 0.7071067811865475*Jf[0]-0.7905694150420947*Jf[4]; 
  ta[2] = 0.6324555320336759*Jf[4]+0.9486832980505137*Jf[2]+0.7071067811865475*Jf[0]; 
  f_no_J[0] = 0.39283710065919303*jinv[2]*ta[2]+0.6285393610547091*jinv[1]*ta[1]+0.39283710065919303*jinv[0]*ta[0]; 
  f_no_J[2] = 0.5270462766947298*jinv[2]*ta[2]-0.5270462766947298*jinv[0]*ta[0]; 
  f_no_J[4] = 0.35136418446315326*jinv[2]*ta[2]-0.7027283689263066*jinv[1]*ta[1]+0.35136418446315326*jinv[0]*ta[0]; 

  // Configuration-space mode 1. 
  ta[0] = 0.6324555320336759*Jf[5]-0.9486832980505137*Jf[3]+0.7071067811865475*Jf[1]; 
  ta[1] = 0.7071067811865475*Jf[1]-0.7905694150420947*Jf[5]; 
  ta[2] = 0.6324555320336759*Jf[5]+0.9486832980505137*Jf[3]+0.7071067811865475*Jf[1]; 
  f_no_J[1] = 0.39283710065919303*jinv[2]*ta[2]+0.6285393610547091*jinv[1]*ta[1]+0.39283710065919303*jinv[0]*ta[0]; 
  f_no_J[3] = 0.5270462766947298*jinv[2]*ta[2]-0.5270462766947298*jinv[0]*ta[0]; 
  f_no_J[5] = 0.35136418446315326*jinv[2]*ta[2]-0.7027283689263066*jinv[1]*ta[1]+0.35136418446315326*jinv[0]*ta[0]; 
} 
