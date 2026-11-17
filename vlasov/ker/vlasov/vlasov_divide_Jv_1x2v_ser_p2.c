#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_divide_Jv_1x2v_ser_p2(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J) 
{ 
  // jacob_vel_gauss: Total velocity-space Jacobian at Gauss-Legendre quadrature points.
  //                  Piecewise constant in the cell for the C^0 linear map of the
  //                  Serendipity basis, so the division is exact coefficient-by-coefficient.
  // Jf:              Input distribution function (J*f).
  // f_no_J:          Output distribution function (f).
 
  const double jacob_vel_inv = 1.0/jacob_vel_gauss[0]; 
  f_no_J[0] = jacob_vel_inv*Jf[0]; 
  f_no_J[1] = jacob_vel_inv*Jf[1]; 
  f_no_J[2] = jacob_vel_inv*Jf[2]; 
  f_no_J[3] = jacob_vel_inv*Jf[3]; 
  f_no_J[4] = jacob_vel_inv*Jf[4]; 
  f_no_J[5] = jacob_vel_inv*Jf[5]; 
  f_no_J[6] = jacob_vel_inv*Jf[6]; 
  f_no_J[7] = jacob_vel_inv*Jf[7]; 
  f_no_J[8] = jacob_vel_inv*Jf[8]; 
  f_no_J[9] = jacob_vel_inv*Jf[9]; 
  f_no_J[10] = jacob_vel_inv*Jf[10]; 
  f_no_J[11] = jacob_vel_inv*Jf[11]; 
  f_no_J[12] = jacob_vel_inv*Jf[12]; 
  f_no_J[13] = jacob_vel_inv*Jf[13]; 
  f_no_J[14] = jacob_vel_inv*Jf[14]; 
  f_no_J[15] = jacob_vel_inv*Jf[15]; 
  f_no_J[16] = jacob_vel_inv*Jf[16]; 
  f_no_J[17] = jacob_vel_inv*Jf[17]; 
  f_no_J[18] = jacob_vel_inv*Jf[18]; 
  f_no_J[19] = jacob_vel_inv*Jf[19]; 
} 
