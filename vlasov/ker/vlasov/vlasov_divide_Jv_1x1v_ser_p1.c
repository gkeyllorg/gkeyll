#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_divide_Jv_1x1v_ser_p1(const double *jacob_vel_gauss, const double *Jf, double* GKYL_RESTRICT f_no_J) 
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
} 
