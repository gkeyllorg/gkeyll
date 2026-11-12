#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_rescale_Jv_1x1v_ser_p1(const double *jacob_vel_gauss, const double *f_no_J, double* GKYL_RESTRICT Jf) 
{ 
  // jacob_vel_gauss: Total velocity-space Jacobian at Gauss-Legendre quadrature points.
  //                  Piecewise constant in the cell for the C^0 linear map of the
  //                  Serendipity basis, so the rescaling is exact coefficient-by-coefficient.
  // f_no_J:          Input distribution function (f).
  // Jf:              Output distribution function (J*f).
 
  Jf[0] = jacob_vel_gauss[0]*f_no_J[0]; 
  Jf[1] = jacob_vel_gauss[0]*f_no_J[1]; 
  Jf[2] = jacob_vel_gauss[0]*f_no_J[2]; 
  Jf[3] = jacob_vel_gauss[0]*f_no_J[3]; 
} 
