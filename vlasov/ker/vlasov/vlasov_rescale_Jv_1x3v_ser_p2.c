#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_rescale_Jv_1x3v_ser_p2(const double *jacob_vel_gauss, const double *f_no_J, double* GKYL_RESTRICT Jf) 
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
  Jf[4] = jacob_vel_gauss[0]*f_no_J[4]; 
  Jf[5] = jacob_vel_gauss[0]*f_no_J[5]; 
  Jf[6] = jacob_vel_gauss[0]*f_no_J[6]; 
  Jf[7] = jacob_vel_gauss[0]*f_no_J[7]; 
  Jf[8] = jacob_vel_gauss[0]*f_no_J[8]; 
  Jf[9] = jacob_vel_gauss[0]*f_no_J[9]; 
  Jf[10] = jacob_vel_gauss[0]*f_no_J[10]; 
  Jf[11] = jacob_vel_gauss[0]*f_no_J[11]; 
  Jf[12] = jacob_vel_gauss[0]*f_no_J[12]; 
  Jf[13] = jacob_vel_gauss[0]*f_no_J[13]; 
  Jf[14] = jacob_vel_gauss[0]*f_no_J[14]; 
  Jf[15] = jacob_vel_gauss[0]*f_no_J[15]; 
  Jf[16] = jacob_vel_gauss[0]*f_no_J[16]; 
  Jf[17] = jacob_vel_gauss[0]*f_no_J[17]; 
  Jf[18] = jacob_vel_gauss[0]*f_no_J[18]; 
  Jf[19] = jacob_vel_gauss[0]*f_no_J[19]; 
  Jf[20] = jacob_vel_gauss[0]*f_no_J[20]; 
  Jf[21] = jacob_vel_gauss[0]*f_no_J[21]; 
  Jf[22] = jacob_vel_gauss[0]*f_no_J[22]; 
  Jf[23] = jacob_vel_gauss[0]*f_no_J[23]; 
  Jf[24] = jacob_vel_gauss[0]*f_no_J[24]; 
  Jf[25] = jacob_vel_gauss[0]*f_no_J[25]; 
  Jf[26] = jacob_vel_gauss[0]*f_no_J[26]; 
  Jf[27] = jacob_vel_gauss[0]*f_no_J[27]; 
  Jf[28] = jacob_vel_gauss[0]*f_no_J[28]; 
  Jf[29] = jacob_vel_gauss[0]*f_no_J[29]; 
  Jf[30] = jacob_vel_gauss[0]*f_no_J[30]; 
  Jf[31] = jacob_vel_gauss[0]*f_no_J[31]; 
  Jf[32] = jacob_vel_gauss[0]*f_no_J[32]; 
  Jf[33] = jacob_vel_gauss[0]*f_no_J[33]; 
  Jf[34] = jacob_vel_gauss[0]*f_no_J[34]; 
  Jf[35] = jacob_vel_gauss[0]*f_no_J[35]; 
  Jf[36] = jacob_vel_gauss[0]*f_no_J[36]; 
  Jf[37] = jacob_vel_gauss[0]*f_no_J[37]; 
  Jf[38] = jacob_vel_gauss[0]*f_no_J[38]; 
  Jf[39] = jacob_vel_gauss[0]*f_no_J[39]; 
  Jf[40] = jacob_vel_gauss[0]*f_no_J[40]; 
  Jf[41] = jacob_vel_gauss[0]*f_no_J[41]; 
  Jf[42] = jacob_vel_gauss[0]*f_no_J[42]; 
  Jf[43] = jacob_vel_gauss[0]*f_no_J[43]; 
  Jf[44] = jacob_vel_gauss[0]*f_no_J[44]; 
  Jf[45] = jacob_vel_gauss[0]*f_no_J[45]; 
  Jf[46] = jacob_vel_gauss[0]*f_no_J[46]; 
  Jf[47] = jacob_vel_gauss[0]*f_no_J[47]; 
} 
