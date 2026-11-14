#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void em_calc_num_ExB_3x_tensor_p2(const double *em, double* GKYL_RESTRICT out) 
{ 
  // em:  Input electromagnetic fields. 
  // out: Output (E x B)_i (E x B velocity numerator) and B_i^2 (for denominator). 
 
  double *num_ExB_x = &out[0]; 
  double *num_ExB_y = &out[27]; 
  double *num_ExB_z = &out[54]; 
  double *B_x_sq  = &out[81]; 
  double *B_y_sq  = &out[108]; 
  double *B_z_sq  = &out[135]; 
 
  const double *E_x = &em[0]; 
  const double *E_y = &em[27]; 
  const double *E_z = &em[54]; 
  const double *B_x = &em[81]; 
  const double *B_y = &em[108]; 
  const double *B_z = &em[135]; 
 
  // Calculate E_i B_j. 
  double E_x_B_y[27] = {0.0}; 
  binop_mul_3d_tensor_p2(E_x, B_y, E_x_B_y); 
 
  double E_x_B_z[27] = {0.0}; 
  binop_mul_3d_tensor_p2(E_x, B_z, E_x_B_z); 
 
  double E_y_B_x[27] = {0.0}; 
  binop_mul_3d_tensor_p2(E_y, B_x, E_y_B_x); 
 
  double E_y_B_z[27] = {0.0}; 
  binop_mul_3d_tensor_p2(E_y, B_z, E_y_B_z); 
 
  double E_z_B_x[27] = {0.0}; 
  binop_mul_3d_tensor_p2(E_z, B_x, E_z_B_x); 
 
  double E_z_B_y[27] = {0.0}; 
  binop_mul_3d_tensor_p2(E_z, B_y, E_z_B_y); 
 
  num_ExB_x[0] = E_y_B_z[0] - E_z_B_y[0]; 
  num_ExB_y[0] = E_z_B_x[0] - E_x_B_z[0]; 
  num_ExB_z[0] = E_x_B_y[0] - E_y_B_x[0]; 
  num_ExB_x[1] = E_y_B_z[1] - E_z_B_y[1]; 
  num_ExB_y[1] = E_z_B_x[1] - E_x_B_z[1]; 
  num_ExB_z[1] = E_x_B_y[1] - E_y_B_x[1]; 
  num_ExB_x[2] = E_y_B_z[2] - E_z_B_y[2]; 
  num_ExB_y[2] = E_z_B_x[2] - E_x_B_z[2]; 
  num_ExB_z[2] = E_x_B_y[2] - E_y_B_x[2]; 
  num_ExB_x[3] = E_y_B_z[3] - E_z_B_y[3]; 
  num_ExB_y[3] = E_z_B_x[3] - E_x_B_z[3]; 
  num_ExB_z[3] = E_x_B_y[3] - E_y_B_x[3]; 
  num_ExB_x[4] = E_y_B_z[4] - E_z_B_y[4]; 
  num_ExB_y[4] = E_z_B_x[4] - E_x_B_z[4]; 
  num_ExB_z[4] = E_x_B_y[4] - E_y_B_x[4]; 
  num_ExB_x[5] = E_y_B_z[5] - E_z_B_y[5]; 
  num_ExB_y[5] = E_z_B_x[5] - E_x_B_z[5]; 
  num_ExB_z[5] = E_x_B_y[5] - E_y_B_x[5]; 
  num_ExB_x[6] = E_y_B_z[6] - E_z_B_y[6]; 
  num_ExB_y[6] = E_z_B_x[6] - E_x_B_z[6]; 
  num_ExB_z[6] = E_x_B_y[6] - E_y_B_x[6]; 
  num_ExB_x[7] = E_y_B_z[7] - E_z_B_y[7]; 
  num_ExB_y[7] = E_z_B_x[7] - E_x_B_z[7]; 
  num_ExB_z[7] = E_x_B_y[7] - E_y_B_x[7]; 
  num_ExB_x[8] = E_y_B_z[8] - E_z_B_y[8]; 
  num_ExB_y[8] = E_z_B_x[8] - E_x_B_z[8]; 
  num_ExB_z[8] = E_x_B_y[8] - E_y_B_x[8]; 
  num_ExB_x[9] = E_y_B_z[9] - E_z_B_y[9]; 
  num_ExB_y[9] = E_z_B_x[9] - E_x_B_z[9]; 
  num_ExB_z[9] = E_x_B_y[9] - E_y_B_x[9]; 
  num_ExB_x[10] = E_y_B_z[10] - E_z_B_y[10]; 
  num_ExB_y[10] = E_z_B_x[10] - E_x_B_z[10]; 
  num_ExB_z[10] = E_x_B_y[10] - E_y_B_x[10]; 
  num_ExB_x[11] = E_y_B_z[11] - E_z_B_y[11]; 
  num_ExB_y[11] = E_z_B_x[11] - E_x_B_z[11]; 
  num_ExB_z[11] = E_x_B_y[11] - E_y_B_x[11]; 
  num_ExB_x[12] = E_y_B_z[12] - E_z_B_y[12]; 
  num_ExB_y[12] = E_z_B_x[12] - E_x_B_z[12]; 
  num_ExB_z[12] = E_x_B_y[12] - E_y_B_x[12]; 
  num_ExB_x[13] = E_y_B_z[13] - E_z_B_y[13]; 
  num_ExB_y[13] = E_z_B_x[13] - E_x_B_z[13]; 
  num_ExB_z[13] = E_x_B_y[13] - E_y_B_x[13]; 
  num_ExB_x[14] = E_y_B_z[14] - E_z_B_y[14]; 
  num_ExB_y[14] = E_z_B_x[14] - E_x_B_z[14]; 
  num_ExB_z[14] = E_x_B_y[14] - E_y_B_x[14]; 
  num_ExB_x[15] = E_y_B_z[15] - E_z_B_y[15]; 
  num_ExB_y[15] = E_z_B_x[15] - E_x_B_z[15]; 
  num_ExB_z[15] = E_x_B_y[15] - E_y_B_x[15]; 
  num_ExB_x[16] = E_y_B_z[16] - E_z_B_y[16]; 
  num_ExB_y[16] = E_z_B_x[16] - E_x_B_z[16]; 
  num_ExB_z[16] = E_x_B_y[16] - E_y_B_x[16]; 
  num_ExB_x[17] = E_y_B_z[17] - E_z_B_y[17]; 
  num_ExB_y[17] = E_z_B_x[17] - E_x_B_z[17]; 
  num_ExB_z[17] = E_x_B_y[17] - E_y_B_x[17]; 
  num_ExB_x[18] = E_y_B_z[18] - E_z_B_y[18]; 
  num_ExB_y[18] = E_z_B_x[18] - E_x_B_z[18]; 
  num_ExB_z[18] = E_x_B_y[18] - E_y_B_x[18]; 
  num_ExB_x[19] = E_y_B_z[19] - E_z_B_y[19]; 
  num_ExB_y[19] = E_z_B_x[19] - E_x_B_z[19]; 
  num_ExB_z[19] = E_x_B_y[19] - E_y_B_x[19]; 
  num_ExB_x[20] = E_y_B_z[20] - E_z_B_y[20]; 
  num_ExB_y[20] = E_z_B_x[20] - E_x_B_z[20]; 
  num_ExB_z[20] = E_x_B_y[20] - E_y_B_x[20]; 
  num_ExB_x[21] = E_y_B_z[21] - E_z_B_y[21]; 
  num_ExB_y[21] = E_z_B_x[21] - E_x_B_z[21]; 
  num_ExB_z[21] = E_x_B_y[21] - E_y_B_x[21]; 
  num_ExB_x[22] = E_y_B_z[22] - E_z_B_y[22]; 
  num_ExB_y[22] = E_z_B_x[22] - E_x_B_z[22]; 
  num_ExB_z[22] = E_x_B_y[22] - E_y_B_x[22]; 
  num_ExB_x[23] = E_y_B_z[23] - E_z_B_y[23]; 
  num_ExB_y[23] = E_z_B_x[23] - E_x_B_z[23]; 
  num_ExB_z[23] = E_x_B_y[23] - E_y_B_x[23]; 
  num_ExB_x[24] = E_y_B_z[24] - E_z_B_y[24]; 
  num_ExB_y[24] = E_z_B_x[24] - E_x_B_z[24]; 
  num_ExB_z[24] = E_x_B_y[24] - E_y_B_x[24]; 
  num_ExB_x[25] = E_y_B_z[25] - E_z_B_y[25]; 
  num_ExB_y[25] = E_z_B_x[25] - E_x_B_z[25]; 
  num_ExB_z[25] = E_x_B_y[25] - E_y_B_x[25]; 
  num_ExB_x[26] = E_y_B_z[26] - E_z_B_y[26]; 
  num_ExB_y[26] = E_z_B_x[26] - E_x_B_z[26]; 
  num_ExB_z[26] = E_x_B_y[26] - E_y_B_x[26]; 
 
  // Calculate B_i^2. 
  binop_mul_3d_tensor_p2(B_x, B_x, B_x_sq); 
  binop_mul_3d_tensor_p2(B_y, B_y, B_y_sq); 
  binop_mul_3d_tensor_p2(B_z, B_z, B_z_sq); 
 
} 
 
