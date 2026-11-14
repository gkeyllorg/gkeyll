#include <gkyl_mat.h> 
#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH void em_copy_ExB_3x_tensor_p2(int count, struct gkyl_nmat *x, const double *em, int* cell_avg_magB2, 
  double* GKYL_RESTRICT ExB, double* GKYL_RESTRICT ExB_surf) 
{ 
  // count:          Integer to indicate which matrix being fetched. 
  // x:              Input solution vector. 
  // em:             Input electromagnetic fields. 
  // cell_avg_magB2: Input flag for cell average if 1/|B|^2 only used cell averages. 
  // ExB:            E x B velocity = E x B/|B|^2. 
 
  struct gkyl_mat x_ExBx = gkyl_nmat_get(x, count); 
  struct gkyl_mat x_ExBy = gkyl_nmat_get(x, count+1); 
  struct gkyl_mat x_ExBz = gkyl_nmat_get(x, count+2); 
  double *ExB_x = &ExB[0]; 
  double *ExB_y = &ExB[27]; 
  double *ExB_z = &ExB[54]; 
 
  ExB_x[0] = gkyl_mat_get(&x_ExBx,0,0); 
  ExB_y[0] = gkyl_mat_get(&x_ExBy,0,0); 
  ExB_z[0] = gkyl_mat_get(&x_ExBz,0,0); 
  ExB_x[1] = gkyl_mat_get(&x_ExBx,1,0); 
  ExB_y[1] = gkyl_mat_get(&x_ExBy,1,0); 
  ExB_z[1] = gkyl_mat_get(&x_ExBz,1,0); 
  ExB_x[2] = gkyl_mat_get(&x_ExBx,2,0); 
  ExB_y[2] = gkyl_mat_get(&x_ExBy,2,0); 
  ExB_z[2] = gkyl_mat_get(&x_ExBz,2,0); 
  ExB_x[3] = gkyl_mat_get(&x_ExBx,3,0); 
  ExB_y[3] = gkyl_mat_get(&x_ExBy,3,0); 
  ExB_z[3] = gkyl_mat_get(&x_ExBz,3,0); 
  ExB_x[4] = gkyl_mat_get(&x_ExBx,4,0); 
  ExB_y[4] = gkyl_mat_get(&x_ExBy,4,0); 
  ExB_z[4] = gkyl_mat_get(&x_ExBz,4,0); 
  ExB_x[5] = gkyl_mat_get(&x_ExBx,5,0); 
  ExB_y[5] = gkyl_mat_get(&x_ExBy,5,0); 
  ExB_z[5] = gkyl_mat_get(&x_ExBz,5,0); 
  ExB_x[6] = gkyl_mat_get(&x_ExBx,6,0); 
  ExB_y[6] = gkyl_mat_get(&x_ExBy,6,0); 
  ExB_z[6] = gkyl_mat_get(&x_ExBz,6,0); 
  ExB_x[7] = gkyl_mat_get(&x_ExBx,7,0); 
  ExB_y[7] = gkyl_mat_get(&x_ExBy,7,0); 
  ExB_z[7] = gkyl_mat_get(&x_ExBz,7,0); 
  ExB_x[8] = gkyl_mat_get(&x_ExBx,8,0); 
  ExB_y[8] = gkyl_mat_get(&x_ExBy,8,0); 
  ExB_z[8] = gkyl_mat_get(&x_ExBz,8,0); 
  ExB_x[9] = gkyl_mat_get(&x_ExBx,9,0); 
  ExB_y[9] = gkyl_mat_get(&x_ExBy,9,0); 
  ExB_z[9] = gkyl_mat_get(&x_ExBz,9,0); 
  ExB_x[10] = gkyl_mat_get(&x_ExBx,10,0); 
  ExB_y[10] = gkyl_mat_get(&x_ExBy,10,0); 
  ExB_z[10] = gkyl_mat_get(&x_ExBz,10,0); 
  ExB_x[11] = gkyl_mat_get(&x_ExBx,11,0); 
  ExB_y[11] = gkyl_mat_get(&x_ExBy,11,0); 
  ExB_z[11] = gkyl_mat_get(&x_ExBz,11,0); 
  ExB_x[12] = gkyl_mat_get(&x_ExBx,12,0); 
  ExB_y[12] = gkyl_mat_get(&x_ExBy,12,0); 
  ExB_z[12] = gkyl_mat_get(&x_ExBz,12,0); 
  ExB_x[13] = gkyl_mat_get(&x_ExBx,13,0); 
  ExB_y[13] = gkyl_mat_get(&x_ExBy,13,0); 
  ExB_z[13] = gkyl_mat_get(&x_ExBz,13,0); 
  ExB_x[14] = gkyl_mat_get(&x_ExBx,14,0); 
  ExB_y[14] = gkyl_mat_get(&x_ExBy,14,0); 
  ExB_z[14] = gkyl_mat_get(&x_ExBz,14,0); 
  ExB_x[15] = gkyl_mat_get(&x_ExBx,15,0); 
  ExB_y[15] = gkyl_mat_get(&x_ExBy,15,0); 
  ExB_z[15] = gkyl_mat_get(&x_ExBz,15,0); 
  ExB_x[16] = gkyl_mat_get(&x_ExBx,16,0); 
  ExB_y[16] = gkyl_mat_get(&x_ExBy,16,0); 
  ExB_z[16] = gkyl_mat_get(&x_ExBz,16,0); 
  ExB_x[17] = gkyl_mat_get(&x_ExBx,17,0); 
  ExB_y[17] = gkyl_mat_get(&x_ExBy,17,0); 
  ExB_z[17] = gkyl_mat_get(&x_ExBz,17,0); 
  ExB_x[18] = gkyl_mat_get(&x_ExBx,18,0); 
  ExB_y[18] = gkyl_mat_get(&x_ExBy,18,0); 
  ExB_z[18] = gkyl_mat_get(&x_ExBz,18,0); 
  ExB_x[19] = gkyl_mat_get(&x_ExBx,19,0); 
  ExB_y[19] = gkyl_mat_get(&x_ExBy,19,0); 
  ExB_z[19] = gkyl_mat_get(&x_ExBz,19,0); 
  ExB_x[20] = gkyl_mat_get(&x_ExBx,20,0); 
  ExB_y[20] = gkyl_mat_get(&x_ExBy,20,0); 
  ExB_z[20] = gkyl_mat_get(&x_ExBz,20,0); 
  ExB_x[21] = gkyl_mat_get(&x_ExBx,21,0); 
  ExB_y[21] = gkyl_mat_get(&x_ExBy,21,0); 
  ExB_z[21] = gkyl_mat_get(&x_ExBz,21,0); 
  ExB_x[22] = gkyl_mat_get(&x_ExBx,22,0); 
  ExB_y[22] = gkyl_mat_get(&x_ExBy,22,0); 
  ExB_z[22] = gkyl_mat_get(&x_ExBz,22,0); 
  ExB_x[23] = gkyl_mat_get(&x_ExBx,23,0); 
  ExB_y[23] = gkyl_mat_get(&x_ExBy,23,0); 
  ExB_z[23] = gkyl_mat_get(&x_ExBz,23,0); 
  ExB_x[24] = gkyl_mat_get(&x_ExBx,24,0); 
  ExB_y[24] = gkyl_mat_get(&x_ExBy,24,0); 
  ExB_z[24] = gkyl_mat_get(&x_ExBz,24,0); 
  ExB_x[25] = gkyl_mat_get(&x_ExBx,25,0); 
  ExB_y[25] = gkyl_mat_get(&x_ExBy,25,0); 
  ExB_z[25] = gkyl_mat_get(&x_ExBz,25,0); 
  ExB_x[26] = gkyl_mat_get(&x_ExBx,26,0); 
  ExB_y[26] = gkyl_mat_get(&x_ExBy,26,0); 
  ExB_z[26] = gkyl_mat_get(&x_ExBz,26,0); 
} 
 
