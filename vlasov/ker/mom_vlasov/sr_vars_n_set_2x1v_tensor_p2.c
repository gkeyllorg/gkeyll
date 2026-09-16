#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_n_set_2x1v_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A:     preallocated LHS matrix. 
  // rhs:   preallocated RHS vector. 
  // M0:    Lab frame density = Gamma*n.
  // M1i:   Lab frame flux = Gamma*n*V_drift_i.

  struct gkyl_mat A0 = gkyl_nmat_get(A, count+0); 
  struct gkyl_mat rhs0 = gkyl_nmat_get(rhs, count+0); 
  gkyl_mat_clear(&A0, 0.0); gkyl_mat_clear(&rhs0, 0.0); 
  const double *M10 = &M1i[0]; 
  gkyl_mat_set(&rhs0,0,0,M10[0]); 
  gkyl_mat_set(&rhs0,1,0,M10[1]); 
  gkyl_mat_set(&rhs0,2,0,M10[2]); 
  gkyl_mat_set(&rhs0,3,0,M10[3]); 
  gkyl_mat_set(&rhs0,4,0,M10[4]); 
  gkyl_mat_set(&rhs0,5,0,M10[5]); 
  gkyl_mat_set(&rhs0,6,0,M10[6]); 
  gkyl_mat_set(&rhs0,7,0,M10[7]); 
  gkyl_mat_set(&rhs0,8,0,M10[8]); 
  gkyl_mat_set(&A0,0,0,0.5*M0[0]); 
  gkyl_mat_set(&A0,0,1,0.5*M0[1]); 
  gkyl_mat_set(&A0,0,2,0.5*M0[2]); 
  gkyl_mat_set(&A0,0,3,0.5*M0[3]); 
  gkyl_mat_set(&A0,0,4,0.5*M0[4]); 
  gkyl_mat_set(&A0,0,5,0.5*M0[5]); 
  gkyl_mat_set(&A0,0,6,0.5*M0[6]); 
  gkyl_mat_set(&A0,0,7,0.5*M0[7]); 
  gkyl_mat_set(&A0,0,8,0.5*M0[8]); 
  gkyl_mat_set(&A0,1,0,0.5*M0[1]); 
  gkyl_mat_set(&A0,1,1,0.4472135954999579*M0[4]+0.5*M0[0]); 
  gkyl_mat_set(&A0,1,2,0.5*M0[3]); 
  gkyl_mat_set(&A0,1,3,0.44721359549995804*M0[6]+0.5*M0[2]); 
  gkyl_mat_set(&A0,1,4,0.4472135954999579*M0[1]); 
  gkyl_mat_set(&A0,1,5,0.5000000000000001*M0[7]); 
  gkyl_mat_set(&A0,1,6,0.44721359549995804*M0[3]); 
  gkyl_mat_set(&A0,1,7,0.44721359549995804*M0[8]+0.5000000000000001*M0[5]); 
  gkyl_mat_set(&A0,1,8,0.44721359549995804*M0[7]); 
  gkyl_mat_set(&A0,2,0,0.5*M0[2]); 
  gkyl_mat_set(&A0,2,1,0.5*M0[3]); 
  gkyl_mat_set(&A0,2,2,0.4472135954999579*M0[5]+0.5*M0[0]); 
  gkyl_mat_set(&A0,2,3,0.44721359549995804*M0[7]+0.5*M0[1]); 
  gkyl_mat_set(&A0,2,4,0.5000000000000001*M0[6]); 
  gkyl_mat_set(&A0,2,5,0.4472135954999579*M0[2]); 
  gkyl_mat_set(&A0,2,6,0.44721359549995804*M0[8]+0.5000000000000001*M0[4]); 
  gkyl_mat_set(&A0,2,7,0.44721359549995804*M0[3]); 
  gkyl_mat_set(&A0,2,8,0.44721359549995804*M0[6]); 
  gkyl_mat_set(&A0,3,0,0.5*M0[3]); 
  gkyl_mat_set(&A0,3,1,0.44721359549995804*M0[6]+0.5*M0[2]); 
  gkyl_mat_set(&A0,3,2,0.44721359549995804*M0[7]+0.5*M0[1]); 
  gkyl_mat_set(&A0,3,3,0.4*M0[8]+0.4472135954999579*M0[5]+0.4472135954999579*M0[4]+0.5*M0[0]); 
  gkyl_mat_set(&A0,3,4,0.4472135954999579*M0[3]); 
  gkyl_mat_set(&A0,3,5,0.4472135954999579*M0[3]); 
  gkyl_mat_set(&A0,3,6,0.4*M0[7]+0.44721359549995804*M0[1]); 
  gkyl_mat_set(&A0,3,7,0.4*M0[6]+0.44721359549995804*M0[2]); 
  gkyl_mat_set(&A0,3,8,0.4*M0[3]); 
  gkyl_mat_set(&A0,4,0,0.5*M0[4]); 
  gkyl_mat_set(&A0,4,1,0.4472135954999579*M0[1]); 
  gkyl_mat_set(&A0,4,2,0.5000000000000001*M0[6]); 
  gkyl_mat_set(&A0,4,3,0.4472135954999579*M0[3]); 
  gkyl_mat_set(&A0,4,4,0.31943828249996997*M0[4]+0.5*M0[0]); 
  gkyl_mat_set(&A0,4,5,0.5*M0[8]); 
  gkyl_mat_set(&A0,4,6,0.31943828249996997*M0[6]+0.5000000000000001*M0[2]); 
  gkyl_mat_set(&A0,4,7,0.4472135954999579*M0[7]); 
  gkyl_mat_set(&A0,4,8,0.31943828249996997*M0[8]+0.5*M0[5]); 
  gkyl_mat_set(&A0,5,0,0.5*M0[5]); 
  gkyl_mat_set(&A0,5,1,0.5000000000000001*M0[7]); 
  gkyl_mat_set(&A0,5,2,0.4472135954999579*M0[2]); 
  gkyl_mat_set(&A0,5,3,0.4472135954999579*M0[3]); 
  gkyl_mat_set(&A0,5,4,0.5*M0[8]); 
  gkyl_mat_set(&A0,5,5,0.31943828249996997*M0[5]+0.5*M0[0]); 
  gkyl_mat_set(&A0,5,6,0.4472135954999579*M0[6]); 
  gkyl_mat_set(&A0,5,7,0.31943828249996997*M0[7]+0.5000000000000001*M0[1]); 
  gkyl_mat_set(&A0,5,8,0.31943828249996997*M0[8]+0.5*M0[4]); 
  gkyl_mat_set(&A0,6,0,0.5*M0[6]); 
  gkyl_mat_set(&A0,6,1,0.44721359549995804*M0[3]); 
  gkyl_mat_set(&A0,6,2,0.44721359549995804*M0[8]+0.5000000000000001*M0[4]); 
  gkyl_mat_set(&A0,6,3,0.4*M0[7]+0.44721359549995804*M0[1]); 
  gkyl_mat_set(&A0,6,4,0.31943828249996997*M0[6]+0.5000000000000001*M0[2]); 
  gkyl_mat_set(&A0,6,5,0.4472135954999579*M0[6]); 
  gkyl_mat_set(&A0,6,6,0.2857142857142857*M0[8]+0.4472135954999579*M0[5]+0.31943828249996997*M0[4]+0.5*M0[0]); 
  gkyl_mat_set(&A0,6,7,0.4*M0[3]); 
  gkyl_mat_set(&A0,6,8,0.2857142857142857*M0[6]+0.44721359549995804*M0[2]); 
  gkyl_mat_set(&A0,7,0,0.5*M0[7]); 
  gkyl_mat_set(&A0,7,1,0.44721359549995804*M0[8]+0.5000000000000001*M0[5]); 
  gkyl_mat_set(&A0,7,2,0.44721359549995804*M0[3]); 
  gkyl_mat_set(&A0,7,3,0.4*M0[6]+0.44721359549995804*M0[2]); 
  gkyl_mat_set(&A0,7,4,0.4472135954999579*M0[7]); 
  gkyl_mat_set(&A0,7,5,0.31943828249996997*M0[7]+0.5000000000000001*M0[1]); 
  gkyl_mat_set(&A0,7,6,0.4*M0[3]); 
  gkyl_mat_set(&A0,7,7,0.2857142857142857*M0[8]+0.31943828249996997*M0[5]+0.4472135954999579*M0[4]+0.5*M0[0]); 
  gkyl_mat_set(&A0,7,8,0.2857142857142857*M0[7]+0.44721359549995804*M0[1]); 
  gkyl_mat_set(&A0,8,0,0.5*M0[8]); 
  gkyl_mat_set(&A0,8,1,0.44721359549995804*M0[7]); 
  gkyl_mat_set(&A0,8,2,0.44721359549995804*M0[6]); 
  gkyl_mat_set(&A0,8,3,0.4*M0[3]); 
  gkyl_mat_set(&A0,8,4,0.31943828249996997*M0[8]+0.5*M0[5]); 
  gkyl_mat_set(&A0,8,5,0.31943828249996997*M0[8]+0.5*M0[4]); 
  gkyl_mat_set(&A0,8,6,0.2857142857142857*M0[6]+0.44721359549995804*M0[2]); 
  gkyl_mat_set(&A0,8,7,0.2857142857142857*M0[7]+0.44721359549995804*M0[1]); 
  gkyl_mat_set(&A0,8,8,0.20408163265306123*M0[8]+0.31943828249996997*M0[5]+0.31943828249996997*M0[4]+0.5*M0[0]); 
 
} 
