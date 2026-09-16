#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_n_set_1x1v_tensor_p3(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i) 
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
  gkyl_mat_set(&A0,0,0,0.7071067811865475*M0[0]); 
  gkyl_mat_set(&A0,0,1,0.7071067811865475*M0[1]); 
  gkyl_mat_set(&A0,0,2,0.7071067811865475*M0[2]); 
  gkyl_mat_set(&A0,0,3,0.7071067811865475*M0[3]); 
  gkyl_mat_set(&A0,1,0,0.7071067811865475*M0[1]); 
  gkyl_mat_set(&A0,1,1,0.6324555320336759*M0[2]+0.7071067811865475*M0[0]); 
  gkyl_mat_set(&A0,1,2,0.6210590034081186*M0[3]+0.6324555320336759*M0[1]); 
  gkyl_mat_set(&A0,1,3,0.6210590034081186*M0[2]); 
  gkyl_mat_set(&A0,2,0,0.7071067811865475*M0[2]); 
  gkyl_mat_set(&A0,2,1,0.6210590034081186*M0[3]+0.6324555320336759*M0[1]); 
  gkyl_mat_set(&A0,2,2,0.45175395145262565*M0[2]+0.7071067811865475*M0[0]); 
  gkyl_mat_set(&A0,2,3,0.421637021355784*M0[3]+0.6210590034081186*M0[1]); 
  gkyl_mat_set(&A0,3,0,0.7071067811865475*M0[3]); 
  gkyl_mat_set(&A0,3,1,0.6210590034081186*M0[2]); 
  gkyl_mat_set(&A0,3,2,0.421637021355784*M0[3]+0.6210590034081186*M0[1]); 
  gkyl_mat_set(&A0,3,3,0.421637021355784*M0[2]+0.7071067811865475*M0[0]); 
 
} 
