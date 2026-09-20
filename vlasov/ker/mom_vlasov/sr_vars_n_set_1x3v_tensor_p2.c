#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_n_set_1x3v_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i) 
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
  gkyl_mat_set(&A0,0,0,0.7071067811865475*M0[0]); 
  gkyl_mat_set(&A0,0,1,0.7071067811865475*M0[1]); 
  gkyl_mat_set(&A0,0,2,0.7071067811865475*M0[2]); 
  gkyl_mat_set(&A0,1,0,0.7071067811865475*M0[1]); 
  gkyl_mat_set(&A0,1,1,0.6324555320336759*M0[2]+0.7071067811865475*M0[0]); 
  gkyl_mat_set(&A0,1,2,0.6324555320336759*M0[1]); 
  gkyl_mat_set(&A0,2,0,0.7071067811865475*M0[2]); 
  gkyl_mat_set(&A0,2,1,0.6324555320336759*M0[1]); 
  gkyl_mat_set(&A0,2,2,0.45175395145262565*M0[2]+0.7071067811865475*M0[0]); 
 
  struct gkyl_mat A1 = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat rhs1 = gkyl_nmat_get(rhs, count+1); 
  gkyl_mat_clear(&A1, 0.0); gkyl_mat_clear(&rhs1, 0.0); 
  const double *M11 = &M1i[3]; 
  gkyl_mat_set(&rhs1,0,0,M11[0]); 
  gkyl_mat_set(&rhs1,1,0,M11[1]); 
  gkyl_mat_set(&rhs1,2,0,M11[2]); 
  gkyl_mat_set(&A1,0,0,0.7071067811865475*M0[0]); 
  gkyl_mat_set(&A1,0,1,0.7071067811865475*M0[1]); 
  gkyl_mat_set(&A1,0,2,0.7071067811865475*M0[2]); 
  gkyl_mat_set(&A1,1,0,0.7071067811865475*M0[1]); 
  gkyl_mat_set(&A1,1,1,0.6324555320336759*M0[2]+0.7071067811865475*M0[0]); 
  gkyl_mat_set(&A1,1,2,0.6324555320336759*M0[1]); 
  gkyl_mat_set(&A1,2,0,0.7071067811865475*M0[2]); 
  gkyl_mat_set(&A1,2,1,0.6324555320336759*M0[1]); 
  gkyl_mat_set(&A1,2,2,0.45175395145262565*M0[2]+0.7071067811865475*M0[0]); 
 
  struct gkyl_mat A2 = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat rhs2 = gkyl_nmat_get(rhs, count+2); 
  gkyl_mat_clear(&A2, 0.0); gkyl_mat_clear(&rhs2, 0.0); 
  const double *M12 = &M1i[6]; 
  gkyl_mat_set(&rhs2,0,0,M12[0]); 
  gkyl_mat_set(&rhs2,1,0,M12[1]); 
  gkyl_mat_set(&rhs2,2,0,M12[2]); 
  gkyl_mat_set(&A2,0,0,0.7071067811865475*M0[0]); 
  gkyl_mat_set(&A2,0,1,0.7071067811865475*M0[1]); 
  gkyl_mat_set(&A2,0,2,0.7071067811865475*M0[2]); 
  gkyl_mat_set(&A2,1,0,0.7071067811865475*M0[1]); 
  gkyl_mat_set(&A2,1,1,0.6324555320336759*M0[2]+0.7071067811865475*M0[0]); 
  gkyl_mat_set(&A2,1,2,0.6324555320336759*M0[1]); 
  gkyl_mat_set(&A2,2,0,0.7071067811865475*M0[2]); 
  gkyl_mat_set(&A2,2,1,0.6324555320336759*M0[1]); 
  gkyl_mat_set(&A2,2,2,0.45175395145262565*M0[2]+0.7071067811865475*M0[0]); 
 
} 
