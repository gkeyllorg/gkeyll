#include <gkyl_mat.h> 
#include <gkyl_binop_div_ser.h> 
 
GKYL_CU_DH void binop_div_set_2d_ser_p0(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *f, const double *g) 
{ 
  // A:       preallocated LHS matrix. 
  // rhs:     preallocated RHS vector. 
  // f:       numerator field (must be a scalar). 
  // g:       denominator field (must be a scalar). 
 
  double lhs[1]; 
  lhs[0] = g[0]; 
  gkyl_mat_set(rhs,0,0,f[0]); 
 
  // Fill LHS matrix. 
  gkyl_mat_set(A,0,0,0.5*lhs[0]); 

} 
