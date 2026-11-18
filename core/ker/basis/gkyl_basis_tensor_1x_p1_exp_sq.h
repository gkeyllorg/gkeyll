GKYL_CU_DH static inline void 
tensor_1x_p1_exp_sq(const double *A, double *ASq) 
{ 
  // A:   Input DG field. 
  // ASq: Output DG field (expansion of A^2). 
 
  const double A0R2 = pow(A[0],2);
  const double A1R2 = pow(A[1],2);

  ASq[0] = 0.7071067811865475*A1R2+0.7071067811865475*A0R2; 
  ASq[1] = 1.4142135623730951*A[0]*A[1]; 
} 
 
