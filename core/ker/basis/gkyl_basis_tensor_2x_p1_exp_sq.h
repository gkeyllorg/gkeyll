GKYL_CU_DH static inline void 
tensor_2x_p1_exp_sq(const double *A, double *ASq) 
{ 
  // A:   Input DG field. 
  // ASq: Output DG field (expansion of A^2). 
 
  const double A0R2 = pow(A[0],2);
  const double A1R2 = pow(A[1],2);
  const double A2R2 = pow(A[2],2);
  const double A3R2 = pow(A[3],2);

  ASq[0] = 0.5*A3R2+0.5*A2R2+0.5*A1R2+0.5*A0R2; 
  ASq[1] = A[2]*A[3]+A[0]*A[1]; 
  ASq[2] = A[1]*A[3]+A[0]*A[2]; 
  ASq[3] = A[0]*A[3]+A[1]*A[2]; 
} 
 
