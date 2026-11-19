GKYL_CU_DH static inline void 
tensor_1x_p2_exp_sq(const double *A, double *ASq) 
{ 
  // A:   Input DG field. 
  // ASq: Output DG field (expansion of A^2). 
 
  const double A0R2 = pow(A[0],2);
  const double A1R2 = pow(A[1],2);
  const double A2R2 = pow(A[2],2);

  ASq[0] = 0.7071067811865475*A2R2+0.7071067811865475*A1R2+0.7071067811865475*A0R2; 
  ASq[1] = 1.264911064067352*A[1]*A[2]+1.4142135623730951*A[0]*A[1]; 
  ASq[2] = 0.45175395145262565*A2R2+1.4142135623730951*A[0]*A[2]+0.6324555320336759*A1R2; 
} 
 
