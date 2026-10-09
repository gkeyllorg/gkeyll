GKYL_CU_DH static inline void 
ser_2x_p0_inv(const double *A, double *A_inv) 
{ 
  // A:     Input DG field. 
  // A_inv: Output DG field (expansion of 1/A). 
 
  A_inv[0] = 1.0/A[0]; 

} 
 
