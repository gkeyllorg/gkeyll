GKYL_CU_DH static inline void 
ser_2x_p1_sqrt(const double *A, double *ASqrt) 
{ 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A)). 
 
  double AOrd[4] = {0.0}; 

  double temp = 0.0; 
  temp = 0.5*A[3]-0.5*(A[2]+A[1])+0.5*A[0]; 
  if (temp < 0.0) { 
  AOrd[0] = 0.0; 
  } else { 
  AOrd[0] = sqrt(temp); 
  } 
  temp = -(0.5*A[3])+0.5*A[2]-0.5*A[1]+0.5*A[0]; 
  if (temp < 0.0) { 
  AOrd[1] = 0.0; 
  } else { 
  AOrd[1] = sqrt(temp); 
  } 
  temp = 0.5*(A[1]+A[0])-0.5*(A[3]+A[2]); 
  if (temp < 0.0) { 
  AOrd[2] = 0.0; 
  } else { 
  AOrd[2] = sqrt(temp); 
  } 
  temp = 0.5*(A[3]+A[2]+A[1]+A[0]); 
  if (temp < 0.0) { 
  AOrd[3] = 0.0; 
  } else { 
  AOrd[3] = sqrt(temp); 
  } 
  ASqrt[0] = 0.5*AOrd[3]+0.5*AOrd[2]+0.5*AOrd[1]+0.5*AOrd[0]; 
  ASqrt[1] = 0.5*AOrd[3]+0.5*AOrd[2]-0.5*AOrd[1]-0.5*AOrd[0]; 
  ASqrt[2] = 0.5*AOrd[3]-0.5*AOrd[2]+0.5*AOrd[1]-0.5*AOrd[0]; 
  ASqrt[3] = 0.5*AOrd[3]-0.5*AOrd[2]-0.5*AOrd[1]+0.5*AOrd[0]; 

} 
 
