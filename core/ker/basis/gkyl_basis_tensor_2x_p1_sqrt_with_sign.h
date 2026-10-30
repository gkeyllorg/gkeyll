GKYL_CU_DH static inline void 
tensor_2x_p1_sqrt_with_sign(const double *ASign, const double *A, double *ASqrt) 
{ 
  // ASign: Input DG field, used to get correct sign of Asqrt. 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A), with sign determined by Asign). 
 
  double AOrd[4] = {0.0}; 

  double temp = 0.0; 
  double temp_sign = 0.0; 
  temp = 0.5*A[3]-0.5*(A[2]+A[1])+0.5*A[0]; 
  temp_sign = 0.5*ASign[3]-0.5*(ASign[2]+ASign[1])+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[0] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[0] = -sqrt(temp); 
  } else { 
  AOrd[0] = sqrt(temp); 
  } 
  temp = -(0.5*A[3])+0.5*A[2]-0.5*A[1]+0.5*A[0]; 
  temp_sign = -(0.5*ASign[3])+0.5*ASign[2]-0.5*ASign[1]+0.5*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[1] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[1] = -sqrt(temp); 
  } else { 
  AOrd[1] = sqrt(temp); 
  } 
  temp = 0.5*(A[1]+A[0])-0.5*(A[3]+A[2]); 
  temp_sign = 0.5*(ASign[1]+ASign[0])-0.5*(ASign[3]+ASign[2]); 
  if (temp < 0.0) { 
  AOrd[2] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[2] = -sqrt(temp); 
  } else { 
  AOrd[2] = sqrt(temp); 
  } 
  temp = 0.5*(A[3]+A[2]+A[1]+A[0]); 
  temp_sign = 0.5*(ASign[3]+ASign[2]+ASign[1]+ASign[0]); 
  if (temp < 0.0) { 
  AOrd[3] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[3] = -sqrt(temp); 
  } else { 
  AOrd[3] = sqrt(temp); 
  } 
  ASqrt[0] = 0.5*AOrd[3]+0.5*AOrd[2]+0.5*AOrd[1]+0.5*AOrd[0]; 
  ASqrt[1] = 0.5*AOrd[3]+0.5*AOrd[2]-0.5*AOrd[1]-0.5*AOrd[0]; 
  ASqrt[2] = 0.5*AOrd[3]-0.5*AOrd[2]+0.5*AOrd[1]-0.5*AOrd[0]; 
  ASqrt[3] = 0.5*AOrd[3]-0.5*AOrd[2]-0.5*AOrd[1]+0.5*AOrd[0]; 

} 
 
