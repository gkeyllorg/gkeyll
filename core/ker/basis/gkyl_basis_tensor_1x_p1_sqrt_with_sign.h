GKYL_CU_DH static inline void 
tensor_1x_p1_sqrt_with_sign(const double *ASign, const double *A, double *ASqrt) 
{ 
  // ASign: Input DG field, used to get correct sign of Asqrt. 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A), with sign determined by Asign). 
 
  double AOrd[2] = {0.0}; 

  double temp = 0.0; 
  double temp_sign = 0.0; 
  temp = 0.7071067811865475*A[0]-0.7071067811865475*A[1]; 
  temp_sign = 0.7071067811865475*ASign[0]-0.7071067811865475*ASign[1]; 
  if (temp < 0.0) { 
  AOrd[0] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[0] = -sqrt(temp); 
  } else { 
  AOrd[0] = sqrt(temp); 
  } 
  temp = 0.7071067811865475*(A[1]+A[0]); 
  temp_sign = 0.7071067811865475*(ASign[1]+ASign[0]); 
  if (temp < 0.0) { 
  AOrd[1] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[1] = -sqrt(temp); 
  } else { 
  AOrd[1] = sqrt(temp); 
  } 
  ASqrt[0] = 0.7071067811865475*AOrd[1]+0.7071067811865475*AOrd[0]; 
  ASqrt[1] = 0.7071067811865475*AOrd[1]-0.7071067811865475*AOrd[0]; 

} 
 
