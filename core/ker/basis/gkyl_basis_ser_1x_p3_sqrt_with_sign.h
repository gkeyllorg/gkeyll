GKYL_CU_DH static inline void 
ser_1x_p3_sqrt_with_sign(const double *ASign, const double *A, double *ASqrt) 
{ 
  // ASign: Input DG field, used to get correct sign of Asqrt. 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A), with sign determined by Asign). 
 
  double AOrd[4] = {0.0}; 

  double temp = 0.0; 
  double temp_sign = 0.0; 
  temp = -(0.5701294036773671*A[3])+0.9681844646844028*A[2]-1.0546722811938851*A[1]+0.7071067811865475*A[0]; 
  temp_sign = -(0.5701294036773671*ASign[3])+0.9681844646844028*ASign[2]-1.0546722811938851*ASign[1]+0.7071067811865475*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[0] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[0] = -sqrt(temp); 
  } else { 
  AOrd[0] = sqrt(temp); 
  } 
  temp = 0.7702725556588816*A[3]-0.5164305132317774*A[2]-0.41639003950091297*A[1]+0.7071067811865475*A[0]; 
  temp_sign = 0.7702725556588816*ASign[3]-0.5164305132317774*ASign[2]-0.41639003950091297*ASign[1]+0.7071067811865475*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[1] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[1] = -sqrt(temp); 
  } else { 
  AOrd[1] = sqrt(temp); 
  } 
  temp = -(0.7702725556588816*A[3])-0.5164305132317774*A[2]+0.41639003950091297*A[1]+0.7071067811865475*A[0]; 
  temp_sign = -(0.7702725556588816*ASign[3])-0.5164305132317774*ASign[2]+0.41639003950091297*ASign[1]+0.7071067811865475*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[2] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[2] = -sqrt(temp); 
  } else { 
  AOrd[2] = sqrt(temp); 
  } 
  temp = 0.5701294036773671*A[3]+0.9681844646844028*A[2]+1.0546722811938851*A[1]+0.7071067811865475*A[0]; 
  temp_sign = 0.5701294036773671*ASign[3]+0.9681844646844028*ASign[2]+1.0546722811938851*ASign[1]+0.7071067811865475*ASign[0]; 
  if (temp < 0.0) { 
  AOrd[3] = 0.0; 
  } else if (temp > 0.0 && temp_sign < 0.0) { 
  AOrd[3] = -sqrt(temp); 
  } else { 
  AOrd[3] = sqrt(temp); 
  } 
  ASqrt[0] = 0.24597051986529*AOrd[3]+0.46113626132125757*AOrd[2]+0.46113626132125757*AOrd[1]+0.24597051986529*AOrd[0]; 
  ASqrt[1] = 0.3668728630454641*AOrd[3]+0.271546746793545*AOrd[2]-0.271546746793545*AOrd[1]-0.3668728630454641*AOrd[0]; 
  ASqrt[2] = 0.3367876570272824*AOrd[3]-0.3367876570272824*AOrd[2]-0.3367876570272824*AOrd[1]+0.3367876570272824*AOrd[0]; 
  ASqrt[3] = 0.19832227542449957*AOrd[3]-0.5023295150965308*AOrd[2]+0.5023295150965308*AOrd[1]-0.19832227542449957*AOrd[0]; 

} 
 
