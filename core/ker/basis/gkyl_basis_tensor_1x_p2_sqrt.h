GKYL_CU_DH static inline void 
tensor_1x_p2_sqrt(const double *A, double *ASqrt) 
{ 
  // A:     Input DG field. 
  // ASqrt: Output DG field (expansion of sqrt(A)). 
 
  double AOrd[3] = {0.0}; 

  double temp = 0.0; 
  temp = 0.6324555320336759*A[2]-0.9486832980505137*A[1]+0.7071067811865475*A[0]; 
  if (temp < 0.0) { 
  AOrd[0] = 0.0; 
  } else { 
  AOrd[0] = sqrt(temp); 
  } 
  temp = 0.7071067811865475*A[0]-0.7905694150420947*A[2]; 
  if (temp < 0.0) { 
  AOrd[1] = 0.0; 
  } else { 
  AOrd[1] = sqrt(temp); 
  } 
  temp = 0.6324555320336759*A[2]+0.9486832980505137*A[1]+0.7071067811865475*A[0]; 
  if (temp < 0.0) { 
  AOrd[2] = 0.0; 
  } else { 
  AOrd[2] = sqrt(temp); 
  } 
  ASqrt[0] = 0.39283710065919303*AOrd[2]+0.6285393610547091*AOrd[1]+0.39283710065919303*AOrd[0]; 
  ASqrt[1] = 0.5270462766947298*AOrd[2]-0.5270462766947298*AOrd[0]; 
  ASqrt[2] = 0.35136418446315326*AOrd[2]-0.7027283689263066*AOrd[1]+0.35136418446315326*AOrd[0]; 

} 
 
