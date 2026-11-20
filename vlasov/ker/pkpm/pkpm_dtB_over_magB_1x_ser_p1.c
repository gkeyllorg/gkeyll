#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH int pkpm_dtB_over_magB_1x_ser_p1(const double *em, const double *dt_em, double* GKYL_RESTRICT out) 
{ 
  // em:    Input physical electromagnetic fields [Ex, Ey, Ez, Bx, By, Bz, phi, psi]. 
  // dt_em: Input time derivative of the physical fields (the Maxwell right-hand side), same layout. 
  // out:   Output d_t B/|B| [x, y, z] at the 2 Gauss-Legendre nodes (where b = B/|B| is defined), 
  //        interpolated to the p=1 basis. 
  // Returns the number of nodes where |B| = 0 (the output is zero there). 
  const double *B0 = &em[6]; 
  const double *dB0 = &dt_em[6]; 
  double *o0 = &out[0]; 
  const double *B1 = &em[8]; 
  const double *dB1 = &dt_em[8]; 
  double *o1 = &out[2]; 
  const double *B2 = &em[10]; 
  const double *dB2 = &dt_em[10]; 
  double *o2 = &out[4]; 
  double q[3][2] = {{0.0}}; 
  int num_null = 0; 
  { 
    const double B[3] = {0.7071067811865475*B0[0]-0.7071067811865475*B0[1], 0.7071067811865475*B1[0]-0.7071067811865475*B1[1], 0.7071067811865475*B2[0]-0.7071067811865475*B2[1]}; 
    const double dB[3] = {0.7071067811865475*dB0[0]-0.7071067811865475*dB0[1], 0.7071067811865475*dB1[0]-0.7071067811865475*dB1[1], 0.7071067811865475*dB2[0]-0.7071067811865475*dB2[1]}; 
    const double magB = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]); 
    if (magB > 0.0) { 
      for (int i = 0; i < 3; ++i) q[i][0] = dB[i]/magB; 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double B[3] = {0.7071067811865475*B0[1]+0.7071067811865475*B0[0], 0.7071067811865475*B1[1]+0.7071067811865475*B1[0], 0.7071067811865475*B2[1]+0.7071067811865475*B2[0]}; 
    const double dB[3] = {0.7071067811865475*dB0[1]+0.7071067811865475*dB0[0], 0.7071067811865475*dB1[1]+0.7071067811865475*dB1[0], 0.7071067811865475*dB2[1]+0.7071067811865475*dB2[0]}; 
    const double magB = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]); 
    if (magB > 0.0) { 
      for (int i = 0; i < 3; ++i) q[i][1] = dB[i]/magB; 
    } else { 
      num_null += 1; 
    } 
  } 
  o0[0] = 0.7071067811865475*q[0][1]+0.7071067811865475*q[0][0]; 
  o0[1] = 0.7071067811865475*q[0][1]-0.7071067811865475*q[0][0]; 
  o1[0] = 0.7071067811865475*q[1][1]+0.7071067811865475*q[1][0]; 
  o1[1] = 0.7071067811865475*q[1][1]-0.7071067811865475*q[1][0]; 
  o2[0] = 0.7071067811865475*q[2][1]+0.7071067811865475*q[2][0]; 
  o2[1] = 0.7071067811865475*q[2][1]-0.7071067811865475*q[2][0]; 
  return num_null; 
} 
