#include <gkyl_canonical_pb_vars_kernels.h>  
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void canonical_pb_vars_pressure_2x2v_tensor_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, const double *nv_i, double* GKYL_RESTRICT d_Jv_P) 
{ 
  // h_ij_inv:         Input volume expansion of the inverse metric tensor.
  //                   [Hxx, Hxy, Hxz, 
  //                    Hxy, Hyy, Hyz, 
  //                    Hxz, Hyz, Hzz] 
  // MEnergy:          Input volume expansion of the MEnergy moment.
  // v_j:              Input volume expansion of V_drift.
  //                   [vx, vy, vz] 
  // nv_i:              Input volume expansion of M1i = N*Vdrift.
  //                   [nvx, nvy, nvz] 
  // d_Jv_P:            Output volume expansion of d*P*Jv = h^{ij}*M2_{ij} - n*h^{ij}*u_i*u_j .

  const double *NVx = &nv_i[0]; 
  const double *NVy = &nv_i[9]; 

  const double *Vx = &v_j[0]; 
  const double *Vy = &v_j[9]; 

  const double *energy = &MEnergy[0]; 

  const double *Hxx = &h_ij_inv[0]; 
  const double *Hxy = &h_ij_inv[9]; 
  const double *Hyy = &h_ij_inv[18]; 

  // h^{ij}*nv_i*v_j 
  double Hxx_M1x[9] = {0.0}; 
  double Hxx_M1x_Vx[9] = {0.0}; 
  binop_mul_2d_tensor_p2(Hxx, NVx, Hxx_M1x); 
  binop_mul_2d_tensor_p2(Hxx_M1x, Vx, Hxx_M1x_Vx); 
 
  double Hxy_M1x[9] = {0.0}; 
  double Hxy_M1x_Vy[9] = {0.0}; 
  binop_mul_2d_tensor_p2(Hxy, NVx, Hxy_M1x); 
  binop_mul_2d_tensor_p2(Hxy_M1x, Vy, Hxy_M1x_Vy); 
 
  double Hyy_M1y[9] = {0.0}; 
  double Hyy_M1y_Vy[9] = {0.0}; 
  binop_mul_2d_tensor_p2(Hyy, NVy, Hyy_M1y); 
  binop_mul_2d_tensor_p2(Hyy_M1y, Vy, Hyy_M1y_Vy); 
 
  d_Jv_P[0] = 2.0*energy[0]; 
  d_Jv_P[0] += - Hxx_M1x_Vx[0]; 
  d_Jv_P[0] += (- Hxy_M1x_Vy[0])*2.0; 
  d_Jv_P[0] +=  - Hyy_M1y_Vy[0]; 
  d_Jv_P[1] = 2.0*energy[1]; 
  d_Jv_P[1] += - Hxx_M1x_Vx[1]; 
  d_Jv_P[1] += (- Hxy_M1x_Vy[1])*2.0; 
  d_Jv_P[1] +=  - Hyy_M1y_Vy[1]; 
  d_Jv_P[2] = 2.0*energy[2]; 
  d_Jv_P[2] += - Hxx_M1x_Vx[2]; 
  d_Jv_P[2] += (- Hxy_M1x_Vy[2])*2.0; 
  d_Jv_P[2] +=  - Hyy_M1y_Vy[2]; 
  d_Jv_P[3] = 2.0*energy[3]; 
  d_Jv_P[3] += - Hxx_M1x_Vx[3]; 
  d_Jv_P[3] += (- Hxy_M1x_Vy[3])*2.0; 
  d_Jv_P[3] +=  - Hyy_M1y_Vy[3]; 
  d_Jv_P[4] = 2.0*energy[4]; 
  d_Jv_P[4] += - Hxx_M1x_Vx[4]; 
  d_Jv_P[4] += (- Hxy_M1x_Vy[4])*2.0; 
  d_Jv_P[4] +=  - Hyy_M1y_Vy[4]; 
  d_Jv_P[5] = 2.0*energy[5]; 
  d_Jv_P[5] += - Hxx_M1x_Vx[5]; 
  d_Jv_P[5] += (- Hxy_M1x_Vy[5])*2.0; 
  d_Jv_P[5] +=  - Hyy_M1y_Vy[5]; 
  d_Jv_P[6] = 2.0*energy[6]; 
  d_Jv_P[6] += - Hxx_M1x_Vx[6]; 
  d_Jv_P[6] += (- Hxy_M1x_Vy[6])*2.0; 
  d_Jv_P[6] +=  - Hyy_M1y_Vy[6]; 
  d_Jv_P[7] = 2.0*energy[7]; 
  d_Jv_P[7] += - Hxx_M1x_Vx[7]; 
  d_Jv_P[7] += (- Hxy_M1x_Vy[7])*2.0; 
  d_Jv_P[7] +=  - Hyy_M1y_Vy[7]; 
  d_Jv_P[8] = 2.0*energy[8]; 
  d_Jv_P[8] += - Hxx_M1x_Vx[8]; 
  d_Jv_P[8] += (- Hxy_M1x_Vy[8])*2.0; 
  d_Jv_P[8] +=  - Hyy_M1y_Vy[8]; 
 
} 
