#include <gkyl_canonical_pb_vars_kernels.h>  
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void canonical_pb_vars_pressure_1x2v_ser_p2(const double *h_ij_inv, const double *MEnergy, const double *v_j, const double *nv_i, double* GKYL_RESTRICT d_Jv_P) 
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
  const double *NVy = &nv_i[3]; 

  const double *Vx = &v_j[0]; 
  const double *Vy = &v_j[3]; 

  const double *energy = &MEnergy[0]; 

  const double *Hxx = &h_ij_inv[0]; 
  const double *Hxy = &h_ij_inv[3]; 
  const double *Hyy = &h_ij_inv[6]; 

  // h^{ij}*nv_i*v_j 
  double Hxx_M1x[3] = {0.0}; 
  double Hxx_M1x_Vx[3] = {0.0}; 
  binop_mul_1d_ser_p2(Hxx, NVx, Hxx_M1x); 
  binop_mul_1d_ser_p2(Hxx_M1x, Vx, Hxx_M1x_Vx); 
 
  double Hxy_M1x[3] = {0.0}; 
  double Hxy_M1x_Vy[3] = {0.0}; 
  binop_mul_1d_ser_p2(Hxy, NVx, Hxy_M1x); 
  binop_mul_1d_ser_p2(Hxy_M1x, Vy, Hxy_M1x_Vy); 
 
  double Hyy_M1y[3] = {0.0}; 
  double Hyy_M1y_Vy[3] = {0.0}; 
  binop_mul_1d_ser_p2(Hyy, NVy, Hyy_M1y); 
  binop_mul_1d_ser_p2(Hyy_M1y, Vy, Hyy_M1y_Vy); 
 
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
 
} 
