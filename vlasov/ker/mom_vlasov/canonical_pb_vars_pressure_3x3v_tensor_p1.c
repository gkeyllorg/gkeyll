#include <gkyl_canonical_pb_vars_kernels.h>  
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void canonical_pb_vars_pressure_3x3v_tensor_p1(const double *h_ij_inv, const double *MEnergy, const double *v_j, const double *nv_i, double* GKYL_RESTRICT d_Jv_P) 
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
  const double *NVy = &nv_i[8]; 
  const double *NVz = &nv_i[16]; 

  const double *Vx = &v_j[0]; 
  const double *Vy = &v_j[8]; 
  const double *Vz = &v_j[16]; 

  const double *energy = &MEnergy[0]; 

  const double *Hxx = &h_ij_inv[0]; 
  const double *Hxy = &h_ij_inv[8]; 
  const double *Hxz = &h_ij_inv[16]; 
  const double *Hyy = &h_ij_inv[24]; 
  const double *Hyz = &h_ij_inv[32]; 
  const double *Hzz = &h_ij_inv[40]; 

  // h^{ij}*nv_i*v_j 
  double Hxx_M1x[8] = {0.0}; 
  double Hxx_M1x_Vx[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxx, NVx, Hxx_M1x); 
  binop_mul_3d_ser_p1(Hxx_M1x, Vx, Hxx_M1x_Vx); 
 
  double Hxy_M1x[8] = {0.0}; 
  double Hxy_M1x_Vy[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxy, NVx, Hxy_M1x); 
  binop_mul_3d_ser_p1(Hxy_M1x, Vy, Hxy_M1x_Vy); 
 
  double Hyy_M1y[8] = {0.0}; 
  double Hyy_M1y_Vy[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hyy, NVy, Hyy_M1y); 
  binop_mul_3d_ser_p1(Hyy_M1y, Vy, Hyy_M1y_Vy); 
 
  double Hxz_M1x[8] = {0.0}; 
  double Hxz_M1x_Vz[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxz, NVx, Hxz_M1x); 
  binop_mul_3d_ser_p1(Hxz_M1x, Vz, Hxz_M1x_Vz); 
 
  double Hyz_M1y[8] = {0.0}; 
  double Hyz_M1y_Vz[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hyz, NVy, Hyz_M1y); 
  binop_mul_3d_ser_p1(Hyz_M1y, Vz, Hyz_M1y_Vz); 
 
  double Hzz_M1z[8] = {0.0}; 
  double Hzz_M1z_Vz[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hzz, NVz, Hzz_M1z); 
  binop_mul_3d_ser_p1(Hzz_M1z, Vz, Hzz_M1z_Vz); 
 
  d_Jv_P[0] = 2.0*energy[0]; 
  d_Jv_P[0] += - Hxx_M1x_Vx[0]; 
  d_Jv_P[0] += (- Hxy_M1x_Vy[0])*2.0; 
  d_Jv_P[0] +=  - Hyy_M1y_Vy[0]; 
  d_Jv_P[0] += (- Hxz_M1x_Vz[0])*2.0; 
  d_Jv_P[0] += (- Hyz_M1y_Vz[0])*2.0; 
  d_Jv_P[0] +=  - Hzz_M1z_Vz[0]; 
  d_Jv_P[1] = 2.0*energy[1]; 
  d_Jv_P[1] += - Hxx_M1x_Vx[1]; 
  d_Jv_P[1] += (- Hxy_M1x_Vy[1])*2.0; 
  d_Jv_P[1] +=  - Hyy_M1y_Vy[1]; 
  d_Jv_P[1] += (- Hxz_M1x_Vz[1])*2.0; 
  d_Jv_P[1] += (- Hyz_M1y_Vz[1])*2.0; 
  d_Jv_P[1] +=  - Hzz_M1z_Vz[1]; 
  d_Jv_P[2] = 2.0*energy[2]; 
  d_Jv_P[2] += - Hxx_M1x_Vx[2]; 
  d_Jv_P[2] += (- Hxy_M1x_Vy[2])*2.0; 
  d_Jv_P[2] +=  - Hyy_M1y_Vy[2]; 
  d_Jv_P[2] += (- Hxz_M1x_Vz[2])*2.0; 
  d_Jv_P[2] += (- Hyz_M1y_Vz[2])*2.0; 
  d_Jv_P[2] +=  - Hzz_M1z_Vz[2]; 
  d_Jv_P[3] = 2.0*energy[3]; 
  d_Jv_P[3] += - Hxx_M1x_Vx[3]; 
  d_Jv_P[3] += (- Hxy_M1x_Vy[3])*2.0; 
  d_Jv_P[3] +=  - Hyy_M1y_Vy[3]; 
  d_Jv_P[3] += (- Hxz_M1x_Vz[3])*2.0; 
  d_Jv_P[3] += (- Hyz_M1y_Vz[3])*2.0; 
  d_Jv_P[3] +=  - Hzz_M1z_Vz[3]; 
  d_Jv_P[4] = 2.0*energy[4]; 
  d_Jv_P[4] += - Hxx_M1x_Vx[4]; 
  d_Jv_P[4] += (- Hxy_M1x_Vy[4])*2.0; 
  d_Jv_P[4] +=  - Hyy_M1y_Vy[4]; 
  d_Jv_P[4] += (- Hxz_M1x_Vz[4])*2.0; 
  d_Jv_P[4] += (- Hyz_M1y_Vz[4])*2.0; 
  d_Jv_P[4] +=  - Hzz_M1z_Vz[4]; 
  d_Jv_P[5] = 2.0*energy[5]; 
  d_Jv_P[5] += - Hxx_M1x_Vx[5]; 
  d_Jv_P[5] += (- Hxy_M1x_Vy[5])*2.0; 
  d_Jv_P[5] +=  - Hyy_M1y_Vy[5]; 
  d_Jv_P[5] += (- Hxz_M1x_Vz[5])*2.0; 
  d_Jv_P[5] += (- Hyz_M1y_Vz[5])*2.0; 
  d_Jv_P[5] +=  - Hzz_M1z_Vz[5]; 
  d_Jv_P[6] = 2.0*energy[6]; 
  d_Jv_P[6] += - Hxx_M1x_Vx[6]; 
  d_Jv_P[6] += (- Hxy_M1x_Vy[6])*2.0; 
  d_Jv_P[6] +=  - Hyy_M1y_Vy[6]; 
  d_Jv_P[6] += (- Hxz_M1x_Vz[6])*2.0; 
  d_Jv_P[6] += (- Hyz_M1y_Vz[6])*2.0; 
  d_Jv_P[6] +=  - Hzz_M1z_Vz[6]; 
  d_Jv_P[7] = 2.0*energy[7]; 
  d_Jv_P[7] += - Hxx_M1x_Vx[7]; 
  d_Jv_P[7] += (- Hxy_M1x_Vy[7])*2.0; 
  d_Jv_P[7] +=  - Hyy_M1y_Vy[7]; 
  d_Jv_P[7] += (- Hxz_M1x_Vz[7])*2.0; 
  d_Jv_P[7] += (- Hyz_M1y_Vz[7])*2.0; 
  d_Jv_P[7] +=  - Hzz_M1z_Vz[7]; 
 
} 
