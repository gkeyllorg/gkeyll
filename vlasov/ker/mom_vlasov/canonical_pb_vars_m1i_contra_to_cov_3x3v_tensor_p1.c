#include <gkyl_canonical_pb_vars_kernels.h>  
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_3x3v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov) 
{ 
  // h_ij:         Input volume expansion of the covariant components of the metric tensor.
  //                   [Hxx, Hxy, Hxz, 
  //                    Hxy, Hyy, Hyz, 
  //                    Hxz, Hyz, Hzz] 
  // v_j:              Input volume expansion of V_drift (contravariant components).
  //                   [vx, vy, vz] 
  // nv_i:              Input volume expansion of M1i = N*Vdrift (contravariant components).
  //                   [nvx, nvy, nvz] 
  // v_j_cov:           Input volume expansion of V_drift (covariant components).
  //                   [vx, vy, vz] 
  // nv_i_cov:          Input volume expansion of M1i = N*Vdrift (covariant components).
  //                   [nvx, nvy, nvz] 

  const double *NVx = &nv_i[0]; 
  const double *NVy = &nv_i[8]; 
  const double *NVz = &nv_i[16]; 

  const double *Vx = &v_j[0]; 
  const double *Vy = &v_j[8]; 
  const double *Vz = &v_j[16]; 

  double *NVx_cov = &nv_i_cov[0]; 
  double *NVy_cov = &nv_i_cov[8]; 
  double *NVz_cov = &nv_i_cov[16]; 

  double *Vx_cov = &v_j_cov[0]; 
  double *Vy_cov = &v_j_cov[8]; 
  double *Vz_cov = &v_j_cov[16]; 

  const double *Hxx = &h_ij[0]; 
  const double *Hxy = &h_ij[8]; 
  const double *Hxz = &h_ij[16]; 
  const double *Hyy = &h_ij[24]; 
  const double *Hyz = &h_ij[32]; 
  const double *Hzz = &h_ij[40]; 

  // h_{ij}v^j 
  double Hxx_Vx[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxx, Vx, Hxx_Vx); 
 
  double Hxy_Vx[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxy, Vx, Hxy_Vx); 
 
  double Hxy_Vy[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxy, Vy, Hxy_Vy); 
 
  double Hyy_Vy[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hyy, Vy, Hyy_Vy); 
 
  double Hxz_Vx[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxz, Vx, Hxz_Vx); 
 
  double Hyz_Vy[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hyz, Vy, Hyz_Vy); 
 
  double Hxz_Vz[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxz, Vz, Hxz_Vz); 
 
  double Hyz_Vz[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hyz, Vz, Hyz_Vz); 
 
  double Hzz_Vz[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hzz, Vz, Hzz_Vz); 
 
  // h_{ij}Jnv^j 
  double Hxx_NVx[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxx, NVx, Hxx_NVx); 
 
  double Hxy_NVx[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxy, NVx, Hxy_NVx); 
 
  double Hxy_NVy[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxy, NVy, Hxy_NVy); 
 
  double Hyy_NVy[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hyy, NVy, Hyy_NVy); 
 
  double Hxz_NVx[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxz, NVx, Hxz_NVx); 
 
  double Hyz_NVy[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hyz, NVy, Hyz_NVy); 
 
  double Hxz_NVz[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hxz, NVz, Hxz_NVz); 
 
  double Hyz_NVz[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hyz, NVz, Hyz_NVz); 
 
  double Hzz_NVz[8] = {0.0}; 
  binop_mul_3d_ser_p1(Hzz, NVz, Hzz_NVz); 
 
  // u_i_cov = h_{ij}v^j 
  Vx_cov[0] = Hxx_Vx[0]; 
  Vx_cov[0] += Hxy_Vy[0]; 
  Vy_cov[0] = Hxy_Vx[0]; 
  Vy_cov[0] += Hyy_Vy[0]; 
  Vx_cov[0] += Hxz_Vz[0]; 
  Vy_cov[0] += Hyz_Vz[0]; 
  Vz_cov[0] = Hxz_Vx[0]; 
  Vz_cov[0] += Hyz_Vy[0]; 
  Vz_cov[0] += Hzz_Vz[0]; 
  Vx_cov[1] = Hxx_Vx[1]; 
  Vx_cov[1] += Hxy_Vy[1]; 
  Vy_cov[1] = Hxy_Vx[1]; 
  Vy_cov[1] += Hyy_Vy[1]; 
  Vx_cov[1] += Hxz_Vz[1]; 
  Vy_cov[1] += Hyz_Vz[1]; 
  Vz_cov[1] = Hxz_Vx[1]; 
  Vz_cov[1] += Hyz_Vy[1]; 
  Vz_cov[1] += Hzz_Vz[1]; 
  Vx_cov[2] = Hxx_Vx[2]; 
  Vx_cov[2] += Hxy_Vy[2]; 
  Vy_cov[2] = Hxy_Vx[2]; 
  Vy_cov[2] += Hyy_Vy[2]; 
  Vx_cov[2] += Hxz_Vz[2]; 
  Vy_cov[2] += Hyz_Vz[2]; 
  Vz_cov[2] = Hxz_Vx[2]; 
  Vz_cov[2] += Hyz_Vy[2]; 
  Vz_cov[2] += Hzz_Vz[2]; 
  Vx_cov[3] = Hxx_Vx[3]; 
  Vx_cov[3] += Hxy_Vy[3]; 
  Vy_cov[3] = Hxy_Vx[3]; 
  Vy_cov[3] += Hyy_Vy[3]; 
  Vx_cov[3] += Hxz_Vz[3]; 
  Vy_cov[3] += Hyz_Vz[3]; 
  Vz_cov[3] = Hxz_Vx[3]; 
  Vz_cov[3] += Hyz_Vy[3]; 
  Vz_cov[3] += Hzz_Vz[3]; 
  Vx_cov[4] = Hxx_Vx[4]; 
  Vx_cov[4] += Hxy_Vy[4]; 
  Vy_cov[4] = Hxy_Vx[4]; 
  Vy_cov[4] += Hyy_Vy[4]; 
  Vx_cov[4] += Hxz_Vz[4]; 
  Vy_cov[4] += Hyz_Vz[4]; 
  Vz_cov[4] = Hxz_Vx[4]; 
  Vz_cov[4] += Hyz_Vy[4]; 
  Vz_cov[4] += Hzz_Vz[4]; 
  Vx_cov[5] = Hxx_Vx[5]; 
  Vx_cov[5] += Hxy_Vy[5]; 
  Vy_cov[5] = Hxy_Vx[5]; 
  Vy_cov[5] += Hyy_Vy[5]; 
  Vx_cov[5] += Hxz_Vz[5]; 
  Vy_cov[5] += Hyz_Vz[5]; 
  Vz_cov[5] = Hxz_Vx[5]; 
  Vz_cov[5] += Hyz_Vy[5]; 
  Vz_cov[5] += Hzz_Vz[5]; 
  Vx_cov[6] = Hxx_Vx[6]; 
  Vx_cov[6] += Hxy_Vy[6]; 
  Vy_cov[6] = Hxy_Vx[6]; 
  Vy_cov[6] += Hyy_Vy[6]; 
  Vx_cov[6] += Hxz_Vz[6]; 
  Vy_cov[6] += Hyz_Vz[6]; 
  Vz_cov[6] = Hxz_Vx[6]; 
  Vz_cov[6] += Hyz_Vy[6]; 
  Vz_cov[6] += Hzz_Vz[6]; 
  Vx_cov[7] = Hxx_Vx[7]; 
  Vx_cov[7] += Hxy_Vy[7]; 
  Vy_cov[7] = Hxy_Vx[7]; 
  Vy_cov[7] += Hyy_Vy[7]; 
  Vx_cov[7] += Hxz_Vz[7]; 
  Vy_cov[7] += Hyz_Vz[7]; 
  Vz_cov[7] = Hxz_Vx[7]; 
  Vz_cov[7] += Hyz_Vy[7]; 
  Vz_cov[7] += Hzz_Vz[7]; 
 
  // Jnu_i_cov = h_{ij}Jnv^j 
  NVx_cov[0] = Hxx_NVx[0]; 
  NVx_cov[0] += Hxy_NVy[0]; 
  NVy_cov[0] = Hxy_NVx[0]; 
  NVy_cov[0] += Hyy_NVy[0]; 
  NVx_cov[0] += Hxz_NVz[0]; 
  NVy_cov[0] += Hyz_NVz[0]; 
  NVz_cov[0] = Hxz_NVx[0]; 
  NVz_cov[0] += Hyz_NVy[0]; 
  NVz_cov[0] += Hzz_NVz[0]; 
  NVx_cov[1] = Hxx_NVx[1]; 
  NVx_cov[1] += Hxy_NVy[1]; 
  NVy_cov[1] = Hxy_NVx[1]; 
  NVy_cov[1] += Hyy_NVy[1]; 
  NVx_cov[1] += Hxz_NVz[1]; 
  NVy_cov[1] += Hyz_NVz[1]; 
  NVz_cov[1] = Hxz_NVx[1]; 
  NVz_cov[1] += Hyz_NVy[1]; 
  NVz_cov[1] += Hzz_NVz[1]; 
  NVx_cov[2] = Hxx_NVx[2]; 
  NVx_cov[2] += Hxy_NVy[2]; 
  NVy_cov[2] = Hxy_NVx[2]; 
  NVy_cov[2] += Hyy_NVy[2]; 
  NVx_cov[2] += Hxz_NVz[2]; 
  NVy_cov[2] += Hyz_NVz[2]; 
  NVz_cov[2] = Hxz_NVx[2]; 
  NVz_cov[2] += Hyz_NVy[2]; 
  NVz_cov[2] += Hzz_NVz[2]; 
  NVx_cov[3] = Hxx_NVx[3]; 
  NVx_cov[3] += Hxy_NVy[3]; 
  NVy_cov[3] = Hxy_NVx[3]; 
  NVy_cov[3] += Hyy_NVy[3]; 
  NVx_cov[3] += Hxz_NVz[3]; 
  NVy_cov[3] += Hyz_NVz[3]; 
  NVz_cov[3] = Hxz_NVx[3]; 
  NVz_cov[3] += Hyz_NVy[3]; 
  NVz_cov[3] += Hzz_NVz[3]; 
  NVx_cov[4] = Hxx_NVx[4]; 
  NVx_cov[4] += Hxy_NVy[4]; 
  NVy_cov[4] = Hxy_NVx[4]; 
  NVy_cov[4] += Hyy_NVy[4]; 
  NVx_cov[4] += Hxz_NVz[4]; 
  NVy_cov[4] += Hyz_NVz[4]; 
  NVz_cov[4] = Hxz_NVx[4]; 
  NVz_cov[4] += Hyz_NVy[4]; 
  NVz_cov[4] += Hzz_NVz[4]; 
  NVx_cov[5] = Hxx_NVx[5]; 
  NVx_cov[5] += Hxy_NVy[5]; 
  NVy_cov[5] = Hxy_NVx[5]; 
  NVy_cov[5] += Hyy_NVy[5]; 
  NVx_cov[5] += Hxz_NVz[5]; 
  NVy_cov[5] += Hyz_NVz[5]; 
  NVz_cov[5] = Hxz_NVx[5]; 
  NVz_cov[5] += Hyz_NVy[5]; 
  NVz_cov[5] += Hzz_NVz[5]; 
  NVx_cov[6] = Hxx_NVx[6]; 
  NVx_cov[6] += Hxy_NVy[6]; 
  NVy_cov[6] = Hxy_NVx[6]; 
  NVy_cov[6] += Hyy_NVy[6]; 
  NVx_cov[6] += Hxz_NVz[6]; 
  NVy_cov[6] += Hyz_NVz[6]; 
  NVz_cov[6] = Hxz_NVx[6]; 
  NVz_cov[6] += Hyz_NVy[6]; 
  NVz_cov[6] += Hzz_NVz[6]; 
  NVx_cov[7] = Hxx_NVx[7]; 
  NVx_cov[7] += Hxy_NVy[7]; 
  NVy_cov[7] = Hxy_NVx[7]; 
  NVy_cov[7] += Hyy_NVy[7]; 
  NVx_cov[7] += Hxz_NVz[7]; 
  NVy_cov[7] += Hyz_NVz[7]; 
  NVz_cov[7] = Hxz_NVx[7]; 
  NVz_cov[7] += Hyz_NVy[7]; 
  NVz_cov[7] += Hzz_NVz[7]; 
 
} 
