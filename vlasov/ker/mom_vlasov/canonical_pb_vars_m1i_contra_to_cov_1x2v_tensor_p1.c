#include <gkyl_canonical_pb_vars_kernels.h>  
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void canonical_pb_vars_m1i_contra_to_cov_1x2v_tensor_p1(const double *h_ij, const double *v_j, const double *nv_i, double* GKYL_RESTRICT v_j_cov, double* GKYL_RESTRICT nv_i_cov) 
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
  const double *NVy = &nv_i[2]; 

  const double *Vx = &v_j[0]; 
  const double *Vy = &v_j[2]; 

  double *NVx_cov = &nv_i_cov[0]; 
  double *NVy_cov = &nv_i_cov[2]; 

  double *Vx_cov = &v_j_cov[0]; 
  double *Vy_cov = &v_j_cov[2]; 

  const double *Hxx = &h_ij[0]; 
  const double *Hxy = &h_ij[2]; 
  const double *Hyy = &h_ij[4]; 

  // h_{ij}v^j 
  double Hxx_Vx[2] = {0.0}; 
  binop_mul_1d_ser_p1(Hxx, Vx, Hxx_Vx); 
 
  double Hxy_Vx[2] = {0.0}; 
  binop_mul_1d_ser_p1(Hxy, Vx, Hxy_Vx); 
 
  double Hxy_Vy[2] = {0.0}; 
  binop_mul_1d_ser_p1(Hxy, Vy, Hxy_Vy); 
 
  double Hyy_Vy[2] = {0.0}; 
  binop_mul_1d_ser_p1(Hyy, Vy, Hyy_Vy); 
 
  // h_{ij}Jnv^j 
  double Hxx_NVx[2] = {0.0}; 
  binop_mul_1d_ser_p1(Hxx, NVx, Hxx_NVx); 
 
  double Hxy_NVx[2] = {0.0}; 
  binop_mul_1d_ser_p1(Hxy, NVx, Hxy_NVx); 
 
  double Hxy_NVy[2] = {0.0}; 
  binop_mul_1d_ser_p1(Hxy, NVy, Hxy_NVy); 
 
  double Hyy_NVy[2] = {0.0}; 
  binop_mul_1d_ser_p1(Hyy, NVy, Hyy_NVy); 
 
  // u_i_cov = h_{ij}v^j 
  Vx_cov[0] = Hxx_Vx[0]; 
  Vx_cov[0] += Hxy_Vy[0]; 
  Vy_cov[0] = Hxy_Vx[0]; 
  Vy_cov[0] += Hyy_Vy[0]; 
  Vx_cov[1] = Hxx_Vx[1]; 
  Vx_cov[1] += Hxy_Vy[1]; 
  Vy_cov[1] = Hxy_Vx[1]; 
  Vy_cov[1] += Hyy_Vy[1]; 
 
  // Jnu_i_cov = h_{ij}Jnv^j 
  NVx_cov[0] = Hxx_NVx[0]; 
  NVx_cov[0] += Hxy_NVy[0]; 
  NVy_cov[0] = Hxy_NVx[0]; 
  NVy_cov[0] += Hyy_NVy[0]; 
  NVx_cov[1] = Hxx_NVx[1]; 
  NVx_cov[1] += Hxy_NVy[1]; 
  NVy_cov[1] = Hxy_NVx[1]; 
  NVy_cov[1] += Hyy_NVy[1]; 
 
} 
