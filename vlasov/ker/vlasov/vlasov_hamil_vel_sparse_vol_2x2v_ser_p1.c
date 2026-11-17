#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_hamil_vel_sparse_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap,
  const double *jacob_pos, const double *jacob_vel, const double *poisson_tensor_conf, const double *hamil,
  const double *f,
  const double *f_no_J, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[1] += 1.5*f[0]*hamil[1]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[5] += 1.5*hamil[1]*f[2]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[6] += 1.5*hamil[1]*f[3]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[8] += 1.5*hamil[1]*f[4]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[11] += 1.5*hamil[1]*f[7]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[12] += 1.5*hamil[1]*f[9]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[13] += 1.5*hamil[1]*f[10]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 
  out[15] += 1.5*hamil[1]*f[14]*dv10*dx10*jacob_cx_inv*jacob_vx_inv; 

  double dx11 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[3]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  const double *jacob_cy = &jacob_pos[2]; 
  const double jacob_cy_inv = 1.0/jacob_cy[0]; 
  out[2] += 1.5*f[0]*hamil[2]*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[5] += 1.5*f[1]*hamil[2]*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[7] += 1.5*hamil[2]*f[3]*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[9] += 1.5*hamil[2]*f[4]*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[11] += 1.5*hamil[2]*f[6]*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[12] += 1.5*hamil[2]*f[8]*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[14] += 1.5*hamil[2]*f[10]*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 
  out[15] += 1.5*hamil[2]*f[13]*dv11*dx11*jacob_cy_inv*jacob_vy_inv; 

} 
