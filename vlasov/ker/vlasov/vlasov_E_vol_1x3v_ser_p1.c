#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_E_vol_1x3v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double *Ex = &qmem[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[2] += dv10*jacob_vx_inv*(1.224744871391589*Ex[1]*f[1]+1.224744871391589*Ex[0]*f[0]); 
  out[5] += dv10*jacob_vx_inv*(1.224744871391589*Ex[0]*f[1]+1.224744871391589*f[0]*Ex[1]); 
  out[7] += dv10*jacob_vx_inv*(1.224744871391589*Ex[1]*f[6]+1.224744871391589*Ex[0]*f[3]); 
  out[9] += dv10*jacob_vx_inv*(1.224744871391589*Ex[1]*f[8]+1.224744871391589*Ex[0]*f[4]); 
  out[11] += dv10*jacob_vx_inv*(1.224744871391589*Ex[0]*f[6]+1.224744871391589*Ex[1]*f[3]); 
  out[12] += dv10*jacob_vx_inv*(1.224744871391589*Ex[0]*f[8]+1.224744871391589*Ex[1]*f[4]); 
  out[14] += dv10*jacob_vx_inv*(1.224744871391589*Ex[1]*f[13]+1.224744871391589*Ex[0]*f[10]); 
  out[15] += dv10*jacob_vx_inv*(1.224744871391589*Ex[0]*f[13]+1.224744871391589*Ex[1]*f[10]); 

  double dv11 = 2.0/dxv[2]; 
  const double *Ey = &qmem[2]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[3] += dv11*jacob_vy_inv*(1.224744871391589*Ey[1]*f[1]+1.224744871391589*Ey[0]*f[0]); 
  out[6] += dv11*jacob_vy_inv*(1.224744871391589*Ey[0]*f[1]+1.224744871391589*f[0]*Ey[1]); 
  out[7] += dv11*jacob_vy_inv*(1.224744871391589*Ey[1]*f[5]+1.224744871391589*Ey[0]*f[2]); 
  out[10] += dv11*jacob_vy_inv*(1.224744871391589*Ey[1]*f[8]+1.224744871391589*Ey[0]*f[4]); 
  out[11] += dv11*jacob_vy_inv*(1.224744871391589*Ey[0]*f[5]+1.224744871391589*Ey[1]*f[2]); 
  out[13] += dv11*jacob_vy_inv*(1.224744871391589*Ey[0]*f[8]+1.224744871391589*Ey[1]*f[4]); 
  out[14] += dv11*jacob_vy_inv*(1.224744871391589*Ey[1]*f[12]+1.224744871391589*Ey[0]*f[9]); 
  out[15] += dv11*jacob_vy_inv*(1.224744871391589*Ey[0]*f[12]+1.224744871391589*Ey[1]*f[9]); 

  double dv12 = 2.0/dxv[3]; 
  const double *Ez = &qmem[4]; 
  const double *jacob_vz = &jacob_vel[4]; 
  const double jacob_vz_inv = 1.0/jacob_vz[0]; 
  out[4] += dv12*jacob_vz_inv*(1.224744871391589*Ez[1]*f[1]+1.224744871391589*Ez[0]*f[0]); 
  out[8] += dv12*jacob_vz_inv*(1.224744871391589*Ez[0]*f[1]+1.224744871391589*f[0]*Ez[1]); 
  out[9] += dv12*jacob_vz_inv*(1.224744871391589*Ez[1]*f[5]+1.224744871391589*Ez[0]*f[2]); 
  out[10] += dv12*jacob_vz_inv*(1.224744871391589*Ez[1]*f[6]+1.224744871391589*Ez[0]*f[3]); 
  out[12] += dv12*jacob_vz_inv*(1.224744871391589*Ez[0]*f[5]+1.224744871391589*Ez[1]*f[2]); 
  out[13] += dv12*jacob_vz_inv*(1.224744871391589*Ez[0]*f[6]+1.224744871391589*Ez[1]*f[3]); 
  out[14] += dv12*jacob_vz_inv*(1.224744871391589*Ez[1]*f[11]+1.224744871391589*Ez[0]*f[7]); 
  out[15] += dv12*jacob_vz_inv*(1.224744871391589*Ez[0]*f[11]+1.224744871391589*Ez[1]*f[7]); 

} 
