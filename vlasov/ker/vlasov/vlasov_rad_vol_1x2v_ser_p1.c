#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_rad_vol_1x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double *rad_vx = &rad[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[2] += (0.8660254037844386*rad_vx[3]*f[6]+0.8660254037844386*rad_vx[2]*f[3]+0.8660254037844386*rad_vx[1]*f[2]+0.8660254037844386*f[0]*rad_vx[0])*dv10*jacob_vx_inv; 
  out[4] += (0.8660254037844386*rad_vx[3]*f[7]+0.8660254037844386*rad_vx[2]*f[5]+0.8660254037844386*rad_vx[1]*f[4]+0.8660254037844386*rad_vx[0]*f[1])*dv10*jacob_vx_inv; 
  out[6] += (0.8660254037844386*rad_vx[1]*f[6]+0.8660254037844386*f[2]*rad_vx[3]+0.8660254037844386*rad_vx[0]*f[3]+0.8660254037844386*f[0]*rad_vx[2])*dv10*jacob_vx_inv; 
  out[7] += (0.8660254037844386*rad_vx[1]*f[7]+0.8660254037844386*rad_vx[0]*f[5]+0.8660254037844386*rad_vx[3]*f[4]+0.8660254037844386*f[1]*rad_vx[2])*dv10*jacob_vx_inv; 

  double dv11 = 2.0/dxv[2]; 
  const double *rad_vy = &rad[4]; 
  const double *jacob_vy = &jacob_vel[2]; 
  const double jacob_vy_inv = 1.0/jacob_vy[0]; 
  out[3] += (0.8660254037844386*rad_vy[3]*f[6]+0.8660254037844386*rad_vy[2]*f[3]+0.8660254037844386*rad_vy[1]*f[2]+0.8660254037844386*f[0]*rad_vy[0])*dv11*jacob_vy_inv; 
  out[5] += (0.8660254037844386*rad_vy[3]*f[7]+0.8660254037844386*rad_vy[2]*f[5]+0.8660254037844386*rad_vy[1]*f[4]+0.8660254037844386*rad_vy[0]*f[1])*dv11*jacob_vy_inv; 
  out[6] += (0.8660254037844386*rad_vy[2]*f[6]+0.8660254037844386*f[3]*rad_vy[3]+0.8660254037844386*rad_vy[0]*f[2]+0.8660254037844386*f[0]*rad_vy[1])*dv11*jacob_vy_inv; 
  out[7] += (0.8660254037844386*rad_vy[2]*f[7]+0.8660254037844386*rad_vy[3]*f[5]+0.8660254037844386*rad_vy[0]*f[4]+0.8660254037844386*f[1]*rad_vy[1])*dv11*jacob_vy_inv; 

} 
