#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_rad_vol_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double *rad_vx = &rad[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[2] += (1.224744871391589*rad_vx[2]*f[5]+1.224744871391589*rad_vx[1]*f[2]+1.224744871391589*f[0]*rad_vx[0])*dv10*jacob_vx_inv; 
  out[3] += (1.224744871391589*rad_vx[2]*f[7]+1.224744871391589*rad_vx[1]*f[3]+1.224744871391589*rad_vx[0]*f[1])*dv10*jacob_vx_inv; 
  out[5] += (2.4494897427831783*rad_vx[1]*f[5]+2.4494897427831783*f[2]*rad_vx[2]+2.7386127875258306*rad_vx[0]*f[2]+2.7386127875258306*f[0]*rad_vx[1])*dv10*jacob_vx_inv; 
  out[6] += (1.224744871391589*rad_vx[1]*f[6]+1.224744871391589*rad_vx[0]*f[4])*dv10*jacob_vx_inv; 
  out[7] += (2.4494897427831783*rad_vx[1]*f[7]+2.4494897427831783*rad_vx[2]*f[3]+2.7386127875258306*rad_vx[0]*f[3]+2.7386127875258306*f[1]*rad_vx[1])*dv10*jacob_vx_inv; 

} 
