#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_rad_vol_1x1v_ser_p3(const double *w, const double *dxv, 
  const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double *rad_vx = &rad[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[2] += (1.224744871391589*rad_vx[3]*f[9]+1.224744871391589*rad_vx[2]*f[5]+1.224744871391589*rad_vx[1]*f[2]+1.224744871391589*f[0]*rad_vx[0])*dv10*jacob_vx_inv; 
  out[3] += (1.2247448713915887*rad_vx[3]*f[11]+1.224744871391589*rad_vx[2]*f[7]+1.224744871391589*rad_vx[1]*f[3]+1.224744871391589*rad_vx[0]*f[1])*dv10*jacob_vx_inv; 
  out[5] += (2.405351177211819*rad_vx[2]*f[9]+2.405351177211819*rad_vx[3]*f[5]+2.4494897427831783*rad_vx[1]*f[5]+2.4494897427831783*f[2]*rad_vx[2]+2.7386127875258306*rad_vx[0]*f[2]+2.7386127875258306*f[0]*rad_vx[1])*dv10*jacob_vx_inv; 
  out[6] += (1.224744871391589*rad_vx[1]*f[6]+1.224744871391589*rad_vx[0]*f[4])*dv10*jacob_vx_inv; 
  out[7] += (2.405351177211819*rad_vx[2]*f[11]+2.405351177211819*rad_vx[3]*f[7]+2.4494897427831783*rad_vx[1]*f[7]+2.4494897427831783*rad_vx[2]*f[3]+2.7386127875258306*rad_vx[0]*f[3]+2.7386127875258306*f[1]*rad_vx[1])*dv10*jacob_vx_inv; 
  out[9] += (4.365266951236265*rad_vx[3]*f[9]+3.6742346141747664*rad_vx[1]*f[9]+4.543441112511214*rad_vx[2]*f[5]+4.183300132670378*rad_vx[0]*f[5]+3.6742346141747664*f[2]*rad_vx[3]+4.183300132670378*f[0]*rad_vx[2]+5.612486080160912*rad_vx[1]*f[2]+1.8708286933869707*f[0]*rad_vx[0])*dv10*jacob_vx_inv; 
  out[10] += (1.224744871391589*rad_vx[1]*f[10]+1.2247448713915887*rad_vx[0]*f[8])*dv10*jacob_vx_inv; 
  out[11] += (4.365266951236265*rad_vx[3]*f[11]+3.6742346141747664*rad_vx[1]*f[11]+4.543441112511214*rad_vx[2]*f[7]+4.183300132670377*rad_vx[0]*f[7]+3.674234614174766*f[3]*rad_vx[3]+5.612486080160911*rad_vx[1]*f[3]+4.183300132670378*f[1]*rad_vx[2]+1.8708286933869707*rad_vx[0]*f[1])*dv10*jacob_vx_inv; 

} 
