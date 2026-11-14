#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_rad_vol_2x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *rad, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[2]; 
  const double *rad_vx = &rad[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[3] += (1.224744871391589*rad_vx[2]*f[9]+1.224744871391589*rad_vx[1]*f[3]+1.224744871391589*f[0]*rad_vx[0])*dv10*jacob_vx_inv; 
  out[5] += (1.224744871391589*rad_vx[2]*f[15]+1.224744871391589*rad_vx[1]*f[5]+1.224744871391589*rad_vx[0]*f[1])*dv10*jacob_vx_inv; 
  out[6] += (1.224744871391589*rad_vx[2]*f[16]+1.224744871391589*rad_vx[1]*f[6]+1.224744871391589*rad_vx[0]*f[2])*dv10*jacob_vx_inv; 
  out[9] += (2.4494897427831783*rad_vx[1]*f[9]+2.4494897427831783*rad_vx[2]*f[3]+2.7386127875258306*rad_vx[0]*f[3]+2.7386127875258306*f[0]*rad_vx[1])*dv10*jacob_vx_inv; 
  out[10] += (1.224744871391589*rad_vx[2]*f[19]+1.224744871391589*rad_vx[1]*f[10]+1.224744871391589*rad_vx[0]*f[4])*dv10*jacob_vx_inv; 
  out[13] += (1.224744871391589*rad_vx[1]*f[13]+1.224744871391589*rad_vx[0]*f[7])*dv10*jacob_vx_inv; 
  out[14] += (1.224744871391589*rad_vx[1]*f[14]+1.224744871391589*rad_vx[0]*f[8])*dv10*jacob_vx_inv; 
  out[15] += (2.4494897427831783*rad_vx[1]*f[15]+2.4494897427831783*rad_vx[2]*f[5]+2.7386127875258306*rad_vx[0]*f[5]+2.7386127875258306*f[1]*rad_vx[1])*dv10*jacob_vx_inv; 
  out[16] += (2.4494897427831783*rad_vx[1]*f[16]+2.4494897427831783*rad_vx[2]*f[6]+2.7386127875258306*rad_vx[0]*f[6]+2.7386127875258306*rad_vx[1]*f[2])*dv10*jacob_vx_inv; 
  out[17] += (1.224744871391589*rad_vx[1]*f[17]+1.224744871391589*rad_vx[0]*f[11])*dv10*jacob_vx_inv; 
  out[18] += (1.224744871391589*rad_vx[1]*f[18]+1.224744871391589*rad_vx[0]*f[12])*dv10*jacob_vx_inv; 
  out[19] += (2.4494897427831783*rad_vx[1]*f[19]+2.4494897427831783*rad_vx[2]*f[10]+2.7386127875258306*rad_vx[0]*f[10]+2.7386127875258306*rad_vx[1]*f[4])*dv10*jacob_vx_inv; 

} 
