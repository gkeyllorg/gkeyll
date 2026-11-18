#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_E_vol_1x1v_ser_p2(const double *w, const double *dxv, 
  const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double *Ex = &qmem[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[2] += dv10*jacob_vx_inv*(1.224744871391589*Ex[2]*f[4]+1.224744871391589*Ex[1]*f[1]+1.224744871391589*Ex[0]*f[0]); 
  out[3] += dv10*jacob_vx_inv*(Ex[1]*(1.0954451150103324*f[4]+1.224744871391589*f[0])+1.0954451150103324*f[1]*Ex[2]+1.224744871391589*Ex[0]*f[1]); 
  out[5] += dv10*jacob_vx_inv*(2.7386127875258306*Ex[2]*f[6]+2.7386127875258306*Ex[1]*f[3]+2.7386127875258306*Ex[0]*f[2]); 
  out[6] += dv10*jacob_vx_inv*(1.224744871391589*Ex[0]*f[4]+Ex[2]*(0.7824607964359517*f[4]+1.224744871391589*f[0])+1.0954451150103324*Ex[1]*f[1]); 
  out[7] += dv10*jacob_vx_inv*(Ex[1]*(2.4494897427831783*f[6]+2.7386127875258306*f[2])+2.4494897427831783*Ex[2]*f[3]+2.7386127875258306*Ex[0]*f[3]); 

} 
