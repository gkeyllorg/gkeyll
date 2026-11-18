#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_E_vol_1x1v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  const double *Ex = &qmem[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[2] += dv10*jacob_vx_inv*(1.224744871391589*Ex[1]*f[1]+1.224744871391589*Ex[0]*f[0]); 
  out[3] += dv10*jacob_vx_inv*(1.224744871391589*Ex[0]*f[1]+1.224744871391589*f[0]*Ex[1]); 

} 
