#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_E_vol_2x1v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[2]; 
  const double *Ex = &qmem[0]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  out[3] += dv10*jacob_vx_inv*(0.8660254037844386*Ex[3]*f[4]+0.8660254037844386*Ex[2]*f[2]+0.8660254037844386*Ex[1]*f[1]+0.8660254037844386*Ex[0]*f[0]); 
  out[5] += dv10*jacob_vx_inv*(0.8660254037844386*Ex[2]*f[4]+0.8660254037844386*f[2]*Ex[3]+0.8660254037844386*Ex[0]*f[1]+0.8660254037844386*f[0]*Ex[1]); 
  out[6] += dv10*jacob_vx_inv*(0.8660254037844386*Ex[1]*f[4]+0.8660254037844386*f[1]*Ex[3]+0.8660254037844386*Ex[0]*f[2]+0.8660254037844386*f[0]*Ex[2]); 
  out[7] += dv10*jacob_vx_inv*(0.8660254037844386*Ex[0]*f[4]+0.8660254037844386*f[0]*Ex[3]+0.8660254037844386*Ex[1]*f[2]+0.8660254037844386*f[1]*Ex[2]); 

} 
