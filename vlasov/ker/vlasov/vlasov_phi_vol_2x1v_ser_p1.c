#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_2x1v_ser_p1(const double *w, const double *dxv,
  const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[2]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[3] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(1.5*f[2]*phi[3])-1.5*f[0]*phi[1]); 
  out[5] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(1.5*phi[3]*f[4])-1.5*f[1]*phi[1]); 
  out[6] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(1.5*f[0]*phi[3])-1.5*phi[1]*f[2]); 
  out[7] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(1.5*phi[1]*f[4])-1.5*f[1]*phi[3]); 

} 
