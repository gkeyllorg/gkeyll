#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_1x1v_ser_p2(const double *w, const double *dxv,
  const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[2] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*f[1]*phi[2])-2.1213203435596424*f[0]*phi[1]); 
  out[3] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(4.242640687119286*f[4])-4.743416490252569*f[0])-2.1213203435596424*f[1]*phi[1]); 
  out[5] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(10.606601717798211*phi[2]*f[3])-4.743416490252569*phi[1]*f[2]); 
  out[6] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[4])-4.242640687119286*f[1]*phi[2]); 
  out[7] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(9.48683298050514*f[6])-10.606601717798213*f[2])-4.743416490252569*phi[1]*f[3]); 

} 
