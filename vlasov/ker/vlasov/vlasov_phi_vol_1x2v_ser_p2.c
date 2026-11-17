#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_1x2v_ser_p2(const double *w, const double *dxv,
  const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[2] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*f[1]*phi[2])-2.1213203435596424*f[0]*phi[1]); 
  out[4] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(4.242640687119286*f[7])-4.743416490252569*f[0])-2.1213203435596424*f[1]*phi[1]); 
  out[6] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*phi[2]*f[5])-2.1213203435596424*phi[1]*f[3]); 
  out[8] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(10.606601717798211*phi[2]*f[4])-4.743416490252569*phi[1]*f[2]); 
  out[10] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(4.242640687119286*f[13])-4.743416490252569*f[3])-2.1213203435596424*phi[1]*f[5]); 
  out[11] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[7])-4.242640687119286*f[1]*phi[2]); 
  out[12] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(9.48683298050514*f[11])-10.606601717798213*f[2])-4.743416490252569*phi[1]*f[4]); 
  out[14] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(10.606601717798213*phi[2]*f[10])-4.743416490252569*phi[1]*f[6]); 
  out[16] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*phi[2]*f[15])-2.1213203435596424*phi[1]*f[9]); 
  out[17] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[13])-4.242640687119286*phi[2]*f[5]); 
  out[18] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(9.48683298050514*f[17])-10.606601717798211*f[6])-4.743416490252569*phi[1]*f[10]); 
  out[19] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[15])-4.743416490252569*phi[2]*f[9]); 

} 
