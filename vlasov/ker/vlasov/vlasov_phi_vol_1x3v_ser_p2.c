#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_phi_vol_1x3v_ser_p2(const double *w, const double *dxv,
  const double *jacob_pos, const double *jacob_vel, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dx10 = 2.0/dxv[0]; 
  double dv10 = 2.0/dxv[1]; 
  const double *jacob_vx = &jacob_vel[0]; 
  const double jacob_vx_inv = 1.0/jacob_vx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  out[2] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*f[1]*phi[2])-2.1213203435596424*f[0]*phi[1]); 
  out[5] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(4.242640687119286*f[11])-4.743416490252569*f[0])-2.1213203435596424*f[1]*phi[1]); 
  out[7] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*phi[2]*f[6])-2.1213203435596424*phi[1]*f[3]); 
  out[9] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*phi[2]*f[8])-2.1213203435596424*phi[1]*f[4]); 
  out[12] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(10.606601717798211*phi[2]*f[5])-4.743416490252569*phi[1]*f[2]); 
  out[15] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(4.242640687119286*f[21])-4.743416490252569*f[3])-2.1213203435596424*phi[1]*f[6]); 
  out[16] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(4.242640687119286*f[25])-4.743416490252569*f[4])-2.1213203435596424*phi[1]*f[8]); 
  out[18] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*phi[2]*f[17])-2.1213203435596424*phi[1]*f[10]); 
  out[19] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[11])-4.242640687119286*f[1]*phi[2]); 
  out[20] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(9.48683298050514*f[19])-10.606601717798213*f[2])-4.743416490252569*phi[1]*f[5]); 
  out[22] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(10.606601717798213*phi[2]*f[15])-4.743416490252569*phi[1]*f[7]); 
  out[24] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*phi[2]*f[23])-2.1213203435596424*phi[1]*f[13]); 
  out[26] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(10.606601717798213*phi[2]*f[16])-4.743416490252569*phi[1]*f[9]); 
  out[29] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*phi[2]*f[28])-2.1213203435596424*phi[1]*f[14]); 
  out[31] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(4.242640687119286*f[37])-4.743416490252569*f[10])-2.1213203435596424*phi[1]*f[17]); 
  out[32] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[21])-4.242640687119286*phi[2]*f[6]); 
  out[33] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(9.48683298050514*f[32])-10.606601717798211*f[7])-4.743416490252569*phi[1]*f[15]); 
  out[34] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[23])-4.743416490252569*phi[2]*f[13]); 
  out[35] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[25])-4.242640687119286*phi[2]*f[8]); 
  out[36] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(9.48683298050514*f[35])-10.606601717798211*f[9])-4.743416490252569*phi[1]*f[16]); 
  out[38] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(10.606601717798211*phi[2]*f[31])-4.743416490252569*phi[1]*f[18]); 
  out[40] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*phi[2]*f[39])-2.1213203435596424*phi[1]*f[27]); 
  out[41] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[28])-4.743416490252569*phi[2]*f[14]); 
  out[43] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(4.743416490252569*phi[2]*f[42])-2.1213203435596424*phi[1]*f[30]); 
  out[44] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[37])-4.242640687119286*phi[2]*f[17]); 
  out[45] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(phi[2]*(-(9.48683298050514*f[44])-10.606601717798213*f[18])-4.743416490252569*phi[1]*f[31]); 
  out[46] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[39])-4.743416490252569*phi[2]*f[27]); 
  out[47] += dv10*dx10*jacob_cx_inv*jacob_vx_inv*(-(2.1213203435596424*phi[1]*f[42])-4.743416490252569*phi[2]*f[30]); 

} 
