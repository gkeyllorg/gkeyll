#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_By_hamil_vel_dense_vol_1x3v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil,
  const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  double dv12 = 2.0/dxv[3]; 
  const double *By = &qmem[8]; 
  const double *jacob_vy = &jacob_vel[2]; 
  out[2] += dv10*dv12*jacob_vy[0]*By[1]*(-(0.75*hamil[7]*f[11])-0.75*f[6]*hamil[6]-0.75*f[5]*hamil[5]-0.75*f[1]*hamil[3]);
  out[2] += dv10*dv12*By[0]*jacob_vy[0]*(-(0.75*f[7]*hamil[7])-0.75*f[3]*hamil[6]-0.75*f[2]*hamil[5]-0.75*f[0]*hamil[3]);
  out[4] += dv10*dv12*jacob_vy[0]*By[1]*(0.75*hamil[7]*f[13]+0.75*hamil[5]*f[8]+0.75*hamil[4]*f[6]+0.75*f[1]*hamil[1]);
  out[4] += dv10*dv12*By[0]*jacob_vy[0]*(0.75*hamil[7]*f[10]+0.75*f[4]*hamil[5]+0.75*f[3]*hamil[4]+0.75*f[0]*hamil[1]);
  out[5] += dv10*dv12*By[0]*jacob_vy[0]*(-(0.75*hamil[7]*f[11])-0.75*f[6]*hamil[6]-0.75*f[5]*hamil[5]-0.75*f[1]*hamil[3]);
  out[5] += dv10*dv12*jacob_vy[0]*By[1]*(-(0.75*f[7]*hamil[7])-0.75*f[3]*hamil[6]-0.75*f[2]*hamil[5]-0.75*f[0]*hamil[3]);
  out[7] += dv10*dv12*jacob_vy[0]*By[1]*(-(0.75*hamil[5]*f[11])-0.75*f[5]*hamil[7]-0.75*f[1]*hamil[6]-0.75*hamil[3]*f[6]);
  out[7] += dv10*dv12*By[0]*jacob_vy[0]*(-(0.75*f[2]*hamil[7])-0.75*hamil[5]*f[7]-0.75*f[0]*hamil[6]-0.75*f[3]*hamil[3]);
  out[8] += dv10*dv12*By[0]*jacob_vy[0]*(0.75*hamil[7]*f[13]+0.75*hamil[5]*f[8]+0.75*hamil[4]*f[6]+0.75*f[1]*hamil[1]);
  out[8] += dv10*dv12*jacob_vy[0]*By[1]*(0.75*hamil[7]*f[10]+0.75*f[4]*hamil[5]+0.75*f[3]*hamil[4]+0.75*f[0]*hamil[1]);
  out[9] += dv10*dv12*jacob_vy[0]*By[1]*(-(0.75*hamil[6]*f[13])+0.75*hamil[4]*f[11]-0.75*hamil[3]*f[8]+0.75*hamil[1]*f[5]);
  out[9] += dv10*dv12*By[0]*jacob_vy[0]*(-(0.75*hamil[6]*f[10])+0.75*hamil[4]*f[7]-0.75*hamil[3]*f[4]+0.75*hamil[1]*f[2]);
  out[10] += dv10*dv12*jacob_vy[0]*By[1]*(0.75*hamil[5]*f[13]+0.75*hamil[7]*f[8]+0.75*hamil[1]*f[6]+0.75*f[1]*hamil[4]);
  out[10] += dv10*dv12*By[0]*jacob_vy[0]*(0.75*hamil[5]*f[10]+0.75*f[4]*hamil[7]+0.75*f[0]*hamil[4]+0.75*hamil[1]*f[3]);
  out[11] += dv10*dv12*By[0]*jacob_vy[0]*(-(0.75*hamil[5]*f[11])-0.75*f[5]*hamil[7]-0.75*f[1]*hamil[6]-0.75*hamil[3]*f[6]);
  out[11] += dv10*dv12*jacob_vy[0]*By[1]*(-(0.75*f[2]*hamil[7])-0.75*hamil[5]*f[7]-0.75*f[0]*hamil[6]-0.75*f[3]*hamil[3]);
  out[12] += dv10*dv12*By[0]*jacob_vy[0]*(-(0.75*hamil[6]*f[13])+0.75*hamil[4]*f[11]-0.75*hamil[3]*f[8]+0.75*hamil[1]*f[5]);
  out[12] += dv10*dv12*jacob_vy[0]*By[1]*(-(0.75*hamil[6]*f[10])+0.75*hamil[4]*f[7]-0.75*hamil[3]*f[4]+0.75*hamil[1]*f[2]);
  out[13] += dv10*dv12*By[0]*jacob_vy[0]*(0.75*hamil[5]*f[13]+0.75*hamil[7]*f[8]+0.75*hamil[1]*f[6]+0.75*f[1]*hamil[4]);
  out[13] += dv10*dv12*jacob_vy[0]*By[1]*(0.75*hamil[5]*f[10]+0.75*f[4]*hamil[7]+0.75*f[0]*hamil[4]+0.75*hamil[1]*f[3]);
  out[14] += dv10*dv12*jacob_vy[0]*By[1]*(-(0.75*hamil[3]*f[13])+0.75*hamil[1]*f[11]-0.75*hamil[6]*f[8]+0.75*hamil[4]*f[5]);
  out[14] += dv10*dv12*By[0]*jacob_vy[0]*(-(0.75*hamil[3]*f[10])+0.75*hamil[1]*f[7]-0.75*f[4]*hamil[6]+0.75*f[2]*hamil[4]);
  out[15] += dv10*dv12*By[0]*jacob_vy[0]*(-(0.75*hamil[3]*f[13])+0.75*hamil[1]*f[11]-0.75*hamil[6]*f[8]+0.75*hamil[4]*f[5]);
  out[15] += dv10*dv12*jacob_vy[0]*By[1]*(-(0.75*hamil[3]*f[10])+0.75*hamil[1]*f[7]-0.75*f[4]*hamil[6]+0.75*f[2]*hamil[4]);

} 
