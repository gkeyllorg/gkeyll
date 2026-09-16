#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_By_hamil_vel_sparse_vol_1x3v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil,
  const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  double dv12 = 2.0/dxv[3]; 
  const double *By = &qmem[8]; 
  const double *jacob_vy = &jacob_vel[2]; 
  out[2] += dv10*dv12*-(0.75*jacob_vy[0]*By[1]*f[1]*hamil[3]);
  out[2] += dv10*dv12*-(0.75*By[0]*f[0]*jacob_vy[0]*hamil[3]);
  out[4] += dv10*dv12*0.75*jacob_vy[0]*By[1]*f[1]*hamil[1];
  out[4] += dv10*dv12*0.75*By[0]*f[0]*jacob_vy[0]*hamil[1];
  out[5] += dv10*dv12*-(0.75*By[0]*jacob_vy[0]*f[1]*hamil[3]);
  out[5] += dv10*dv12*-(0.75*f[0]*jacob_vy[0]*By[1]*hamil[3]);
  out[7] += dv10*dv12*-(0.75*jacob_vy[0]*By[1]*hamil[3]*f[6]);
  out[7] += dv10*dv12*-(0.75*By[0]*jacob_vy[0]*f[3]*hamil[3]);
  out[8] += dv10*dv12*0.75*By[0]*jacob_vy[0]*f[1]*hamil[1];
  out[8] += dv10*dv12*0.75*f[0]*jacob_vy[0]*By[1]*hamil[1];
  out[9] += dv10*dv12*jacob_vy[0]*By[1]*(0.75*hamil[1]*f[5]-0.75*hamil[3]*f[8]);
  out[9] += dv10*dv12*By[0]*jacob_vy[0]*(0.75*hamil[1]*f[2]-0.75*hamil[3]*f[4]);
  out[10] += dv10*dv12*0.75*jacob_vy[0]*By[1]*hamil[1]*f[6];
  out[10] += dv10*dv12*0.75*By[0]*jacob_vy[0]*hamil[1]*f[3];
  out[11] += dv10*dv12*-(0.75*By[0]*jacob_vy[0]*hamil[3]*f[6]);
  out[11] += dv10*dv12*-(0.75*jacob_vy[0]*By[1]*f[3]*hamil[3]);
  out[12] += dv10*dv12*By[0]*jacob_vy[0]*(0.75*hamil[1]*f[5]-0.75*hamil[3]*f[8]);
  out[12] += dv10*dv12*jacob_vy[0]*By[1]*(0.75*hamil[1]*f[2]-0.75*hamil[3]*f[4]);
  out[13] += dv10*dv12*0.75*By[0]*jacob_vy[0]*hamil[1]*f[6];
  out[13] += dv10*dv12*0.75*jacob_vy[0]*By[1]*hamil[1]*f[3];
  out[14] += dv10*dv12*jacob_vy[0]*By[1]*(0.75*hamil[1]*f[11]-0.75*hamil[3]*f[13]);
  out[14] += dv10*dv12*By[0]*jacob_vy[0]*(0.75*hamil[1]*f[7]-0.75*hamil[3]*f[10]);
  out[15] += dv10*dv12*By[0]*jacob_vy[0]*(0.75*hamil[1]*f[11]-0.75*hamil[3]*f[13]);
  out[15] += dv10*dv12*jacob_vy[0]*By[1]*(0.75*hamil[1]*f[7]-0.75*hamil[3]*f[10]);

} 
