#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_Bx_hamil_vel_sparse_vol_1x3v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil,
  const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv11 = 2.0/dxv[2]; 
  double dv12 = 2.0/dxv[3]; 
  const double *Bx = &qmem[6]; 
  const double *jacob_vx = &jacob_vel[0]; 
  out[3] += dv11*dv12*0.75*jacob_vx[0]*Bx[1]*f[1]*hamil[3];
  out[3] += dv11*dv12*0.75*Bx[0]*f[0]*jacob_vx[0]*hamil[3];
  out[4] += dv11*dv12*-(0.75*jacob_vx[0]*Bx[1]*f[1]*hamil[2]);
  out[4] += dv11*dv12*-(0.75*Bx[0]*f[0]*jacob_vx[0]*hamil[2]);
  out[6] += dv11*dv12*0.75*Bx[0]*jacob_vx[0]*f[1]*hamil[3];
  out[6] += dv11*dv12*0.75*f[0]*jacob_vx[0]*Bx[1]*hamil[3];
  out[7] += dv11*dv12*0.75*jacob_vx[0]*Bx[1]*hamil[3]*f[5];
  out[7] += dv11*dv12*0.75*Bx[0]*jacob_vx[0]*f[2]*hamil[3];
  out[8] += dv11*dv12*-(0.75*Bx[0]*jacob_vx[0]*f[1]*hamil[2]);
  out[8] += dv11*dv12*-(0.75*f[0]*jacob_vx[0]*Bx[1]*hamil[2]);
  out[9] += dv11*dv12*-(0.75*jacob_vx[0]*Bx[1]*hamil[2]*f[5]);
  out[9] += dv11*dv12*-(0.75*Bx[0]*jacob_vx[0]*f[2]*hamil[2]);
  out[10] += dv11*dv12*jacob_vx[0]*Bx[1]*(0.75*hamil[3]*f[8]-0.75*hamil[2]*f[6]);
  out[10] += dv11*dv12*Bx[0]*jacob_vx[0]*(0.75*hamil[3]*f[4]-0.75*hamil[2]*f[3]);
  out[11] += dv11*dv12*0.75*Bx[0]*jacob_vx[0]*hamil[3]*f[5];
  out[11] += dv11*dv12*0.75*jacob_vx[0]*Bx[1]*f[2]*hamil[3];
  out[12] += dv11*dv12*-(0.75*Bx[0]*jacob_vx[0]*hamil[2]*f[5]);
  out[12] += dv11*dv12*-(0.75*jacob_vx[0]*Bx[1]*f[2]*hamil[2]);
  out[13] += dv11*dv12*Bx[0]*jacob_vx[0]*(0.75*hamil[3]*f[8]-0.75*hamil[2]*f[6]);
  out[13] += dv11*dv12*jacob_vx[0]*Bx[1]*(0.75*hamil[3]*f[4]-0.75*hamil[2]*f[3]);
  out[14] += dv11*dv12*jacob_vx[0]*Bx[1]*(0.75*hamil[3]*f[12]-0.75*hamil[2]*f[11]);
  out[14] += dv11*dv12*Bx[0]*jacob_vx[0]*(0.75*hamil[3]*f[9]-0.75*hamil[2]*f[7]);
  out[15] += dv11*dv12*Bx[0]*jacob_vx[0]*(0.75*hamil[3]*f[12]-0.75*hamil[2]*f[11]);
  out[15] += dv11*dv12*jacob_vx[0]*Bx[1]*(0.75*hamil[3]*f[9]-0.75*hamil[2]*f[7]);

} 
