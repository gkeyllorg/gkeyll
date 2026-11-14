#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_Bz_hamil_vel_dense_vol_2x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil,
  const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[2]; 
  double dv11 = 2.0/dxv[3]; 
  const double *Bz = &qmem[20]; 
  out[3] += dv10*dv11*Bz[3]*(0.75*hamil[3]*f[11]+0.75*hamil[2]*f[5]);
  out[3] += dv10*dv11*Bz[2]*(0.75*hamil[3]*f[7]+0.75*f[2]*hamil[2]);
  out[3] += dv10*dv11*Bz[1]*(0.75*hamil[3]*f[6]+0.75*f[1]*hamil[2]);
  out[3] += dv10*dv11*Bz[0]*(0.75*f[3]*hamil[3]+0.75*f[0]*hamil[2]);
  out[4] += dv10*dv11*Bz[3]*(-(0.75*hamil[3]*f[12])-0.75*hamil[1]*f[5]);
  out[4] += dv10*dv11*Bz[2]*(-(0.75*hamil[3]*f[9])-0.75*hamil[1]*f[2]);
  out[4] += dv10*dv11*Bz[1]*(-(0.75*hamil[3]*f[8])-0.75*f[1]*hamil[1]);
  out[4] += dv10*dv11*Bz[0]*(-(0.75*hamil[3]*f[4])-0.75*f[0]*hamil[1]);
  out[6] += dv10*dv11*Bz[2]*(0.75*hamil[3]*f[11]+0.75*hamil[2]*f[5]);
  out[6] += dv10*dv11*Bz[3]*(0.75*hamil[3]*f[7]+0.75*f[2]*hamil[2]);
  out[6] += dv10*dv11*Bz[0]*(0.75*hamil[3]*f[6]+0.75*f[1]*hamil[2]);
  out[6] += dv10*dv11*Bz[1]*(0.75*f[3]*hamil[3]+0.75*f[0]*hamil[2]);
  out[7] += dv10*dv11*Bz[1]*(0.75*hamil[3]*f[11]+0.75*hamil[2]*f[5]);
  out[7] += dv10*dv11*Bz[0]*(0.75*hamil[3]*f[7]+0.75*f[2]*hamil[2]);
  out[7] += dv10*dv11*Bz[3]*(0.75*hamil[3]*f[6]+0.75*f[1]*hamil[2]);
  out[7] += dv10*dv11*Bz[2]*(0.75*f[3]*hamil[3]+0.75*f[0]*hamil[2]);
  out[8] += dv10*dv11*Bz[2]*(-(0.75*hamil[3]*f[12])-0.75*hamil[1]*f[5]);
  out[8] += dv10*dv11*Bz[3]*(-(0.75*hamil[3]*f[9])-0.75*hamil[1]*f[2]);
  out[8] += dv10*dv11*Bz[0]*(-(0.75*hamil[3]*f[8])-0.75*f[1]*hamil[1]);
  out[8] += dv10*dv11*Bz[1]*(-(0.75*hamil[3]*f[4])-0.75*f[0]*hamil[1]);
  out[9] += dv10*dv11*Bz[1]*(-(0.75*hamil[3]*f[12])-0.75*hamil[1]*f[5]);
  out[9] += dv10*dv11*Bz[0]*(-(0.75*hamil[3]*f[9])-0.75*hamil[1]*f[2]);
  out[9] += dv10*dv11*Bz[3]*(-(0.75*hamil[3]*f[8])-0.75*f[1]*hamil[1]);
  out[9] += dv10*dv11*Bz[2]*(-(0.75*hamil[3]*f[4])-0.75*f[0]*hamil[1]);
  out[10] += dv10*dv11*Bz[3]*(0.75*hamil[2]*f[12]-0.75*hamil[1]*f[11]);
  out[10] += dv10*dv11*Bz[2]*(0.75*hamil[2]*f[9]-0.75*hamil[1]*f[7]);
  out[10] += dv10*dv11*Bz[1]*(0.75*hamil[2]*f[8]-0.75*hamil[1]*f[6]);
  out[10] += dv10*dv11*Bz[0]*(0.75*hamil[2]*f[4]-0.75*hamil[1]*f[3]);
  out[11] += dv10*dv11*Bz[0]*(0.75*hamil[3]*f[11]+0.75*hamil[2]*f[5]);
  out[11] += dv10*dv11*Bz[1]*(0.75*hamil[3]*f[7]+0.75*f[2]*hamil[2]);
  out[11] += dv10*dv11*Bz[2]*(0.75*hamil[3]*f[6]+0.75*f[1]*hamil[2]);
  out[11] += dv10*dv11*Bz[3]*(0.75*f[3]*hamil[3]+0.75*f[0]*hamil[2]);
  out[12] += dv10*dv11*Bz[0]*(-(0.75*hamil[3]*f[12])-0.75*hamil[1]*f[5]);
  out[12] += dv10*dv11*Bz[1]*(-(0.75*hamil[3]*f[9])-0.75*hamil[1]*f[2]);
  out[12] += dv10*dv11*Bz[2]*(-(0.75*hamil[3]*f[8])-0.75*f[1]*hamil[1]);
  out[12] += dv10*dv11*Bz[3]*(-(0.75*hamil[3]*f[4])-0.75*f[0]*hamil[1]);
  out[13] += dv10*dv11*Bz[2]*(0.75*hamil[2]*f[12]-0.75*hamil[1]*f[11]);
  out[13] += dv10*dv11*Bz[3]*(0.75*hamil[2]*f[9]-0.75*hamil[1]*f[7]);
  out[13] += dv10*dv11*Bz[0]*(0.75*hamil[2]*f[8]-0.75*hamil[1]*f[6]);
  out[13] += dv10*dv11*Bz[1]*(0.75*hamil[2]*f[4]-0.75*hamil[1]*f[3]);
  out[14] += dv10*dv11*Bz[1]*(0.75*hamil[2]*f[12]-0.75*hamil[1]*f[11]);
  out[14] += dv10*dv11*Bz[0]*(0.75*hamil[2]*f[9]-0.75*hamil[1]*f[7]);
  out[14] += dv10*dv11*Bz[3]*(0.75*hamil[2]*f[8]-0.75*hamil[1]*f[6]);
  out[14] += dv10*dv11*Bz[2]*(0.75*hamil[2]*f[4]-0.75*hamil[1]*f[3]);
  out[15] += dv10*dv11*Bz[0]*(0.75*hamil[2]*f[12]-0.75*hamil[1]*f[11]);
  out[15] += dv10*dv11*Bz[1]*(0.75*hamil[2]*f[9]-0.75*hamil[1]*f[7]);
  out[15] += dv10*dv11*Bz[2]*(0.75*hamil[2]*f[8]-0.75*hamil[1]*f[6]);
  out[15] += dv10*dv11*Bz[3]*(0.75*hamil[2]*f[4]-0.75*hamil[1]*f[3]);

} 
