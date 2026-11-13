#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_1x3v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil,
  const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Bz = &qmem[10]; 
  const double *jacob_vz = &jacob_vel[4]; 
  out[2] += dv10*dv11*0.75*jacob_vz[0]*Bz[1]*f[1]*hamil[2];
  out[2] += dv10*dv11*0.75*Bz[0]*f[0]*jacob_vz[0]*hamil[2];
  out[3] += dv10*dv11*-(0.75*jacob_vz[0]*Bz[1]*f[1]*hamil[1]);
  out[3] += dv10*dv11*-(0.75*Bz[0]*f[0]*jacob_vz[0]*hamil[1]);
  out[5] += dv10*dv11*0.75*Bz[0]*jacob_vz[0]*f[1]*hamil[2];
  out[5] += dv10*dv11*0.75*f[0]*jacob_vz[0]*Bz[1]*hamil[2];
  out[6] += dv10*dv11*-(0.75*Bz[0]*jacob_vz[0]*f[1]*hamil[1]);
  out[6] += dv10*dv11*-(0.75*f[0]*jacob_vz[0]*Bz[1]*hamil[1]);
  out[7] += dv10*dv11*jacob_vz[0]*Bz[1]*(0.75*hamil[2]*f[6]-0.75*hamil[1]*f[5]);
  out[7] += dv10*dv11*Bz[0]*jacob_vz[0]*(0.75*hamil[2]*f[3]-0.75*hamil[1]*f[2]);
  out[9] += dv10*dv11*0.75*jacob_vz[0]*Bz[1]*hamil[2]*f[8];
  out[9] += dv10*dv11*0.75*Bz[0]*jacob_vz[0]*hamil[2]*f[4];
  out[10] += dv10*dv11*-(0.75*jacob_vz[0]*Bz[1]*hamil[1]*f[8]);
  out[10] += dv10*dv11*-(0.75*Bz[0]*jacob_vz[0]*hamil[1]*f[4]);
  out[11] += dv10*dv11*Bz[0]*jacob_vz[0]*(0.75*hamil[2]*f[6]-0.75*hamil[1]*f[5]);
  out[11] += dv10*dv11*jacob_vz[0]*Bz[1]*(0.75*hamil[2]*f[3]-0.75*hamil[1]*f[2]);
  out[12] += dv10*dv11*0.75*Bz[0]*jacob_vz[0]*hamil[2]*f[8];
  out[12] += dv10*dv11*0.75*jacob_vz[0]*Bz[1]*hamil[2]*f[4];
  out[13] += dv10*dv11*-(0.75*Bz[0]*jacob_vz[0]*hamil[1]*f[8]);
  out[13] += dv10*dv11*-(0.75*jacob_vz[0]*Bz[1]*hamil[1]*f[4]);
  out[14] += dv10*dv11*jacob_vz[0]*Bz[1]*(0.75*hamil[2]*f[13]-0.75*hamil[1]*f[12]);
  out[14] += dv10*dv11*Bz[0]*jacob_vz[0]*(0.75*hamil[2]*f[10]-0.75*hamil[1]*f[9]);
  out[15] += dv10*dv11*Bz[0]*jacob_vz[0]*(0.75*hamil[2]*f[13]-0.75*hamil[1]*f[12]);
  out[15] += dv10*dv11*jacob_vz[0]*Bz[1]*(0.75*hamil[2]*f[10]-0.75*hamil[1]*f[9]);

} 
