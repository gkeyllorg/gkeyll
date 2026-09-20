#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_1x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil,
  const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Bz = &qmem[10]; 
  out[2] += dv10*dv11*Bz[0]*(0.75*f[4]*hamil[7]+0.75*f[2]*hamil[6]+0.75*f[1]*hamil[5]+0.75*f[0]*hamil[3]);
  out[2] += dv10*dv11*Bz[1]*(0.75*f[2]*hamil[7]+0.75*f[4]*hamil[6]+0.75*f[0]*hamil[5]+0.75*f[1]*hamil[3]);
  out[3] += dv10*dv11*Bz[0]*(-(0.75*f[5]*hamil[7])-0.75*f[3]*hamil[6]-0.75*f[1]*hamil[4]-0.75*f[0]*hamil[2]);
  out[3] += dv10*dv11*Bz[1]*(-(0.75*f[3]*hamil[7])-0.75*f[5]*hamil[6]-0.75*f[0]*hamil[4]-0.75*f[1]*hamil[2]);
  out[4] += dv10*dv11*Bz[1]*(1.35*f[4]*hamil[7]+0.75*f[2]*hamil[6]+1.35*f[1]*hamil[5]+0.75*f[0]*hamil[3]);
  out[4] += dv10*dv11*Bz[0]*(0.75*f[2]*hamil[7]+0.75*f[4]*hamil[6]+0.75*f[0]*hamil[5]+0.75*f[1]*hamil[3]);
  out[5] += dv10*dv11*Bz[1]*(-(1.35*f[5]*hamil[7])-0.75*f[3]*hamil[6]-1.35*f[1]*hamil[4]-0.75*f[0]*hamil[2]);
  out[5] += dv10*dv11*Bz[0]*(-(0.75*f[3]*hamil[7])-0.75*f[5]*hamil[6]-0.75*f[0]*hamil[4]-0.75*f[1]*hamil[2]);
  out[6] += dv10*dv11*Bz[0]*(0.75*f[5]*hamil[5]-0.75*f[4]*hamil[4]+0.75*f[3]*hamil[3]-0.75*f[2]*hamil[2]);
  out[6] += dv10*dv11*Bz[1]*(0.75*f[3]*hamil[5]+0.75*hamil[3]*f[5]-0.75*f[2]*hamil[4]-0.75*hamil[2]*f[4]);
  out[7] += dv10*dv11*Bz[1]*(1.35*f[5]*hamil[5]-1.35*f[4]*hamil[4]+0.75*f[3]*hamil[3]-0.75*f[2]*hamil[2]);
  out[7] += dv10*dv11*Bz[0]*(0.75*f[3]*hamil[5]+0.75*hamil[3]*f[5]-0.75*f[2]*hamil[4]-0.75*hamil[2]*f[4]);

} 
