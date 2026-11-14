#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_Bz_hamil_vel_sparse_vol_1x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil,
  const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[1]; 
  double dv11 = 2.0/dxv[2]; 
  const double *Bz = &qmem[10]; 
  out[2] += dv10*dv11*1.060660171779821*Bz[1]*f[1]*hamil[2];
  out[2] += dv10*dv11*1.060660171779821*Bz[0]*f[0]*hamil[2];
  out[3] += dv10*dv11*-(1.060660171779821*Bz[1]*f[1]*hamil[1]);
  out[3] += dv10*dv11*-(1.060660171779821*Bz[0]*f[0]*hamil[1]);
  out[4] += dv10*dv11*1.060660171779821*Bz[0]*f[1]*hamil[2];
  out[4] += dv10*dv11*1.060660171779821*f[0]*Bz[1]*hamil[2];
  out[5] += dv10*dv11*-(1.060660171779821*Bz[0]*f[1]*hamil[1]);
  out[5] += dv10*dv11*-(1.060660171779821*f[0]*Bz[1]*hamil[1]);
  out[6] += dv10*dv11*Bz[1]*(1.060660171779821*hamil[2]*f[5]-1.060660171779821*hamil[1]*f[4]);
  out[6] += dv10*dv11*Bz[0]*(1.060660171779821*hamil[2]*f[3]-1.060660171779821*hamil[1]*f[2]);
  out[7] += dv10*dv11*Bz[0]*(1.060660171779821*hamil[2]*f[5]-1.060660171779821*hamil[1]*f[4]);
  out[7] += dv10*dv11*Bz[1]*(1.060660171779821*hamil[2]*f[3]-1.060660171779821*hamil[1]*f[2]);

} 
