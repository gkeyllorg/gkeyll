#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_Bz_hamil_phase_vol_2x2v_ser_p1(const double *w, const double *dxv, 
  const double *jacob_vel, const double *hamil,
  const double *qmem, const double *f, double* GKYL_RESTRICT out) 
{ 
  double dv10 = 2.0/dxv[2]; 
  double dv11 = 2.0/dxv[3]; 
  const double *Bz = &qmem[20]; 
  out[3] += dv10*dv11*Bz[0]*(0.375*f[11]*hamil[15]+0.375*f[7]*hamil[14]+0.375*f[6]*hamil[13]+0.375*f[5]*hamil[12]+0.375*f[3]*hamil[10]+0.375*f[2]*hamil[9]+0.375*f[1]*hamil[8]+0.375*f[0]*hamil[4]);
  out[3] += dv10*dv11*Bz[1]*(0.375*f[7]*hamil[15]+0.375*f[11]*hamil[14]+0.375*f[3]*hamil[13]+0.375*f[2]*hamil[12]+0.375*f[6]*hamil[10]+0.375*f[5]*hamil[9]+0.375*f[0]*hamil[8]+0.375*f[1]*hamil[4]);
  out[3] += dv10*dv11*Bz[2]*(0.375*f[6]*hamil[15]+0.375*f[3]*hamil[14]+0.375*f[11]*hamil[13]+0.375*f[1]*hamil[12]+0.375*f[7]*hamil[10]+0.375*f[0]*hamil[9]+0.375*f[5]*hamil[8]+0.375*f[2]*hamil[4]);
  out[3] += dv10*dv11*Bz[3]*(0.375*f[3]*hamil[15]+0.375*f[6]*hamil[14]+0.375*f[7]*hamil[13]+0.375*f[0]*hamil[12]+0.375*hamil[10]*f[11]+0.375*f[1]*hamil[9]+0.375*f[2]*hamil[8]+0.375*hamil[4]*f[5]);
  out[4] += dv10*dv11*Bz[0]*(-(0.375*f[12]*hamil[15])-0.375*f[9]*hamil[14]-0.375*f[8]*hamil[13]-0.375*f[5]*hamil[11]-0.375*f[4]*hamil[10]-0.375*f[2]*hamil[7]-0.375*f[1]*hamil[6]-0.375*f[0]*hamil[3]);
  out[4] += dv10*dv11*Bz[1]*(-(0.375*f[9]*hamil[15])-0.375*f[12]*hamil[14]-0.375*f[4]*hamil[13]-0.375*f[2]*hamil[11]-0.375*f[8]*hamil[10]-0.375*f[5]*hamil[7]-0.375*f[0]*hamil[6]-0.375*f[1]*hamil[3]);
  out[4] += dv10*dv11*Bz[2]*(-(0.375*f[8]*hamil[15])-0.375*f[4]*hamil[14]-0.375*f[12]*hamil[13]-0.375*f[1]*hamil[11]-0.375*f[9]*hamil[10]-0.375*f[0]*hamil[7]-0.375*f[5]*hamil[6]-0.375*f[2]*hamil[3]);
  out[4] += dv10*dv11*Bz[3]*(-(0.375*f[4]*hamil[15])-0.375*f[8]*hamil[14]-0.375*f[9]*hamil[13]-0.375*hamil[10]*f[12]-0.375*f[0]*hamil[11]-0.375*f[1]*hamil[7]-0.375*f[2]*hamil[6]-0.375*hamil[3]*f[5]);
  out[6] += dv10*dv11*Bz[1]*(0.675*f[11]*hamil[15]+0.375*f[7]*hamil[14]+0.675*f[6]*hamil[13]+0.675*f[5]*hamil[12]+0.375*f[3]*hamil[10]+0.375*f[2]*hamil[9]+0.675*f[1]*hamil[8]+0.375*f[0]*hamil[4]);
  out[6] += dv10*dv11*Bz[0]*(0.375*f[7]*hamil[15]+0.375*f[11]*hamil[14]+0.375*f[3]*hamil[13]+0.375*f[2]*hamil[12]+0.375*f[6]*hamil[10]+0.375*f[5]*hamil[9]+0.375*f[0]*hamil[8]+0.375*f[1]*hamil[4]);
  out[6] += dv10*dv11*Bz[3]*(0.675*f[6]*hamil[15]+0.375*f[3]*hamil[14]+0.675*f[11]*hamil[13]+0.675*f[1]*hamil[12]+0.375*f[7]*hamil[10]+0.375*f[0]*hamil[9]+0.675*f[5]*hamil[8]+0.375*f[2]*hamil[4]);
  out[6] += dv10*dv11*Bz[2]*(0.375*f[3]*hamil[15]+0.375*f[6]*hamil[14]+0.375*f[7]*hamil[13]+0.375*f[0]*hamil[12]+0.375*hamil[10]*f[11]+0.375*f[1]*hamil[9]+0.375*f[2]*hamil[8]+0.375*hamil[4]*f[5]);
  out[7] += dv10*dv11*Bz[2]*(0.675*f[11]*hamil[15]+0.675*f[7]*hamil[14]+0.375*f[6]*hamil[13]+0.675*f[5]*hamil[12]+0.375*f[3]*hamil[10]+0.675*f[2]*hamil[9]+0.375*f[1]*hamil[8]+0.375*f[0]*hamil[4]);
  out[7] += dv10*dv11*Bz[3]*(0.675*f[7]*hamil[15]+0.675*f[11]*hamil[14]+0.375*f[3]*hamil[13]+0.675*f[2]*hamil[12]+0.375*f[6]*hamil[10]+0.675*f[5]*hamil[9]+0.375*f[0]*hamil[8]+0.375*f[1]*hamil[4]);
  out[7] += dv10*dv11*Bz[0]*(0.375*f[6]*hamil[15]+0.375*f[3]*hamil[14]+0.375*f[11]*hamil[13]+0.375*f[1]*hamil[12]+0.375*f[7]*hamil[10]+0.375*f[0]*hamil[9]+0.375*f[5]*hamil[8]+0.375*f[2]*hamil[4]);
  out[7] += dv10*dv11*Bz[1]*(0.375*f[3]*hamil[15]+0.375*f[6]*hamil[14]+0.375*f[7]*hamil[13]+0.375*f[0]*hamil[12]+0.375*hamil[10]*f[11]+0.375*f[1]*hamil[9]+0.375*f[2]*hamil[8]+0.375*hamil[4]*f[5]);
  out[8] += dv10*dv11*Bz[1]*(-(0.675*f[12]*hamil[15])-0.375*f[9]*hamil[14]-0.675*f[8]*hamil[13]-0.675*f[5]*hamil[11]-0.375*f[4]*hamil[10]-0.375*f[2]*hamil[7]-0.675*f[1]*hamil[6]-0.375*f[0]*hamil[3]);
  out[8] += dv10*dv11*Bz[0]*(-(0.375*f[9]*hamil[15])-0.375*f[12]*hamil[14]-0.375*f[4]*hamil[13]-0.375*f[2]*hamil[11]-0.375*f[8]*hamil[10]-0.375*f[5]*hamil[7]-0.375*f[0]*hamil[6]-0.375*f[1]*hamil[3]);
  out[8] += dv10*dv11*Bz[3]*(-(0.675*f[8]*hamil[15])-0.375*f[4]*hamil[14]-0.675*f[12]*hamil[13]-0.675*f[1]*hamil[11]-0.375*f[9]*hamil[10]-0.375*f[0]*hamil[7]-0.675*f[5]*hamil[6]-0.375*f[2]*hamil[3]);
  out[8] += dv10*dv11*Bz[2]*(-(0.375*f[4]*hamil[15])-0.375*f[8]*hamil[14]-0.375*f[9]*hamil[13]-0.375*hamil[10]*f[12]-0.375*f[0]*hamil[11]-0.375*f[1]*hamil[7]-0.375*f[2]*hamil[6]-0.375*hamil[3]*f[5]);
  out[9] += dv10*dv11*Bz[2]*(-(0.675*f[12]*hamil[15])-0.675*f[9]*hamil[14]-0.375*f[8]*hamil[13]-0.675*f[5]*hamil[11]-0.375*f[4]*hamil[10]-0.675*f[2]*hamil[7]-0.375*f[1]*hamil[6]-0.375*f[0]*hamil[3]);
  out[9] += dv10*dv11*Bz[3]*(-(0.675*f[9]*hamil[15])-0.675*f[12]*hamil[14]-0.375*f[4]*hamil[13]-0.675*f[2]*hamil[11]-0.375*f[8]*hamil[10]-0.675*f[5]*hamil[7]-0.375*f[0]*hamil[6]-0.375*f[1]*hamil[3]);
  out[9] += dv10*dv11*Bz[0]*(-(0.375*f[8]*hamil[15])-0.375*f[4]*hamil[14]-0.375*f[12]*hamil[13]-0.375*f[1]*hamil[11]-0.375*f[9]*hamil[10]-0.375*f[0]*hamil[7]-0.375*f[5]*hamil[6]-0.375*f[2]*hamil[3]);
  out[9] += dv10*dv11*Bz[1]*(-(0.375*f[4]*hamil[15])-0.375*f[8]*hamil[14]-0.375*f[9]*hamil[13]-0.375*hamil[10]*f[12]-0.375*f[0]*hamil[11]-0.375*f[1]*hamil[7]-0.375*f[2]*hamil[6]-0.375*hamil[3]*f[5]);
  out[10] += dv10*dv11*Bz[0]*(0.375*f[12]*hamil[12]-0.375*f[11]*hamil[11]+0.375*f[9]*hamil[9]+0.375*f[8]*hamil[8]-0.375*f[7]*hamil[7]-0.375*f[6]*hamil[6]+0.375*f[4]*hamil[4]-0.375*f[3]*hamil[3]);
  out[10] += dv10*dv11*Bz[1]*(0.375*f[9]*hamil[12]+0.375*hamil[9]*f[12]-0.375*f[7]*hamil[11]-0.375*hamil[7]*f[11]+0.375*f[4]*hamil[8]+0.375*hamil[4]*f[8]-0.375*f[3]*hamil[6]-0.375*hamil[3]*f[6]);
  out[10] += dv10*dv11*Bz[2]*(0.375*f[8]*hamil[12]+0.375*hamil[8]*f[12]-0.375*f[6]*hamil[11]-0.375*hamil[6]*f[11]+0.375*f[4]*hamil[9]+0.375*hamil[4]*f[9]-0.375*f[3]*hamil[7]-0.375*hamil[3]*f[7]);
  out[10] += dv10*dv11*Bz[3]*(0.375*f[4]*hamil[12]+0.375*hamil[4]*f[12]-0.375*f[3]*hamil[11]-0.375*hamil[3]*f[11]+0.375*f[8]*hamil[9]+0.375*hamil[8]*f[9]-0.375*f[6]*hamil[7]-0.375*hamil[6]*f[7]);
  out[11] += dv10*dv11*Bz[3]*(1.215*f[11]*hamil[15]+0.675*f[7]*hamil[14]+0.675*f[6]*hamil[13]+1.215*f[5]*hamil[12]+0.375*f[3]*hamil[10]+0.675*f[2]*hamil[9]+0.675*f[1]*hamil[8]+0.375*f[0]*hamil[4]);
  out[11] += dv10*dv11*Bz[2]*(0.675*f[7]*hamil[15]+0.675*f[11]*hamil[14]+0.375*f[3]*hamil[13]+0.675*f[2]*hamil[12]+0.375*f[6]*hamil[10]+0.675*f[5]*hamil[9]+0.375*f[0]*hamil[8]+0.375*f[1]*hamil[4]);
  out[11] += dv10*dv11*Bz[1]*(0.675*f[6]*hamil[15]+0.375*f[3]*hamil[14]+0.675*f[11]*hamil[13]+0.675*f[1]*hamil[12]+0.375*f[7]*hamil[10]+0.375*f[0]*hamil[9]+0.675*f[5]*hamil[8]+0.375*f[2]*hamil[4]);
  out[11] += dv10*dv11*Bz[0]*(0.375*f[3]*hamil[15]+0.375*f[6]*hamil[14]+0.375*f[7]*hamil[13]+0.375*f[0]*hamil[12]+0.375*hamil[10]*f[11]+0.375*f[1]*hamil[9]+0.375*f[2]*hamil[8]+0.375*hamil[4]*f[5]);
  out[12] += dv10*dv11*Bz[3]*(-(1.215*f[12]*hamil[15])-0.675*f[9]*hamil[14]-0.675*f[8]*hamil[13]-1.215*f[5]*hamil[11]-0.375*f[4]*hamil[10]-0.675*f[2]*hamil[7]-0.675*f[1]*hamil[6]-0.375*f[0]*hamil[3]);
  out[12] += dv10*dv11*Bz[2]*(-(0.675*f[9]*hamil[15])-0.675*f[12]*hamil[14]-0.375*f[4]*hamil[13]-0.675*f[2]*hamil[11]-0.375*f[8]*hamil[10]-0.675*f[5]*hamil[7]-0.375*f[0]*hamil[6]-0.375*f[1]*hamil[3]);
  out[12] += dv10*dv11*Bz[1]*(-(0.675*f[8]*hamil[15])-0.375*f[4]*hamil[14]-0.675*f[12]*hamil[13]-0.675*f[1]*hamil[11]-0.375*f[9]*hamil[10]-0.375*f[0]*hamil[7]-0.675*f[5]*hamil[6]-0.375*f[2]*hamil[3]);
  out[12] += dv10*dv11*Bz[0]*(-(0.375*f[4]*hamil[15])-0.375*f[8]*hamil[14]-0.375*f[9]*hamil[13]-0.375*hamil[10]*f[12]-0.375*f[0]*hamil[11]-0.375*f[1]*hamil[7]-0.375*f[2]*hamil[6]-0.375*hamil[3]*f[5]);
  out[13] += dv10*dv11*Bz[1]*(0.675*f[12]*hamil[12]-0.675*f[11]*hamil[11]+0.375*f[9]*hamil[9]+0.675*f[8]*hamil[8]-0.375*f[7]*hamil[7]-0.675*f[6]*hamil[6]+0.375*f[4]*hamil[4]-0.375*f[3]*hamil[3]);
  out[13] += dv10*dv11*Bz[0]*(0.375*f[9]*hamil[12]+0.375*hamil[9]*f[12]-0.375*f[7]*hamil[11]-0.375*hamil[7]*f[11]+0.375*f[4]*hamil[8]+0.375*hamil[4]*f[8]-0.375*f[3]*hamil[6]-0.375*hamil[3]*f[6]);
  out[13] += dv10*dv11*Bz[3]*(0.675*f[8]*hamil[12]+0.675*hamil[8]*f[12]-0.675*f[6]*hamil[11]-0.675*hamil[6]*f[11]+0.375*f[4]*hamil[9]+0.375*hamil[4]*f[9]-0.375*f[3]*hamil[7]-0.375*hamil[3]*f[7]);
  out[13] += dv10*dv11*Bz[2]*(0.375*f[4]*hamil[12]+0.375*hamil[4]*f[12]-0.375*f[3]*hamil[11]-0.375*hamil[3]*f[11]+0.375*f[8]*hamil[9]+0.375*hamil[8]*f[9]-0.375*f[6]*hamil[7]-0.375*hamil[6]*f[7]);
  out[14] += dv10*dv11*Bz[2]*(0.675*f[12]*hamil[12]-0.675*f[11]*hamil[11]+0.675*f[9]*hamil[9]+0.375*f[8]*hamil[8]-0.675*f[7]*hamil[7]-0.375*f[6]*hamil[6]+0.375*f[4]*hamil[4]-0.375*f[3]*hamil[3]);
  out[14] += dv10*dv11*Bz[3]*(0.675*f[9]*hamil[12]+0.675*hamil[9]*f[12]-0.675*f[7]*hamil[11]-0.675*hamil[7]*f[11]+0.375*f[4]*hamil[8]+0.375*hamil[4]*f[8]-0.375*f[3]*hamil[6]-0.375*hamil[3]*f[6]);
  out[14] += dv10*dv11*Bz[0]*(0.375*f[8]*hamil[12]+0.375*hamil[8]*f[12]-0.375*f[6]*hamil[11]-0.375*hamil[6]*f[11]+0.375*f[4]*hamil[9]+0.375*hamil[4]*f[9]-0.375*f[3]*hamil[7]-0.375*hamil[3]*f[7]);
  out[14] += dv10*dv11*Bz[1]*(0.375*f[4]*hamil[12]+0.375*hamil[4]*f[12]-0.375*f[3]*hamil[11]-0.375*hamil[3]*f[11]+0.375*f[8]*hamil[9]+0.375*hamil[8]*f[9]-0.375*f[6]*hamil[7]-0.375*hamil[6]*f[7]);
  out[15] += dv10*dv11*Bz[3]*(1.215*f[12]*hamil[12]-1.215*f[11]*hamil[11]+0.675*f[9]*hamil[9]+0.675*f[8]*hamil[8]-0.675*f[7]*hamil[7]-0.675*f[6]*hamil[6]+0.375*f[4]*hamil[4]-0.375*f[3]*hamil[3]);
  out[15] += dv10*dv11*Bz[2]*(0.675*f[9]*hamil[12]+0.675*hamil[9]*f[12]-0.675*f[7]*hamil[11]-0.675*hamil[7]*f[11]+0.375*f[4]*hamil[8]+0.375*hamil[4]*f[8]-0.375*f[3]*hamil[6]-0.375*hamil[3]*f[6]);
  out[15] += dv10*dv11*Bz[1]*(0.675*f[8]*hamil[12]+0.675*hamil[8]*f[12]-0.675*f[6]*hamil[11]-0.675*hamil[6]*f[11]+0.375*f[4]*hamil[9]+0.375*hamil[4]*f[9]-0.375*f[3]*hamil[7]-0.375*hamil[3]*f[7]);
  out[15] += dv10*dv11*Bz[0]*(0.375*f[4]*hamil[12]+0.375*hamil[4]*f[12]-0.375*f[3]*hamil[11]-0.375*hamil[3]*f[11]+0.375*f[8]*hamil[9]+0.375*hamil[8]*f[9]-0.375*f[6]*hamil[7]-0.375*hamil[6]*f[7]);

} 
