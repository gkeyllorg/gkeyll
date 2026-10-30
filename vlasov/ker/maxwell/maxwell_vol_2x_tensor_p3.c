#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_vol_2x_tensor_p3(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[16]; 
  const double *ez = &q[32]; 
  const double *bx = &q[48]; 
  const double *by = &q[64]; 
  const double *bz = &q[80]; 
  const double *ph = &q[96]; 
  const double *ps = &q[112]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[16]; 
  double *outEz = &out[32]; 
  double *outBx = &out[48]; 
  double *outBy = &out[64]; 
  double *outBz = &out[80]; 
  double *outPh = &out[96]; 
  double *outPs = &out[112]; 
 
  double dx0 = 2.0/dx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  double dx1 = 2.0/dx[1]; 
  const double *jacob_cy = &jacob_pos[4]; 
  const double jacob_cy_inv = 1.0/jacob_cy[0]; 

  outEx[1] += 1.7320508075688772*ph[0]*c2chi*dx0*jacob_cx_inv; 
  outEx[2] += -(1.7320508075688772*bz[0]*c2*dx1*jacob_cy_inv); 
  outEx[3] += 1.7320508075688772*ph[2]*c2chi*dx0*jacob_cx_inv-1.7320508075688772*bz[1]*c2*dx1*jacob_cy_inv; 
  outEx[4] += 3.872983346207417*ph[1]*c2chi*dx0*jacob_cx_inv; 
  outEx[5] += -(3.872983346207417*bz[2]*c2*dx1*jacob_cy_inv); 
  outEx[6] += 3.872983346207417*ph[3]*c2chi*dx0*jacob_cx_inv-1.7320508075688774*bz[4]*c2*dx1*jacob_cy_inv; 
  outEx[7] += 1.7320508075688774*ph[5]*c2chi*dx0*jacob_cx_inv-3.872983346207417*bz[3]*c2*dx1*jacob_cy_inv; 
  outEx[8] += 5.916079783099617*ph[4]*c2chi*dx0*jacob_cx_inv+2.6457513110645907*ph[0]*c2chi*dx0*jacob_cx_inv; 
  outEx[9] += -(5.916079783099617*bz[5]*c2*dx1*jacob_cy_inv)-2.6457513110645907*bz[0]*c2*dx1*jacob_cy_inv; 
  outEx[10] += 3.872983346207417*ph[7]*c2chi*dx0*jacob_cx_inv-3.872983346207417*bz[6]*c2*dx1*jacob_cy_inv; 
  outEx[11] += -(1.7320508075688772*bz[8]*c2*dx1*jacob_cy_inv)+5.916079783099615*ph[6]*c2chi*dx0*jacob_cx_inv+2.6457513110645907*ph[2]*c2chi*dx0*jacob_cx_inv; 
  outEx[12] += -(5.916079783099615*bz[7]*c2*dx1*jacob_cy_inv)-2.6457513110645907*bz[1]*c2*dx1*jacob_cy_inv+1.7320508075688772*ph[9]*c2chi*dx0*jacob_cx_inv; 
  outEx[13] += -(3.8729833462074166*bz[11]*c2*dx1*jacob_cy_inv)+5.916079783099616*ph[10]*c2chi*dx0*jacob_cx_inv+2.6457513110645907*ph[5]*c2chi*dx0*jacob_cx_inv; 
  outEx[14] += -(5.916079783099616*bz[10]*c2*dx1*jacob_cy_inv)-2.6457513110645907*bz[4]*c2*dx1*jacob_cy_inv+3.8729833462074166*ph[12]*c2chi*dx0*jacob_cx_inv; 
  outEx[15] += -(5.916079783099616*bz[13]*c2*dx1*jacob_cy_inv)-2.6457513110645907*bz[8]*c2*dx1*jacob_cy_inv+5.916079783099616*ph[14]*c2chi*dx0*jacob_cx_inv+2.6457513110645907*ph[9]*c2chi*dx0*jacob_cx_inv; 

  outEy[1] += 1.7320508075688772*bz[0]*c2*dx0*jacob_cx_inv; 
  outEy[2] += 1.7320508075688772*ph[0]*c2chi*dx1*jacob_cy_inv; 
  outEy[3] += 1.7320508075688772*ph[1]*c2chi*dx1*jacob_cy_inv+1.7320508075688772*bz[2]*c2*dx0*jacob_cx_inv; 
  outEy[4] += 3.872983346207417*bz[1]*c2*dx0*jacob_cx_inv; 
  outEy[5] += 3.872983346207417*ph[2]*c2chi*dx1*jacob_cy_inv; 
  outEy[6] += 1.7320508075688774*ph[4]*c2chi*dx1*jacob_cy_inv+3.872983346207417*bz[3]*c2*dx0*jacob_cx_inv; 
  outEy[7] += 3.872983346207417*ph[3]*c2chi*dx1*jacob_cy_inv+1.7320508075688774*bz[5]*c2*dx0*jacob_cx_inv; 
  outEy[8] += 5.916079783099617*bz[4]*c2*dx0*jacob_cx_inv+2.6457513110645907*bz[0]*c2*dx0*jacob_cx_inv; 
  outEy[9] += 5.916079783099617*ph[5]*c2chi*dx1*jacob_cy_inv+2.6457513110645907*ph[0]*c2chi*dx1*jacob_cy_inv; 
  outEy[10] += 3.872983346207417*ph[6]*c2chi*dx1*jacob_cy_inv+3.872983346207417*bz[7]*c2*dx0*jacob_cx_inv; 
  outEy[11] += 1.7320508075688772*ph[8]*c2chi*dx1*jacob_cy_inv+5.916079783099615*bz[6]*c2*dx0*jacob_cx_inv+2.6457513110645907*bz[2]*c2*dx0*jacob_cx_inv; 
  outEy[12] += 5.916079783099615*ph[7]*c2chi*dx1*jacob_cy_inv+2.6457513110645907*ph[1]*c2chi*dx1*jacob_cy_inv+1.7320508075688772*bz[9]*c2*dx0*jacob_cx_inv; 
  outEy[13] += 3.8729833462074166*ph[11]*c2chi*dx1*jacob_cy_inv+5.916079783099616*bz[10]*c2*dx0*jacob_cx_inv+2.6457513110645907*bz[5]*c2*dx0*jacob_cx_inv; 
  outEy[14] += 5.916079783099616*ph[10]*c2chi*dx1*jacob_cy_inv+2.6457513110645907*ph[4]*c2chi*dx1*jacob_cy_inv+3.8729833462074166*bz[12]*c2*dx0*jacob_cx_inv; 
  outEy[15] += 5.916079783099616*ph[13]*c2chi*dx1*jacob_cy_inv+2.6457513110645907*ph[8]*c2chi*dx1*jacob_cy_inv+5.916079783099616*bz[14]*c2*dx0*jacob_cx_inv+2.6457513110645907*bz[9]*c2*dx0*jacob_cx_inv; 

  outEz[1] += -(1.7320508075688772*by[0]*c2*dx0*jacob_cx_inv); 
  outEz[2] += 1.7320508075688772*bx[0]*c2*dx1*jacob_cy_inv; 
  outEz[3] += 1.7320508075688772*bx[1]*c2*dx1*jacob_cy_inv-1.7320508075688772*by[2]*c2*dx0*jacob_cx_inv; 
  outEz[4] += -(3.872983346207417*by[1]*c2*dx0*jacob_cx_inv); 
  outEz[5] += 3.872983346207417*bx[2]*c2*dx1*jacob_cy_inv; 
  outEz[6] += 1.7320508075688774*bx[4]*c2*dx1*jacob_cy_inv-3.872983346207417*by[3]*c2*dx0*jacob_cx_inv; 
  outEz[7] += 3.872983346207417*bx[3]*c2*dx1*jacob_cy_inv-1.7320508075688774*by[5]*c2*dx0*jacob_cx_inv; 
  outEz[8] += -(5.916079783099617*by[4]*c2*dx0*jacob_cx_inv)-2.6457513110645907*by[0]*c2*dx0*jacob_cx_inv; 
  outEz[9] += 5.916079783099617*bx[5]*c2*dx1*jacob_cy_inv+2.6457513110645907*bx[0]*c2*dx1*jacob_cy_inv; 
  outEz[10] += 3.872983346207417*bx[6]*c2*dx1*jacob_cy_inv-3.872983346207417*by[7]*c2*dx0*jacob_cx_inv; 
  outEz[11] += 1.7320508075688772*bx[8]*c2*dx1*jacob_cy_inv-5.916079783099615*by[6]*c2*dx0*jacob_cx_inv-2.6457513110645907*by[2]*c2*dx0*jacob_cx_inv; 
  outEz[12] += 5.916079783099615*bx[7]*c2*dx1*jacob_cy_inv+2.6457513110645907*bx[1]*c2*dx1*jacob_cy_inv-1.7320508075688772*by[9]*c2*dx0*jacob_cx_inv; 
  outEz[13] += 3.8729833462074166*bx[11]*c2*dx1*jacob_cy_inv-5.916079783099616*by[10]*c2*dx0*jacob_cx_inv-2.6457513110645907*by[5]*c2*dx0*jacob_cx_inv; 
  outEz[14] += 5.916079783099616*bx[10]*c2*dx1*jacob_cy_inv+2.6457513110645907*bx[4]*c2*dx1*jacob_cy_inv-3.8729833462074166*by[12]*c2*dx0*jacob_cx_inv; 
  outEz[15] += 5.916079783099616*bx[13]*c2*dx1*jacob_cy_inv+2.6457513110645907*bx[8]*c2*dx1*jacob_cy_inv-5.916079783099616*by[14]*c2*dx0*jacob_cx_inv-2.6457513110645907*by[9]*c2*dx0*jacob_cx_inv; 

  outBx[1] += 1.7320508075688772*ps[0]*dx0*jacob_cx_inv*gamma; 
  outBx[2] += 1.7320508075688772*ez[0]*dx1*jacob_cy_inv; 
  outBx[3] += 1.7320508075688772*ps[2]*dx0*jacob_cx_inv*gamma+1.7320508075688772*ez[1]*dx1*jacob_cy_inv; 
  outBx[4] += 3.872983346207417*ps[1]*dx0*jacob_cx_inv*gamma; 
  outBx[5] += 3.872983346207417*ez[2]*dx1*jacob_cy_inv; 
  outBx[6] += 3.872983346207417*ps[3]*dx0*jacob_cx_inv*gamma+1.7320508075688774*ez[4]*dx1*jacob_cy_inv; 
  outBx[7] += 1.7320508075688774*ps[5]*dx0*jacob_cx_inv*gamma+3.872983346207417*ez[3]*dx1*jacob_cy_inv; 
  outBx[8] += 5.916079783099617*ps[4]*dx0*jacob_cx_inv*gamma+2.6457513110645907*ps[0]*dx0*jacob_cx_inv*gamma; 
  outBx[9] += 5.916079783099617*ez[5]*dx1*jacob_cy_inv+2.6457513110645907*ez[0]*dx1*jacob_cy_inv; 
  outBx[10] += 3.872983346207417*ps[7]*dx0*jacob_cx_inv*gamma+3.872983346207417*ez[6]*dx1*jacob_cy_inv; 
  outBx[11] += 5.916079783099615*ps[6]*dx0*jacob_cx_inv*gamma+2.6457513110645907*ps[2]*dx0*jacob_cx_inv*gamma+1.7320508075688772*ez[8]*dx1*jacob_cy_inv; 
  outBx[12] += 1.7320508075688772*ps[9]*dx0*jacob_cx_inv*gamma+5.916079783099615*ez[7]*dx1*jacob_cy_inv+2.6457513110645907*ez[1]*dx1*jacob_cy_inv; 
  outBx[13] += 5.916079783099616*ps[10]*dx0*jacob_cx_inv*gamma+2.6457513110645907*ps[5]*dx0*jacob_cx_inv*gamma+3.8729833462074166*ez[11]*dx1*jacob_cy_inv; 
  outBx[14] += 3.8729833462074166*ps[12]*dx0*jacob_cx_inv*gamma+5.916079783099616*ez[10]*dx1*jacob_cy_inv+2.6457513110645907*ez[4]*dx1*jacob_cy_inv; 
  outBx[15] += 5.916079783099616*ps[14]*dx0*jacob_cx_inv*gamma+2.6457513110645907*ps[9]*dx0*jacob_cx_inv*gamma+5.916079783099616*ez[13]*dx1*jacob_cy_inv+2.6457513110645907*ez[8]*dx1*jacob_cy_inv; 

  outBy[1] += -(1.7320508075688772*ez[0]*dx0*jacob_cx_inv); 
  outBy[2] += 1.7320508075688772*ps[0]*dx1*jacob_cy_inv*gamma; 
  outBy[3] += 1.7320508075688772*ps[1]*dx1*jacob_cy_inv*gamma-1.7320508075688772*ez[2]*dx0*jacob_cx_inv; 
  outBy[4] += -(3.872983346207417*ez[1]*dx0*jacob_cx_inv); 
  outBy[5] += 3.872983346207417*ps[2]*dx1*jacob_cy_inv*gamma; 
  outBy[6] += 1.7320508075688774*ps[4]*dx1*jacob_cy_inv*gamma-3.872983346207417*ez[3]*dx0*jacob_cx_inv; 
  outBy[7] += 3.872983346207417*ps[3]*dx1*jacob_cy_inv*gamma-1.7320508075688774*ez[5]*dx0*jacob_cx_inv; 
  outBy[8] += -(5.916079783099617*ez[4]*dx0*jacob_cx_inv)-2.6457513110645907*ez[0]*dx0*jacob_cx_inv; 
  outBy[9] += 5.916079783099617*ps[5]*dx1*jacob_cy_inv*gamma+2.6457513110645907*ps[0]*dx1*jacob_cy_inv*gamma; 
  outBy[10] += 3.872983346207417*ps[6]*dx1*jacob_cy_inv*gamma-3.872983346207417*ez[7]*dx0*jacob_cx_inv; 
  outBy[11] += 1.7320508075688772*ps[8]*dx1*jacob_cy_inv*gamma-5.916079783099615*ez[6]*dx0*jacob_cx_inv-2.6457513110645907*ez[2]*dx0*jacob_cx_inv; 
  outBy[12] += 5.916079783099615*ps[7]*dx1*jacob_cy_inv*gamma+2.6457513110645907*ps[1]*dx1*jacob_cy_inv*gamma-1.7320508075688772*ez[9]*dx0*jacob_cx_inv; 
  outBy[13] += 3.8729833462074166*ps[11]*dx1*jacob_cy_inv*gamma-5.916079783099616*ez[10]*dx0*jacob_cx_inv-2.6457513110645907*ez[5]*dx0*jacob_cx_inv; 
  outBy[14] += 5.916079783099616*ps[10]*dx1*jacob_cy_inv*gamma+2.6457513110645907*ps[4]*dx1*jacob_cy_inv*gamma-3.8729833462074166*ez[12]*dx0*jacob_cx_inv; 
  outBy[15] += 5.916079783099616*ps[13]*dx1*jacob_cy_inv*gamma+2.6457513110645907*ps[8]*dx1*jacob_cy_inv*gamma-5.916079783099616*ez[14]*dx0*jacob_cx_inv-2.6457513110645907*ez[9]*dx0*jacob_cx_inv; 

  outBz[1] += 1.7320508075688772*ey[0]*dx0*jacob_cx_inv; 
  outBz[2] += -(1.7320508075688772*ex[0]*dx1*jacob_cy_inv); 
  outBz[3] += 1.7320508075688772*ey[2]*dx0*jacob_cx_inv-1.7320508075688772*ex[1]*dx1*jacob_cy_inv; 
  outBz[4] += 3.872983346207417*ey[1]*dx0*jacob_cx_inv; 
  outBz[5] += -(3.872983346207417*ex[2]*dx1*jacob_cy_inv); 
  outBz[6] += 3.872983346207417*ey[3]*dx0*jacob_cx_inv-1.7320508075688774*ex[4]*dx1*jacob_cy_inv; 
  outBz[7] += 1.7320508075688774*ey[5]*dx0*jacob_cx_inv-3.872983346207417*ex[3]*dx1*jacob_cy_inv; 
  outBz[8] += 5.916079783099617*ey[4]*dx0*jacob_cx_inv+2.6457513110645907*ey[0]*dx0*jacob_cx_inv; 
  outBz[9] += -(5.916079783099617*ex[5]*dx1*jacob_cy_inv)-2.6457513110645907*ex[0]*dx1*jacob_cy_inv; 
  outBz[10] += 3.872983346207417*ey[7]*dx0*jacob_cx_inv-3.872983346207417*ex[6]*dx1*jacob_cy_inv; 
  outBz[11] += -(1.7320508075688772*ex[8]*dx1*jacob_cy_inv)+5.916079783099615*ey[6]*dx0*jacob_cx_inv+2.6457513110645907*ey[2]*dx0*jacob_cx_inv; 
  outBz[12] += -(5.916079783099615*ex[7]*dx1*jacob_cy_inv)-2.6457513110645907*ex[1]*dx1*jacob_cy_inv+1.7320508075688772*ey[9]*dx0*jacob_cx_inv; 
  outBz[13] += -(3.8729833462074166*ex[11]*dx1*jacob_cy_inv)+5.916079783099616*ey[10]*dx0*jacob_cx_inv+2.6457513110645907*ey[5]*dx0*jacob_cx_inv; 
  outBz[14] += -(5.916079783099616*ex[10]*dx1*jacob_cy_inv)-2.6457513110645907*ex[4]*dx1*jacob_cy_inv+3.8729833462074166*ey[12]*dx0*jacob_cx_inv; 
  outBz[15] += -(5.916079783099616*ex[13]*dx1*jacob_cy_inv)-2.6457513110645907*ex[8]*dx1*jacob_cy_inv+5.916079783099616*ey[14]*dx0*jacob_cx_inv+2.6457513110645907*ey[9]*dx0*jacob_cx_inv; 

  outPh[1] += 1.7320508075688772*ex[0]*chi*dx0*jacob_cx_inv; 
  outPh[2] += 1.7320508075688772*ey[0]*chi*dx1*jacob_cy_inv; 
  outPh[3] += 1.7320508075688772*ey[1]*chi*dx1*jacob_cy_inv+1.7320508075688772*ex[2]*chi*dx0*jacob_cx_inv; 
  outPh[4] += 3.872983346207417*ex[1]*chi*dx0*jacob_cx_inv; 
  outPh[5] += 3.872983346207417*ey[2]*chi*dx1*jacob_cy_inv; 
  outPh[6] += 1.7320508075688774*ey[4]*chi*dx1*jacob_cy_inv+3.872983346207417*ex[3]*chi*dx0*jacob_cx_inv; 
  outPh[7] += 3.872983346207417*ey[3]*chi*dx1*jacob_cy_inv+1.7320508075688774*ex[5]*chi*dx0*jacob_cx_inv; 
  outPh[8] += 5.916079783099617*ex[4]*chi*dx0*jacob_cx_inv+2.6457513110645907*ex[0]*chi*dx0*jacob_cx_inv; 
  outPh[9] += 5.916079783099617*ey[5]*chi*dx1*jacob_cy_inv+2.6457513110645907*ey[0]*chi*dx1*jacob_cy_inv; 
  outPh[10] += 3.872983346207417*ey[6]*chi*dx1*jacob_cy_inv+3.872983346207417*ex[7]*chi*dx0*jacob_cx_inv; 
  outPh[11] += 1.7320508075688772*ey[8]*chi*dx1*jacob_cy_inv+5.916079783099615*ex[6]*chi*dx0*jacob_cx_inv+2.6457513110645907*ex[2]*chi*dx0*jacob_cx_inv; 
  outPh[12] += 5.916079783099615*ey[7]*chi*dx1*jacob_cy_inv+2.6457513110645907*ey[1]*chi*dx1*jacob_cy_inv+1.7320508075688772*ex[9]*chi*dx0*jacob_cx_inv; 
  outPh[13] += 3.8729833462074166*ey[11]*chi*dx1*jacob_cy_inv+5.916079783099616*ex[10]*chi*dx0*jacob_cx_inv+2.6457513110645907*ex[5]*chi*dx0*jacob_cx_inv; 
  outPh[14] += 5.916079783099616*ey[10]*chi*dx1*jacob_cy_inv+2.6457513110645907*ey[4]*chi*dx1*jacob_cy_inv+3.8729833462074166*ex[12]*chi*dx0*jacob_cx_inv; 
  outPh[15] += 5.916079783099616*ey[13]*chi*dx1*jacob_cy_inv+2.6457513110645907*ey[8]*chi*dx1*jacob_cy_inv+5.916079783099616*ex[14]*chi*dx0*jacob_cx_inv+2.6457513110645907*ex[9]*chi*dx0*jacob_cx_inv; 

  outPs[1] += 1.7320508075688772*bx[0]*c2gamma*dx0*jacob_cx_inv; 
  outPs[2] += 1.7320508075688772*by[0]*c2gamma*dx1*jacob_cy_inv; 
  outPs[3] += 1.7320508075688772*by[1]*c2gamma*dx1*jacob_cy_inv+1.7320508075688772*bx[2]*c2gamma*dx0*jacob_cx_inv; 
  outPs[4] += 3.872983346207417*bx[1]*c2gamma*dx0*jacob_cx_inv; 
  outPs[5] += 3.872983346207417*by[2]*c2gamma*dx1*jacob_cy_inv; 
  outPs[6] += 1.7320508075688774*by[4]*c2gamma*dx1*jacob_cy_inv+3.872983346207417*bx[3]*c2gamma*dx0*jacob_cx_inv; 
  outPs[7] += 3.872983346207417*by[3]*c2gamma*dx1*jacob_cy_inv+1.7320508075688774*bx[5]*c2gamma*dx0*jacob_cx_inv; 
  outPs[8] += 5.916079783099617*bx[4]*c2gamma*dx0*jacob_cx_inv+2.6457513110645907*bx[0]*c2gamma*dx0*jacob_cx_inv; 
  outPs[9] += 5.916079783099617*by[5]*c2gamma*dx1*jacob_cy_inv+2.6457513110645907*by[0]*c2gamma*dx1*jacob_cy_inv; 
  outPs[10] += 3.872983346207417*by[6]*c2gamma*dx1*jacob_cy_inv+3.872983346207417*bx[7]*c2gamma*dx0*jacob_cx_inv; 
  outPs[11] += 1.7320508075688772*by[8]*c2gamma*dx1*jacob_cy_inv+5.916079783099615*bx[6]*c2gamma*dx0*jacob_cx_inv+2.6457513110645907*bx[2]*c2gamma*dx0*jacob_cx_inv; 
  outPs[12] += 5.916079783099615*by[7]*c2gamma*dx1*jacob_cy_inv+2.6457513110645907*by[1]*c2gamma*dx1*jacob_cy_inv+1.7320508075688772*bx[9]*c2gamma*dx0*jacob_cx_inv; 
  outPs[13] += 3.8729833462074166*by[11]*c2gamma*dx1*jacob_cy_inv+5.916079783099616*bx[10]*c2gamma*dx0*jacob_cx_inv+2.6457513110645907*bx[5]*c2gamma*dx0*jacob_cx_inv; 
  outPs[14] += 5.916079783099616*by[10]*c2gamma*dx1*jacob_cy_inv+2.6457513110645907*by[4]*c2gamma*dx1*jacob_cy_inv+3.8729833462074166*bx[12]*c2gamma*dx0*jacob_cx_inv; 
  outPs[15] += 5.916079783099616*by[13]*c2gamma*dx1*jacob_cy_inv+2.6457513110645907*by[8]*c2gamma*dx1*jacob_cy_inv+5.916079783099616*bx[14]*c2gamma*dx0*jacob_cx_inv+2.6457513110645907*bx[9]*c2gamma*dx0*jacob_cx_inv; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]*jacob_cx_inv; 
  cflFreq += meq->c/dx[1]*jacob_cy_inv; 
  return 7.0*cflFreq; 
} 
