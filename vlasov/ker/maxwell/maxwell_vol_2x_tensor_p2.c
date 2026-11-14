#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_vol_2x_tensor_p2(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[9]; 
  const double *ez = &q[18]; 
  const double *bx = &q[27]; 
  const double *by = &q[36]; 
  const double *bz = &q[45]; 
  const double *ph = &q[54]; 
  const double *ps = &q[63]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[9]; 
  double *outEz = &out[18]; 
  double *outBx = &out[27]; 
  double *outBy = &out[36]; 
  double *outBz = &out[45]; 
  double *outPh = &out[54]; 
  double *outPs = &out[63]; 
 
  double dx0 = 2.0/dx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  double dx1 = 2.0/dx[1]; 
  const double *jacob_cy = &jacob_pos[3]; 
  const double jacob_cy_inv = 1.0/jacob_cy[0]; 

  outEx[1] += 1.7320508075688772*ph[0]*c2chi*dx0*jacob_cx_inv; 
  outEx[2] += -(1.7320508075688772*bz[0]*c2*dx1*jacob_cy_inv); 
  outEx[3] += 1.7320508075688772*ph[2]*c2chi*dx0*jacob_cx_inv-1.7320508075688772*bz[1]*c2*dx1*jacob_cy_inv; 
  outEx[4] += 3.872983346207417*ph[1]*c2chi*dx0*jacob_cx_inv; 
  outEx[5] += -(3.872983346207417*bz[2]*c2*dx1*jacob_cy_inv); 
  outEx[6] += 3.872983346207417*ph[3]*c2chi*dx0*jacob_cx_inv-1.7320508075688774*bz[4]*c2*dx1*jacob_cy_inv; 
  outEx[7] += 1.7320508075688774*ph[5]*c2chi*dx0*jacob_cx_inv-3.872983346207417*bz[3]*c2*dx1*jacob_cy_inv; 
  outEx[8] += 3.872983346207417*ph[7]*c2chi*dx0*jacob_cx_inv-3.872983346207417*bz[6]*c2*dx1*jacob_cy_inv; 

  outEy[1] += 1.7320508075688772*bz[0]*c2*dx0*jacob_cx_inv; 
  outEy[2] += 1.7320508075688772*ph[0]*c2chi*dx1*jacob_cy_inv; 
  outEy[3] += 1.7320508075688772*ph[1]*c2chi*dx1*jacob_cy_inv+1.7320508075688772*bz[2]*c2*dx0*jacob_cx_inv; 
  outEy[4] += 3.872983346207417*bz[1]*c2*dx0*jacob_cx_inv; 
  outEy[5] += 3.872983346207417*ph[2]*c2chi*dx1*jacob_cy_inv; 
  outEy[6] += 1.7320508075688774*ph[4]*c2chi*dx1*jacob_cy_inv+3.872983346207417*bz[3]*c2*dx0*jacob_cx_inv; 
  outEy[7] += 3.872983346207417*ph[3]*c2chi*dx1*jacob_cy_inv+1.7320508075688774*bz[5]*c2*dx0*jacob_cx_inv; 
  outEy[8] += 3.872983346207417*ph[6]*c2chi*dx1*jacob_cy_inv+3.872983346207417*bz[7]*c2*dx0*jacob_cx_inv; 

  outEz[1] += -(1.7320508075688772*by[0]*c2*dx0*jacob_cx_inv); 
  outEz[2] += 1.7320508075688772*bx[0]*c2*dx1*jacob_cy_inv; 
  outEz[3] += 1.7320508075688772*bx[1]*c2*dx1*jacob_cy_inv-1.7320508075688772*by[2]*c2*dx0*jacob_cx_inv; 
  outEz[4] += -(3.872983346207417*by[1]*c2*dx0*jacob_cx_inv); 
  outEz[5] += 3.872983346207417*bx[2]*c2*dx1*jacob_cy_inv; 
  outEz[6] += 1.7320508075688774*bx[4]*c2*dx1*jacob_cy_inv-3.872983346207417*by[3]*c2*dx0*jacob_cx_inv; 
  outEz[7] += 3.872983346207417*bx[3]*c2*dx1*jacob_cy_inv-1.7320508075688774*by[5]*c2*dx0*jacob_cx_inv; 
  outEz[8] += 3.872983346207417*bx[6]*c2*dx1*jacob_cy_inv-3.872983346207417*by[7]*c2*dx0*jacob_cx_inv; 

  outBx[1] += 1.7320508075688772*ps[0]*dx0*jacob_cx_inv*gamma; 
  outBx[2] += 1.7320508075688772*ez[0]*dx1*jacob_cy_inv; 
  outBx[3] += 1.7320508075688772*ps[2]*dx0*jacob_cx_inv*gamma+1.7320508075688772*ez[1]*dx1*jacob_cy_inv; 
  outBx[4] += 3.872983346207417*ps[1]*dx0*jacob_cx_inv*gamma; 
  outBx[5] += 3.872983346207417*ez[2]*dx1*jacob_cy_inv; 
  outBx[6] += 3.872983346207417*ps[3]*dx0*jacob_cx_inv*gamma+1.7320508075688774*ez[4]*dx1*jacob_cy_inv; 
  outBx[7] += 1.7320508075688774*ps[5]*dx0*jacob_cx_inv*gamma+3.872983346207417*ez[3]*dx1*jacob_cy_inv; 
  outBx[8] += 3.872983346207417*ps[7]*dx0*jacob_cx_inv*gamma+3.872983346207417*ez[6]*dx1*jacob_cy_inv; 

  outBy[1] += -(1.7320508075688772*ez[0]*dx0*jacob_cx_inv); 
  outBy[2] += 1.7320508075688772*ps[0]*dx1*jacob_cy_inv*gamma; 
  outBy[3] += 1.7320508075688772*ps[1]*dx1*jacob_cy_inv*gamma-1.7320508075688772*ez[2]*dx0*jacob_cx_inv; 
  outBy[4] += -(3.872983346207417*ez[1]*dx0*jacob_cx_inv); 
  outBy[5] += 3.872983346207417*ps[2]*dx1*jacob_cy_inv*gamma; 
  outBy[6] += 1.7320508075688774*ps[4]*dx1*jacob_cy_inv*gamma-3.872983346207417*ez[3]*dx0*jacob_cx_inv; 
  outBy[7] += 3.872983346207417*ps[3]*dx1*jacob_cy_inv*gamma-1.7320508075688774*ez[5]*dx0*jacob_cx_inv; 
  outBy[8] += 3.872983346207417*ps[6]*dx1*jacob_cy_inv*gamma-3.872983346207417*ez[7]*dx0*jacob_cx_inv; 

  outBz[1] += 1.7320508075688772*ey[0]*dx0*jacob_cx_inv; 
  outBz[2] += -(1.7320508075688772*ex[0]*dx1*jacob_cy_inv); 
  outBz[3] += 1.7320508075688772*ey[2]*dx0*jacob_cx_inv-1.7320508075688772*ex[1]*dx1*jacob_cy_inv; 
  outBz[4] += 3.872983346207417*ey[1]*dx0*jacob_cx_inv; 
  outBz[5] += -(3.872983346207417*ex[2]*dx1*jacob_cy_inv); 
  outBz[6] += 3.872983346207417*ey[3]*dx0*jacob_cx_inv-1.7320508075688774*ex[4]*dx1*jacob_cy_inv; 
  outBz[7] += 1.7320508075688774*ey[5]*dx0*jacob_cx_inv-3.872983346207417*ex[3]*dx1*jacob_cy_inv; 
  outBz[8] += 3.872983346207417*ey[7]*dx0*jacob_cx_inv-3.872983346207417*ex[6]*dx1*jacob_cy_inv; 

  outPh[1] += 1.7320508075688772*ex[0]*chi*dx0*jacob_cx_inv; 
  outPh[2] += 1.7320508075688772*ey[0]*chi*dx1*jacob_cy_inv; 
  outPh[3] += 1.7320508075688772*ey[1]*chi*dx1*jacob_cy_inv+1.7320508075688772*ex[2]*chi*dx0*jacob_cx_inv; 
  outPh[4] += 3.872983346207417*ex[1]*chi*dx0*jacob_cx_inv; 
  outPh[5] += 3.872983346207417*ey[2]*chi*dx1*jacob_cy_inv; 
  outPh[6] += 1.7320508075688774*ey[4]*chi*dx1*jacob_cy_inv+3.872983346207417*ex[3]*chi*dx0*jacob_cx_inv; 
  outPh[7] += 3.872983346207417*ey[3]*chi*dx1*jacob_cy_inv+1.7320508075688774*ex[5]*chi*dx0*jacob_cx_inv; 
  outPh[8] += 3.872983346207417*ey[6]*chi*dx1*jacob_cy_inv+3.872983346207417*ex[7]*chi*dx0*jacob_cx_inv; 

  outPs[1] += 1.7320508075688772*bx[0]*c2gamma*dx0*jacob_cx_inv; 
  outPs[2] += 1.7320508075688772*by[0]*c2gamma*dx1*jacob_cy_inv; 
  outPs[3] += 1.7320508075688772*by[1]*c2gamma*dx1*jacob_cy_inv+1.7320508075688772*bx[2]*c2gamma*dx0*jacob_cx_inv; 
  outPs[4] += 3.872983346207417*bx[1]*c2gamma*dx0*jacob_cx_inv; 
  outPs[5] += 3.872983346207417*by[2]*c2gamma*dx1*jacob_cy_inv; 
  outPs[6] += 1.7320508075688774*by[4]*c2gamma*dx1*jacob_cy_inv+3.872983346207417*bx[3]*c2gamma*dx0*jacob_cx_inv; 
  outPs[7] += 3.872983346207417*by[3]*c2gamma*dx1*jacob_cy_inv+1.7320508075688774*bx[5]*c2gamma*dx0*jacob_cx_inv; 
  outPs[8] += 3.872983346207417*by[6]*c2gamma*dx1*jacob_cy_inv+3.872983346207417*bx[7]*c2gamma*dx0*jacob_cx_inv; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]*jacob_cx_inv; 
  cflFreq += meq->c/dx[1]*jacob_cy_inv; 
  return 5.0*cflFreq; 
} 
