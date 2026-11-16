#include <gkyl_maxwell_kernels.h> 
GKYL_CU_DH double maxwell_vol_2x_ser_p1(const gkyl_maxwell_inp *meq, const double *w, const double *dx, const double *jacob_pos, const double *q, double* GKYL_RESTRICT out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[4]; 
  const double *ez = &q[8]; 
  const double *bx = &q[12]; 
  const double *by = &q[16]; 
  const double *bz = &q[20]; 
  const double *ph = &q[24]; 
  const double *ps = &q[28]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[4]; 
  double *outEz = &out[8]; 
  double *outBx = &out[12]; 
  double *outBy = &out[16]; 
  double *outBz = &out[20]; 
  double *outPh = &out[24]; 
  double *outPs = &out[28]; 
 
  double dx0 = 2.0/dx[0]; 
  const double *jacob_cx = &jacob_pos[0]; 
  const double jacob_cx_inv = 1.0/jacob_cx[0]; 
  double dx1 = 2.0/dx[1]; 
  const double *jacob_cy = &jacob_pos[2]; 
  const double jacob_cy_inv = 1.0/jacob_cy[0]; 

  outEx[1] += 1.7320508075688772*ph[0]*c2chi*dx0*jacob_cx_inv; 
  outEx[2] += -(1.7320508075688772*bz[0]*c2*dx1*jacob_cy_inv); 
  outEx[3] += 1.7320508075688772*ph[2]*c2chi*dx0*jacob_cx_inv-1.7320508075688772*bz[1]*c2*dx1*jacob_cy_inv; 

  outEy[1] += 1.7320508075688772*bz[0]*c2*dx0*jacob_cx_inv; 
  outEy[2] += 1.7320508075688772*ph[0]*c2chi*dx1*jacob_cy_inv; 
  outEy[3] += 1.7320508075688772*ph[1]*c2chi*dx1*jacob_cy_inv+1.7320508075688772*bz[2]*c2*dx0*jacob_cx_inv; 

  outEz[1] += -(1.7320508075688772*by[0]*c2*dx0*jacob_cx_inv); 
  outEz[2] += 1.7320508075688772*bx[0]*c2*dx1*jacob_cy_inv; 
  outEz[3] += 1.7320508075688772*bx[1]*c2*dx1*jacob_cy_inv-1.7320508075688772*by[2]*c2*dx0*jacob_cx_inv; 

  outBx[1] += 1.7320508075688772*ps[0]*dx0*jacob_cx_inv*gamma; 
  outBx[2] += 1.7320508075688772*ez[0]*dx1*jacob_cy_inv; 
  outBx[3] += 1.7320508075688772*ps[2]*dx0*jacob_cx_inv*gamma+1.7320508075688772*ez[1]*dx1*jacob_cy_inv; 

  outBy[1] += -(1.7320508075688772*ez[0]*dx0*jacob_cx_inv); 
  outBy[2] += 1.7320508075688772*ps[0]*dx1*jacob_cy_inv*gamma; 
  outBy[3] += 1.7320508075688772*ps[1]*dx1*jacob_cy_inv*gamma-1.7320508075688772*ez[2]*dx0*jacob_cx_inv; 

  outBz[1] += 1.7320508075688772*ey[0]*dx0*jacob_cx_inv; 
  outBz[2] += -(1.7320508075688772*ex[0]*dx1*jacob_cy_inv); 
  outBz[3] += 1.7320508075688772*ey[2]*dx0*jacob_cx_inv-1.7320508075688772*ex[1]*dx1*jacob_cy_inv; 

  outPh[1] += 1.7320508075688772*ex[0]*chi*dx0*jacob_cx_inv; 
  outPh[2] += 1.7320508075688772*ey[0]*chi*dx1*jacob_cy_inv; 
  outPh[3] += 1.7320508075688772*ey[1]*chi*dx1*jacob_cy_inv+1.7320508075688772*ex[2]*chi*dx0*jacob_cx_inv; 

  outPs[1] += 1.7320508075688772*bx[0]*c2gamma*dx0*jacob_cx_inv; 
  outPs[2] += 1.7320508075688772*by[0]*c2gamma*dx1*jacob_cy_inv; 
  outPs[3] += 1.7320508075688772*by[1]*c2gamma*dx1*jacob_cy_inv+1.7320508075688772*bx[2]*c2gamma*dx0*jacob_cx_inv; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]*jacob_cx_inv; 
  cflFreq += meq->c/dx[1]*jacob_cy_inv; 
  return 3.0*cflFreq; 
} 
