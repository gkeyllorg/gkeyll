#include <gkyl_advection_kernels.h> 
GKYL_CU_DH double advection_vol_1x_ser_p3(const double *w, const double *dxv, const double *u, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // u[NDIM]:   Advection velocity.
  // f:         Input function.
  // out:       Incremented output.
  const double rdx2 = 2.0/dxv[0]; 
  double cflFreq_mid = 0.0; 
  cflFreq_mid += fabs((2.474873734152916*u[0]-2.7669929526473314*u[2])*rdx2); 

  out[1] += 1.224744871391589*(f[3]*u[3]+f[2]*u[2]+f[1]*u[1]+f[0]*u[0])*rdx2; 
  out[2] += (2.405351177211819*(f[2]*u[3]+u[2]*f[3])+2.4494897427831783*(f[1]*u[2]+u[1]*f[2])+2.7386127875258306*(f[0]*u[1]+u[0]*f[1]))*rdx2; 
  out[3] += (4.365266951236265*f[3]*u[3]+3.6742346141747664*(f[1]*u[3]+u[1]*f[3])+4.543441112511214*f[2]*u[2]+4.183300132670378*(f[0]*u[2]+u[0]*f[2])+5.612486080160912*f[1]*u[1]+1.8708286933869707*f[0]*u[0])*rdx2; 

  return cflFreq_mid; 
} 
