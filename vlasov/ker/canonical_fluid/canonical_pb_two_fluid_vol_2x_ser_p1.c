#include <gkyl_canonical_fluid_kernels.h> 
double canonical_pb_two_fluid_vol_2x_ser_p1(const double *w, const double *dxv, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // phi: Potential in fluid system given by grad^2 phi = f where f is (one of) the evolved quantities.
  // f: input state vector in center cell.
  // out: output increment in center cell.
  double dxdvInv = 4.0/(dxv[0]*dxv[1]); 
  const double *f1 = &f[0]; 
  const double *f2 = &f[4]; 
  double *out1 = &out[0]; 
  double *out2 = &out[4]; 
  out1[1] += (1.5*f1[1]*phi[3]+1.5*f1[0]*phi[2])*dxdvInv; 
  out2[1] += (1.5*f2[1]*phi[3]+1.5*f2[0]*phi[2])*dxdvInv; 
  out1[2] += (-(1.5*f1[2]*phi[3])-1.5*f1[0]*phi[1])*dxdvInv; 
  out2[2] += (-(1.5*f2[2]*phi[3])-1.5*f2[0]*phi[1])*dxdvInv; 
  out1[3] += (1.5*f1[2]*phi[2]-1.5*f1[1]*phi[1])*dxdvInv; 
  out2[3] += (1.5*f2[2]*phi[2]-1.5*f2[1]*phi[1])*dxdvInv; 
  return 0.; 
} 
