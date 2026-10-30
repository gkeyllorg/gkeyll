#include <gkyl_canonical_fluid_kernels.h> 
double canonical_pb_vol_2x_ser_p1(const double *w, const double *dxv, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // phi: Potential in fluid system given by grad^2 phi = f where f is (one of) the evolved quantities.
  // f: input state vector in center cell.
  // out: output increment in center cell.
  double dxdvInv = 4.0/(dxv[0]*dxv[1]); 
  out[1] += (1.5*f[1]*phi[3]+1.5*f[0]*phi[2])*dxdvInv; 
  out[2] += (-(1.5*f[2]*phi[3])-1.5*f[0]*phi[1])*dxdvInv; 
  out[3] += (1.5*f[2]*phi[2]-1.5*f[1]*phi[1])*dxdvInv; 
  return 0.; 
} 
