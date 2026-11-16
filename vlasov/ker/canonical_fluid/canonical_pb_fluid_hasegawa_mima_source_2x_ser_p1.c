#include <gkyl_canonical_fluid_kernels.h> 
void canonical_pb_fluid_hasegawa_mima_source_2x_ser_p1(const double *dxv, double alpha, const double *phi, const double *n0, const double *f, double* GKYL_RESTRICT out) 
{ 
  // dxv[NDIM]: cell length.
  // alpha: Adiabaticity parameter for adiabatic coupling of vorticity and density (zero for Hasegawa-Mima).
  // phi: Potential in fluid system given by grad^2 phi = f where f is (one of) the evolved quantities.
  // n0: Background density gradient.
  // f: input state vector in center cell.
  // out: output increment in center cell.
  double dx = 2.0/dxv[0]; 
  double dy = 2.0/dxv[1]; 
  out[0] += 1.5*phi[1]*n0[2]*dx*dy-1.5*n0[1]*phi[2]*dx*dy; 
  out[1] += 1.5*phi[1]*n0[3]*dx*dy-1.5*n0[1]*phi[3]*dx*dy; 
  out[2] += 1.5*n0[2]*phi[3]*dx*dy-1.5*phi[2]*n0[3]*dx*dy; 
} 
