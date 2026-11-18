#include <gkyl_canonical_fluid_kernels.h> 
void canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p1(const double *dxv, double alpha, const double *phi, const double *n0, const double *adiabatic_coupling_phi_n, double* GKYL_RESTRICT out) 
{ 
  // dxv[NDIM]: cell length.
  // alpha: Adiabaticity parameter for adiabatic coupling of vorticity and density (zero for Hasegawa-Mima).
  // phi: Potential in fluid system given by grad^2 phi = f where f is (one of) the evolved quantities.
  // n0: Background density gradient.
  // adiabatic_coupling_phi_n: (phi, n) array for adiabatic coupling, with potentially zonal component subtracted out.
  // out: output increment in center cell.
  double dx = 2.0/dxv[0]; 
  double dy = 2.0/dxv[1]; 
  const double *phi_adiabat = &adiabatic_coupling_phi_n[0]; 
  const double *n_adiabat = &adiabatic_coupling_phi_n[4]; 
  double *out_zeta = &out[0]; 
  double *out_n = &out[4]; 
  out_zeta[0] += phi_adiabat[0]*alpha-1.0*n_adiabat[0]*alpha; 
  out_zeta[1] += phi_adiabat[1]*alpha-1.0*n_adiabat[1]*alpha; 
  out_zeta[2] += phi_adiabat[2]*alpha-1.0*n_adiabat[2]*alpha; 
  out_zeta[3] += phi_adiabat[3]*alpha-1.0*n_adiabat[3]*alpha; 
  out_n[0] += -(1.5*n0[1]*phi[2]*dx*dy)+1.5*phi[1]*n0[2]*dx*dy+phi_adiabat[0]*alpha-1.0*n_adiabat[0]*alpha; 
  out_n[1] += -(1.5*n0[1]*phi[3]*dx*dy)+1.5*phi[1]*n0[3]*dx*dy+phi_adiabat[1]*alpha-1.0*n_adiabat[1]*alpha; 
  out_n[2] += 1.5*n0[2]*phi[3]*dx*dy-1.5*phi[2]*n0[3]*dx*dy+phi_adiabat[2]*alpha-1.0*n_adiabat[2]*alpha; 
  out_n[3] += phi_adiabat[3]*alpha-1.0*n_adiabat[3]*alpha; 
} 
