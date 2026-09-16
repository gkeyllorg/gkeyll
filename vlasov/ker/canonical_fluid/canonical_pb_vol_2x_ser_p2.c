#include <gkyl_canonical_fluid_kernels.h> 
double canonical_pb_vol_2x_ser_p2(const double *w, const double *dxv, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // phi: Potential in fluid system given by grad^2 phi = f where f is (one of) the evolved quantities.
  // f: input state vector in center cell.
  // out: output increment in center cell.
  double dxdvInv = 4.0/(dxv[0]*dxv[1]); 
  out[1] += (3.3541019662496843*f[3]*phi[7]+1.5*f[4]*phi[6]+3.3541019662496847*f[2]*phi[5]+1.5*f[1]*phi[3]+1.5*f[0]*phi[2])*dxdvInv; 
  out[2] += (-(1.5*f[5]*phi[7])-3.3541019662496843*f[3]*phi[6]-3.3541019662496847*f[1]*phi[4]-1.5*f[2]*phi[3]-1.5*f[0]*phi[1])*dxdvInv; 
  out[3] += (1.5*f[7]*phi[7]+3.3541019662496843*f[1]*phi[7]-1.5*f[6]*phi[6]-3.3541019662496843*f[2]*phi[6]+3.0*f[5]*phi[5]+3.3541019662496847*f[0]*phi[5]-3.0*f[4]*phi[4]-3.3541019662496847*f[0]*phi[4]+1.5*f[2]*phi[2]-1.5*f[1]*phi[1])*dxdvInv; 
  out[4] += (6.708203932499369*f[6]*phi[7]+7.500000000000001*f[2]*phi[7]+3.0*f[1]*phi[6]+7.5*f[3]*phi[5]+3.0*phi[3]*f[4]+3.3541019662496847*f[0]*phi[3]+3.3541019662496847*f[1]*phi[2])*dxdvInv; 
  out[5] += (-(3.0*f[2]*phi[7])-6.708203932499369*phi[6]*f[7]-7.500000000000001*f[1]*phi[6]-3.0*phi[3]*f[5]-7.5*f[3]*phi[4]-3.3541019662496847*f[0]*phi[3]-3.3541019662496847*phi[1]*f[2])*dxdvInv; 
  out[6] += (6.708203932499369*f[5]*phi[7]+6.708203932499369*f[4]*phi[7]+7.5*f[0]*phi[7]+6.708203932499369*phi[5]*f[7]+1.5*phi[3]*f[6]+7.500000000000001*f[1]*phi[5]-3.0*f[1]*phi[4]-1.5*phi[1]*f[4]+3.3541019662496843*f[2]*phi[3]+3.3541019662496843*phi[2]*f[3])*dxdvInv; 
  out[7] += (-(1.5*phi[3]*f[7])-6.708203932499369*f[5]*phi[6]-6.708203932499369*f[4]*phi[6]-7.5*f[0]*phi[6]-6.708203932499369*phi[4]*f[6]+3.0*f[2]*phi[5]+1.5*phi[2]*f[5]-7.500000000000001*f[2]*phi[4]-3.3541019662496843*f[1]*phi[3]-3.3541019662496843*phi[1]*f[3])*dxdvInv; 
  return 0.; 
} 
