#include <math.h>
#include <gkyl_moment_prim_sr_euler.h>

#define TU 2
#define TV 3
#define TW 4

double
gkyl_sr_euler_max_abs_speed(double gas_gamma, const double q[5])
{
  double v[5];
  gkyl_sr_euler_prim_vars(gas_gamma, q, v);

  double pr = v[1];
  double gamma = 1. / sqrt(1. - (v[2]*v[2] + v[3]*v[3] + v[4]*v[4]));
  double fac0 = q[1] + pr;
  double v4 = gamma*gamma*pr/fac0;
  double fac1 = 1 - gas_gamma*v4;
  double fac2 = gas_gamma*pr*fac1*(fac0 - q[2]*q[2] / fac0) + gas_gamma*gas_gamma*pr*pr;
  
  return (fac1*q[2] + sqrt(fac2)) / (fac1*fac0 + gas_gamma*pr);
}

void
gkyl_sr_euler_flux(double gas_gamma, const double q[5], double flux[5])
{
  double v[5];
  gkyl_sr_euler_prim_vars(gas_gamma, q, v);
  double pr = v[1];

  double fac0 = q[1] + pr;
  flux[0] = q[0]*q[TU]/fac0; // gamma*rho*u
  flux[1] = q[TU]; //gamma^2*rho*h*u
  flux[TU] = q[TU]*q[TU]/fac0 + pr; // gamma^2*rho*h*u*u + pr
  flux[TV] = q[TU]*q[TV]/fac0; // gamma^2*rho*h*u*v
  flux[TW] = q[TU]*q[TW]/fac0; // gamma^2*rho*h*u*w
}
