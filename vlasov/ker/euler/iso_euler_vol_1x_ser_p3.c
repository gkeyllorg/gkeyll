#include <gkyl_iso_euler_kernels.h> 
GKYL_CU_DH double iso_euler_vol_1x_ser_p3(const double *w, const double *dxv, double vt, 
    const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // vt:        Thermal velocity.
  // u:         Input flow velocity [ux, uy, uz].
  // p:         Input pressure (vt^2 rho).
  // fluid:     [rho, rho ux, rho uy, rho uz], Fluid input state vector.
  // out:       Incremented output.

  double dx10 = 2./dxv[0]; 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[4]; 
  const double *rhouy = &fluid[8]; 
  const double *rhouz = &fluid[12]; 

  const double *ux = &u[0]; 
  const double *uy = &u[4]; 
  const double *uz = &u[8]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[4]; 
  double *outrhouy = &out[8]; 
  double *outrhouz = &out[12]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*7.0*dx10*(fabs(0.7071067811865475*ux[0]-0.7905694150420947*ux[2]) + vt); 

  outrho[1] += 1.7320508075688772*rhoux[0]*dx10; 
  outrho[2] += 3.872983346207417*rhoux[1]*dx10; 
  outrho[3] += 5.916079783099617*rhoux[2]*dx10+2.6457513110645907*rhoux[0]*dx10; 

  outrhoux[1] += 1.224744871391589*rhoux[3]*ux[3]*dx10+1.224744871391589*rhoux[2]*ux[2]*dx10+1.224744871391589*rhoux[1]*ux[1]*dx10+1.224744871391589*rhoux[0]*ux[0]*dx10+1.7320508075688772*p[0]*dx10; 
  outrhoux[2] += 2.405351177211819*rhoux[2]*ux[3]*dx10+2.405351177211819*ux[2]*rhoux[3]*dx10+2.4494897427831783*rhoux[1]*ux[2]*dx10+2.4494897427831783*ux[1]*rhoux[2]*dx10+2.7386127875258306*rhoux[0]*ux[1]*dx10+2.7386127875258306*ux[0]*rhoux[1]*dx10+3.872983346207417*p[1]*dx10; 
  outrhoux[3] += 4.365266951236265*rhoux[3]*ux[3]*dx10+3.6742346141747664*rhoux[1]*ux[3]*dx10+3.6742346141747664*ux[1]*rhoux[3]*dx10+4.543441112511214*rhoux[2]*ux[2]*dx10+4.183300132670378*rhoux[0]*ux[2]*dx10+4.183300132670378*ux[0]*rhoux[2]*dx10+5.916079783099617*p[2]*dx10+5.612486080160912*rhoux[1]*ux[1]*dx10+1.8708286933869707*rhoux[0]*ux[0]*dx10+2.6457513110645907*p[0]*dx10; 

  outrhouy[1] += 1.224744871391589*rhouy[3]*ux[3]*dx10+1.224744871391589*rhouy[2]*ux[2]*dx10+1.224744871391589*rhouy[1]*ux[1]*dx10+1.224744871391589*rhouy[0]*ux[0]*dx10; 
  outrhouy[2] += 2.405351177211819*rhouy[2]*ux[3]*dx10+2.405351177211819*ux[2]*rhouy[3]*dx10+2.4494897427831783*rhouy[1]*ux[2]*dx10+2.4494897427831783*ux[1]*rhouy[2]*dx10+2.7386127875258306*rhouy[0]*ux[1]*dx10+2.7386127875258306*ux[0]*rhouy[1]*dx10; 
  outrhouy[3] += 4.365266951236265*rhouy[3]*ux[3]*dx10+3.6742346141747664*rhouy[1]*ux[3]*dx10+3.6742346141747664*ux[1]*rhouy[3]*dx10+4.543441112511214*rhouy[2]*ux[2]*dx10+4.183300132670378*rhouy[0]*ux[2]*dx10+4.183300132670378*ux[0]*rhouy[2]*dx10+5.612486080160912*rhouy[1]*ux[1]*dx10+1.8708286933869707*rhouy[0]*ux[0]*dx10; 

  outrhouz[1] += 1.224744871391589*rhouz[3]*ux[3]*dx10+1.224744871391589*rhouz[2]*ux[2]*dx10+1.224744871391589*rhouz[1]*ux[1]*dx10+1.224744871391589*rhouz[0]*ux[0]*dx10; 
  outrhouz[2] += 2.405351177211819*rhouz[2]*ux[3]*dx10+2.405351177211819*ux[2]*rhouz[3]*dx10+2.4494897427831783*rhouz[1]*ux[2]*dx10+2.4494897427831783*ux[1]*rhouz[2]*dx10+2.7386127875258306*rhouz[0]*ux[1]*dx10+2.7386127875258306*ux[0]*rhouz[1]*dx10; 
  outrhouz[3] += 4.365266951236265*rhouz[3]*ux[3]*dx10+3.6742346141747664*rhouz[1]*ux[3]*dx10+3.6742346141747664*ux[1]*rhouz[3]*dx10+4.543441112511214*rhouz[2]*ux[2]*dx10+4.183300132670378*rhouz[0]*ux[2]*dx10+4.183300132670378*ux[0]*rhouz[2]*dx10+5.612486080160912*rhouz[1]*ux[1]*dx10+1.8708286933869707*rhouz[0]*ux[0]*dx10; 

  return cflFreq_mid; 
} 
