#include <gkyl_iso_euler_kernels.h> 
GKYL_CU_DH double iso_euler_vol_1x_ser_p1(const double *w, const double *dxv, double vt, 
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
  const double *rhoux = &fluid[2]; 
  const double *rhouy = &fluid[4]; 
  const double *rhouz = &fluid[6]; 

  const double *ux = &u[0]; 
  const double *uy = &u[2]; 
  const double *uz = &u[4]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[2]; 
  double *outrhouy = &out[4]; 
  double *outrhouz = &out[6]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*3.0*dx10*(fabs(0.7071067811865475*ux[0]) + vt); 

  outrho[1] += 1.7320508075688772*rhoux[0]*dx10; 

  outrhoux[1] += 1.224744871391589*rhoux[1]*ux[1]*dx10+1.224744871391589*rhoux[0]*ux[0]*dx10+1.7320508075688772*p[0]*dx10; 

  outrhouy[1] += 1.224744871391589*rhouy[1]*ux[1]*dx10+1.224744871391589*rhouy[0]*ux[0]*dx10; 

  outrhouz[1] += 1.224744871391589*rhouz[1]*ux[1]*dx10+1.224744871391589*rhouz[0]*ux[0]*dx10; 

  return cflFreq_mid; 
} 
