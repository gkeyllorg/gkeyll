#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_vol_1x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, double gas_gamma, 
    const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing (computational).
  // jacob_pos: Position-map Jacobian per configuration direction (cdim blocks of 2 equal
  //            entries): the physical size of the cell in direction d is J_d dxv[d].
  // gas_gamma: Adiabatic index.
  // u:         Input flow velocity [ux, uy, uz].
  // p:         Input pressure .
  // fluid:     [rho, rho ux, rho uy, rho uz, E], Fluid input state vector.
  // out:       Incremented output.

  double dx10 = 2.0/(dxv[0]*jacob_pos[0]); 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[2]; 
  const double *rhouy = &fluid[4]; 
  const double *rhouz = &fluid[6]; 
  const double *energy = &fluid[8]; 

  const double *ux = &u[0]; 
  const double *uy = &u[2]; 
  const double *uz = &u[4]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[2]; 
  double *outrhouy = &out[4]; 
  double *outrhouz = &out[6]; 
  double *outenergy = &out[8]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*3.0*dx10*(fabs(0.7071067811865475*ux[0]) + sqrt(gas_gamma*0.7071067811865475*p[0]/(0.7071067811865475*rho[0]))); 

  outrho[1] += 1.7320508075688772*rhoux[0]*dx10; 

  outrhoux[1] += 1.224744871391589*rhoux[1]*ux[1]*dx10+1.224744871391589*rhoux[0]*ux[0]*dx10+1.7320508075688772*p[0]*dx10; 

  outrhouy[1] += 1.224744871391589*rhouy[1]*ux[1]*dx10+1.224744871391589*rhouy[0]*ux[0]*dx10; 

  outrhouz[1] += 1.224744871391589*rhouz[1]*ux[1]*dx10+1.224744871391589*rhouz[0]*ux[0]*dx10; 

  outenergy[1] += 1.224744871391589*p[1]*ux[1]*dx10+1.224744871391589*energy[1]*ux[1]*dx10+1.224744871391589*p[0]*ux[0]*dx10+1.224744871391589*energy[0]*ux[0]*dx10; 

  return cflFreq_mid; 
} 
