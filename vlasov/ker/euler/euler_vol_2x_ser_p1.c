#include <gkyl_euler_kernels.h> 
GKYL_CU_DH double euler_vol_2x_ser_p1(const double *w, const double *dxv, const double *jacob_pos, double gas_gamma, 
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
  double dx11 = 2.0/(dxv[1]*jacob_pos[2]); 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[4]; 
  const double *rhouy = &fluid[8]; 
  const double *rhouz = &fluid[12]; 
  const double *energy = &fluid[16]; 

  const double *ux = &u[0]; 
  const double *uy = &u[4]; 
  const double *uz = &u[8]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[4]; 
  double *outrhouy = &out[8]; 
  double *outrhouz = &out[12]; 
  double *outenergy = &out[16]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*3.0*dx10*(fabs(0.5*ux[0]) + sqrt(gas_gamma*0.5*p[0]/(0.5*rho[0]))); 
  cflFreq_mid += 0.5*3.0*dx11*(fabs(0.5*uy[0]) + sqrt(gas_gamma*0.5*p[0]/(0.5*rho[0]))); 

  outrho[1] += 1.7320508075688772*rhoux[0]*dx10; 
  outrho[2] += 1.7320508075688772*rhouy[0]*dx11; 
  outrho[3] += 1.7320508075688772*rhouy[1]*dx11+1.7320508075688772*rhoux[2]*dx10; 

  outrhoux[1] += 0.8660254037844386*rhoux[3]*ux[3]*dx10+0.8660254037844386*rhoux[2]*ux[2]*dx10+0.8660254037844386*rhoux[1]*ux[1]*dx10+0.8660254037844386*rhoux[0]*ux[0]*dx10+1.7320508075688772*p[0]*dx10; 
  outrhoux[2] += 0.8660254037844386*rhoux[3]*uy[3]*dx11+0.8660254037844386*rhoux[2]*uy[2]*dx11+0.8660254037844386*rhoux[1]*uy[1]*dx11+0.8660254037844386*rhoux[0]*uy[0]*dx11; 
  outrhoux[3] += 0.8660254037844386*rhoux[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhoux[3]*dx11+0.8660254037844386*rhoux[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhoux[1]*dx11+0.8660254037844386*rhoux[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhoux[3]*dx10+0.8660254037844386*rhoux[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhoux[2]*dx10+1.7320508075688772*p[2]*dx10; 

  outrhouy[1] += 0.8660254037844386*rhouy[3]*ux[3]*dx10+0.8660254037844386*rhouy[2]*ux[2]*dx10+0.8660254037844386*rhouy[1]*ux[1]*dx10+0.8660254037844386*rhouy[0]*ux[0]*dx10; 
  outrhouy[2] += 0.8660254037844386*rhouy[3]*uy[3]*dx11+0.8660254037844386*rhouy[2]*uy[2]*dx11+0.8660254037844386*rhouy[1]*uy[1]*dx11+0.8660254037844386*rhouy[0]*uy[0]*dx11+1.7320508075688772*p[0]*dx11; 
  outrhouy[3] += 0.8660254037844386*rhouy[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhouy[3]*dx11+0.8660254037844386*rhouy[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhouy[1]*dx11+1.7320508075688772*p[1]*dx11+0.8660254037844386*rhouy[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhouy[3]*dx10+0.8660254037844386*rhouy[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhouy[2]*dx10; 

  outrhouz[1] += 0.8660254037844386*rhouz[3]*ux[3]*dx10+0.8660254037844386*rhouz[2]*ux[2]*dx10+0.8660254037844386*rhouz[1]*ux[1]*dx10+0.8660254037844386*rhouz[0]*ux[0]*dx10; 
  outrhouz[2] += 0.8660254037844386*rhouz[3]*uy[3]*dx11+0.8660254037844386*rhouz[2]*uy[2]*dx11+0.8660254037844386*rhouz[1]*uy[1]*dx11+0.8660254037844386*rhouz[0]*uy[0]*dx11; 
  outrhouz[3] += 0.8660254037844386*rhouz[2]*uy[3]*dx11+0.8660254037844386*uy[2]*rhouz[3]*dx11+0.8660254037844386*rhouz[0]*uy[1]*dx11+0.8660254037844386*uy[0]*rhouz[1]*dx11+0.8660254037844386*rhouz[1]*ux[3]*dx10+0.8660254037844386*ux[1]*rhouz[3]*dx10+0.8660254037844386*rhouz[0]*ux[2]*dx10+0.8660254037844386*ux[0]*rhouz[2]*dx10; 

  outenergy[1] += 0.8660254037844386*p[3]*ux[3]*dx10+0.8660254037844386*energy[3]*ux[3]*dx10+0.8660254037844386*p[2]*ux[2]*dx10+0.8660254037844386*energy[2]*ux[2]*dx10+0.8660254037844386*p[1]*ux[1]*dx10+0.8660254037844386*energy[1]*ux[1]*dx10+0.8660254037844386*p[0]*ux[0]*dx10+0.8660254037844386*energy[0]*ux[0]*dx10; 
  outenergy[2] += 0.8660254037844386*p[3]*uy[3]*dx11+0.8660254037844386*energy[3]*uy[3]*dx11+0.8660254037844386*p[2]*uy[2]*dx11+0.8660254037844386*energy[2]*uy[2]*dx11+0.8660254037844386*p[1]*uy[1]*dx11+0.8660254037844386*energy[1]*uy[1]*dx11+0.8660254037844386*p[0]*uy[0]*dx11+0.8660254037844386*energy[0]*uy[0]*dx11; 
  outenergy[3] += 0.8660254037844386*p[2]*uy[3]*dx11+0.8660254037844386*energy[2]*uy[3]*dx11+0.8660254037844386*uy[2]*p[3]*dx11+0.8660254037844386*uy[2]*energy[3]*dx11+0.8660254037844386*p[0]*uy[1]*dx11+0.8660254037844386*energy[0]*uy[1]*dx11+0.8660254037844386*uy[0]*p[1]*dx11+0.8660254037844386*uy[0]*energy[1]*dx11+0.8660254037844386*p[1]*ux[3]*dx10+0.8660254037844386*energy[1]*ux[3]*dx10+0.8660254037844386*ux[1]*p[3]*dx10+0.8660254037844386*ux[1]*energy[3]*dx10+0.8660254037844386*p[0]*ux[2]*dx10+0.8660254037844386*energy[0]*ux[2]*dx10+0.8660254037844386*ux[0]*p[2]*dx10+0.8660254037844386*ux[0]*energy[2]*dx10; 

  return cflFreq_mid; 
} 
