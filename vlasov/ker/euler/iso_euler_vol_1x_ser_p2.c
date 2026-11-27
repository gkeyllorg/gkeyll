#include <gkyl_iso_euler_kernels.h> 
GKYL_CU_DH double iso_euler_vol_1x_ser_p2(const double *w, const double *dxv, const double *jacob_pos, double vt, 
    const double *u, const double *p, const double *fluid, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing (computational).
  // jacob_pos: Position-map Jacobian per configuration direction (cdim blocks of 3 equal
  //            entries): the physical size of the cell in direction d is J_d dxv[d].
  // vt:        Thermal velocity.
  // u:         Input flow velocity [ux, uy, uz].
  // p:         Input pressure (vt^2 rho).
  // fluid:     [rho, rho ux, rho uy, rho uz], Fluid input state vector.
  // out:       Incremented output.

  double dx10 = 2.0/(dxv[0]*jacob_pos[0]); 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[3]; 
  const double *rhouy = &fluid[6]; 
  const double *rhouz = &fluid[9]; 

  const double *ux = &u[0]; 
  const double *uy = &u[3]; 
  const double *uz = &u[6]; 

  double *outrho = &out[0]; 
  double *outrhoux = &out[3]; 
  double *outrhouy = &out[6]; 
  double *outrhouz = &out[9]; 

  double cflFreq_mid = 0.0; 
  cflFreq_mid += 0.5*5.0*dx10*(fabs(0.7071067811865475*ux[0]-0.7905694150420947*ux[2]) + vt); 

  outrho[1] += 1.7320508075688772*rhoux[0]*dx10; 
  outrho[2] += 3.872983346207417*rhoux[1]*dx10; 

  outrhoux[1] += 1.224744871391589*rhoux[2]*ux[2]*dx10+1.224744871391589*rhoux[1]*ux[1]*dx10+1.224744871391589*rhoux[0]*ux[0]*dx10+1.7320508075688772*p[0]*dx10; 
  outrhoux[2] += 2.4494897427831783*rhoux[1]*ux[2]*dx10+2.4494897427831783*ux[1]*rhoux[2]*dx10+2.7386127875258306*rhoux[0]*ux[1]*dx10+2.7386127875258306*ux[0]*rhoux[1]*dx10+3.872983346207417*p[1]*dx10; 

  outrhouy[1] += 1.224744871391589*rhouy[2]*ux[2]*dx10+1.224744871391589*rhouy[1]*ux[1]*dx10+1.224744871391589*rhouy[0]*ux[0]*dx10; 
  outrhouy[2] += 2.4494897427831783*rhouy[1]*ux[2]*dx10+2.4494897427831783*ux[1]*rhouy[2]*dx10+2.7386127875258306*rhouy[0]*ux[1]*dx10+2.7386127875258306*ux[0]*rhouy[1]*dx10; 

  outrhouz[1] += 1.224744871391589*rhouz[2]*ux[2]*dx10+1.224744871391589*rhouz[1]*ux[1]*dx10+1.224744871391589*rhouz[0]*ux[0]*dx10; 
  outrhouz[2] += 2.4494897427831783*rhouz[1]*ux[2]*dx10+2.4494897427831783*ux[1]*rhouz[2]*dx10+2.7386127875258306*rhouz[0]*ux[1]*dx10+2.7386127875258306*ux[0]*rhouz[1]*dx10; 

  return cflFreq_mid; 
} 
