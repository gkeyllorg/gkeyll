#include <gkyl_iso_euler_kernels.h> 
#include <math.h> 
GKYL_CU_DH double fluid_vars_iso_em_source_1x_ser_p1(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out) 
{ 
  // qmem:       [ax, ay, az, q/m Bx, q/m By, q/m Bz], acceleration (q/m E plus any applied
  //             acceleration) and q/m times the magnetic field.
  // omega_p_sq: Plasma frequency squared of the cell (one value; all the coupled species).
  // fluid:      [rho, rho ux, rho uy, rho uz], Fluid input state vector.
  // out:        Output increment
  // Returns the rate of the source for the time-step (from the cell averages).

  const double *ax = &qmem[0]; 
  const double *ay = &qmem[2]; 
  const double *az = &qmem[4]; 
  const double *Bx = &qmem[6]; 
  const double *By = &qmem[8]; 
  const double *Bz = &qmem[10]; 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[2]; 
  const double *rhouy = &fluid[4]; 
  const double *rhouz = &fluid[6]; 

  double *out_rhoux = &out[2]; 
  double *out_rhouy = &out[4]; 
  double *out_rhouz = &out[6]; 

  out_rhoux[0] += -(0.7071067811865475*By[1]*rhouz[1])+0.7071067811865475*Bz[1]*rhouy[1]+0.7071067811865475*ax[1]*rho[1]-0.7071067811865475*By[0]*rhouz[0]+0.7071067811865475*Bz[0]*rhouy[0]+0.7071067811865475*ax[0]*rho[0]; 
  out_rhoux[1] += -(0.7071067811865475*By[0]*rhouz[1])+0.7071067811865475*Bz[0]*rhouy[1]+0.7071067811865475*ax[0]*rho[1]+0.7071067811865475*rho[0]*ax[1]+0.7071067811865475*rhouy[0]*Bz[1]-0.7071067811865475*rhouz[0]*By[1]; 

  out_rhouy[0] += 0.7071067811865475*Bx[1]*rhouz[1]-0.7071067811865475*Bz[1]*rhoux[1]+0.7071067811865475*ay[1]*rho[1]+0.7071067811865475*Bx[0]*rhouz[0]-0.7071067811865475*Bz[0]*rhoux[0]+0.7071067811865475*ay[0]*rho[0]; 
  out_rhouy[1] += 0.7071067811865475*Bx[0]*rhouz[1]-0.7071067811865475*Bz[0]*rhoux[1]+0.7071067811865475*ay[0]*rho[1]+0.7071067811865475*rho[0]*ay[1]-0.7071067811865475*rhoux[0]*Bz[1]+0.7071067811865475*rhouz[0]*Bx[1]; 

  out_rhouz[0] += -(0.7071067811865475*Bx[1]*rhouy[1])+0.7071067811865475*By[1]*rhoux[1]+0.7071067811865475*az[1]*rho[1]-0.7071067811865475*Bx[0]*rhouy[0]+0.7071067811865475*By[0]*rhoux[0]+0.7071067811865475*az[0]*rho[0]; 
  out_rhouz[1] += -(0.7071067811865475*Bx[0]*rhouy[1])+0.7071067811865475*By[0]*rhoux[1]+0.7071067811865475*az[0]*rho[1]+0.7071067811865475*rho[0]*az[1]+0.7071067811865475*rhoux[0]*By[1]-0.7071067811865475*rhouy[0]*Bx[1]; 

  const double avg = 0.7071067811865475; 
  const double Bx_a = avg*Bx[0], By_a = avg*By[0], Bz_a = avg*Bz[0]; 
  return sqrt(omega_p_sq[0]) + sqrt(Bx_a*Bx_a + By_a*By_a + Bz_a*Bz_a); 

} 
