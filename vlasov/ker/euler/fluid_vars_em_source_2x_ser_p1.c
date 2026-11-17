#include <gkyl_euler_kernels.h> 
#include <math.h> 
GKYL_CU_DH double fluid_vars_em_source_2x_ser_p1(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out) 
{ 
  // qmem:       [ax, ay, az, q/m Bx, q/m By, q/m Bz], acceleration (q/m E plus any applied
  //             acceleration) and q/m times the magnetic field.
  // omega_p_sq: Plasma frequency squared of the cell (one value; all the coupled species).
  // fluid:      [rho, rho ux, rho uy, rho uz, energy], Fluid input state vector.
  // out:        Output increment
  // Returns the rate of the source for the time-step (from the cell averages).

  const double *ax = &qmem[0]; 
  const double *ay = &qmem[4]; 
  const double *az = &qmem[8]; 
  const double *Bx = &qmem[12]; 
  const double *By = &qmem[16]; 
  const double *Bz = &qmem[20]; 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[4]; 
  const double *rhouy = &fluid[8]; 
  const double *rhouz = &fluid[12]; 

  double *out_rhoux = &out[4]; 
  double *out_rhouy = &out[8]; 
  double *out_rhouz = &out[12]; 
  double *out_energy = &out[16]; 

  out_rhoux[0] += -(0.5*By[3]*rhouz[3])+0.5*Bz[3]*rhouy[3]+0.5*ax[3]*rho[3]-0.5*By[2]*rhouz[2]+0.5*Bz[2]*rhouy[2]+0.5*ax[2]*rho[2]-0.5*By[1]*rhouz[1]+0.5*Bz[1]*rhouy[1]+0.5*ax[1]*rho[1]-0.5*By[0]*rhouz[0]+0.5*Bz[0]*rhouy[0]+0.5*ax[0]*rho[0]; 
  out_rhoux[1] += -(0.5*By[2]*rhouz[3])+0.5*Bz[2]*rhouy[3]+0.5*ax[2]*rho[3]+0.5*rho[2]*ax[3]+0.5*rhouy[2]*Bz[3]-0.5*rhouz[2]*By[3]-0.5*By[0]*rhouz[1]+0.5*Bz[0]*rhouy[1]+0.5*ax[0]*rho[1]+0.5*rho[0]*ax[1]+0.5*rhouy[0]*Bz[1]-0.5*rhouz[0]*By[1]; 
  out_rhoux[2] += -(0.5*By[1]*rhouz[3])+0.5*Bz[1]*rhouy[3]+0.5*ax[1]*rho[3]+0.5*rho[1]*ax[3]+0.5*rhouy[1]*Bz[3]-0.5*rhouz[1]*By[3]-0.5*By[0]*rhouz[2]+0.5*Bz[0]*rhouy[2]+0.5*ax[0]*rho[2]+0.5*rho[0]*ax[2]+0.5*rhouy[0]*Bz[2]-0.5*rhouz[0]*By[2]; 
  out_rhoux[3] += -(0.5*By[0]*rhouz[3])+0.5*Bz[0]*rhouy[3]+0.5*ax[0]*rho[3]+0.5*rho[0]*ax[3]+0.5*rhouy[0]*Bz[3]-0.5*rhouz[0]*By[3]-0.5*By[1]*rhouz[2]+0.5*Bz[1]*rhouy[2]+0.5*ax[1]*rho[2]+0.5*rho[1]*ax[2]+0.5*rhouy[1]*Bz[2]-0.5*rhouz[1]*By[2]; 

  out_rhouy[0] += 0.5*Bx[3]*rhouz[3]-0.5*Bz[3]*rhoux[3]+0.5*ay[3]*rho[3]+0.5*Bx[2]*rhouz[2]-0.5*Bz[2]*rhoux[2]+0.5*ay[2]*rho[2]+0.5*Bx[1]*rhouz[1]-0.5*Bz[1]*rhoux[1]+0.5*ay[1]*rho[1]+0.5*Bx[0]*rhouz[0]-0.5*Bz[0]*rhoux[0]+0.5*ay[0]*rho[0]; 
  out_rhouy[1] += 0.5*Bx[2]*rhouz[3]-0.5*Bz[2]*rhoux[3]+0.5*ay[2]*rho[3]+0.5*rho[2]*ay[3]-0.5*rhoux[2]*Bz[3]+0.5*rhouz[2]*Bx[3]+0.5*Bx[0]*rhouz[1]-0.5*Bz[0]*rhoux[1]+0.5*ay[0]*rho[1]+0.5*rho[0]*ay[1]-0.5*rhoux[0]*Bz[1]+0.5*rhouz[0]*Bx[1]; 
  out_rhouy[2] += 0.5*Bx[1]*rhouz[3]-0.5*Bz[1]*rhoux[3]+0.5*ay[1]*rho[3]+0.5*rho[1]*ay[3]-0.5*rhoux[1]*Bz[3]+0.5*rhouz[1]*Bx[3]+0.5*Bx[0]*rhouz[2]-0.5*Bz[0]*rhoux[2]+0.5*ay[0]*rho[2]+0.5*rho[0]*ay[2]-0.5*rhoux[0]*Bz[2]+0.5*rhouz[0]*Bx[2]; 
  out_rhouy[3] += 0.5*Bx[0]*rhouz[3]-0.5*Bz[0]*rhoux[3]+0.5*ay[0]*rho[3]+0.5*rho[0]*ay[3]-0.5*rhoux[0]*Bz[3]+0.5*rhouz[0]*Bx[3]+0.5*Bx[1]*rhouz[2]-0.5*Bz[1]*rhoux[2]+0.5*ay[1]*rho[2]+0.5*rho[1]*ay[2]-0.5*rhoux[1]*Bz[2]+0.5*rhouz[1]*Bx[2]; 

  out_rhouz[0] += -(0.5*Bx[3]*rhouy[3])+0.5*By[3]*rhoux[3]+0.5*az[3]*rho[3]-0.5*Bx[2]*rhouy[2]+0.5*By[2]*rhoux[2]+0.5*az[2]*rho[2]-0.5*Bx[1]*rhouy[1]+0.5*By[1]*rhoux[1]+0.5*az[1]*rho[1]-0.5*Bx[0]*rhouy[0]+0.5*By[0]*rhoux[0]+0.5*az[0]*rho[0]; 
  out_rhouz[1] += -(0.5*Bx[2]*rhouy[3])+0.5*By[2]*rhoux[3]+0.5*az[2]*rho[3]+0.5*rho[2]*az[3]+0.5*rhoux[2]*By[3]-0.5*rhouy[2]*Bx[3]-0.5*Bx[0]*rhouy[1]+0.5*By[0]*rhoux[1]+0.5*az[0]*rho[1]+0.5*rho[0]*az[1]+0.5*rhoux[0]*By[1]-0.5*rhouy[0]*Bx[1]; 
  out_rhouz[2] += -(0.5*Bx[1]*rhouy[3])+0.5*By[1]*rhoux[3]+0.5*az[1]*rho[3]+0.5*rho[1]*az[3]+0.5*rhoux[1]*By[3]-0.5*rhouy[1]*Bx[3]-0.5*Bx[0]*rhouy[2]+0.5*By[0]*rhoux[2]+0.5*az[0]*rho[2]+0.5*rho[0]*az[2]+0.5*rhoux[0]*By[2]-0.5*rhouy[0]*Bx[2]; 
  out_rhouz[3] += -(0.5*Bx[0]*rhouy[3])+0.5*By[0]*rhoux[3]+0.5*az[0]*rho[3]+0.5*rho[0]*az[3]+0.5*rhoux[0]*By[3]-0.5*rhouy[0]*Bx[3]-0.5*Bx[1]*rhouy[2]+0.5*By[1]*rhoux[2]+0.5*az[1]*rho[2]+0.5*rho[1]*az[2]+0.5*rhoux[1]*By[2]-0.5*rhouy[1]*Bx[2]; 

  out_energy[0] += 0.5*az[3]*rhouz[3]+0.5*ay[3]*rhouy[3]+0.5*ax[3]*rhoux[3]+0.5*az[2]*rhouz[2]+0.5*ay[2]*rhouy[2]+0.5*ax[2]*rhoux[2]+0.5*az[1]*rhouz[1]+0.5*ay[1]*rhouy[1]+0.5*ax[1]*rhoux[1]+0.5*az[0]*rhouz[0]+0.5*ay[0]*rhouy[0]+0.5*ax[0]*rhoux[0]; 
  out_energy[1] += 0.5*az[2]*rhouz[3]+0.5*ay[2]*rhouy[3]+0.5*ax[2]*rhoux[3]+0.5*rhouz[2]*az[3]+0.5*rhouy[2]*ay[3]+0.5*rhoux[2]*ax[3]+0.5*az[0]*rhouz[1]+0.5*ay[0]*rhouy[1]+0.5*ax[0]*rhoux[1]+0.5*rhouz[0]*az[1]+0.5*rhouy[0]*ay[1]+0.5*rhoux[0]*ax[1]; 
  out_energy[2] += 0.5*az[1]*rhouz[3]+0.5*ay[1]*rhouy[3]+0.5*ax[1]*rhoux[3]+0.5*rhouz[1]*az[3]+0.5*rhouy[1]*ay[3]+0.5*rhoux[1]*ax[3]+0.5*az[0]*rhouz[2]+0.5*ay[0]*rhouy[2]+0.5*ax[0]*rhoux[2]+0.5*rhouz[0]*az[2]+0.5*rhouy[0]*ay[2]+0.5*rhoux[0]*ax[2]; 
  out_energy[3] += 0.5*az[0]*rhouz[3]+0.5*ay[0]*rhouy[3]+0.5*ax[0]*rhoux[3]+0.5*rhouz[0]*az[3]+0.5*rhouy[0]*ay[3]+0.5*rhoux[0]*ax[3]+0.5*az[1]*rhouz[2]+0.5*ay[1]*rhouy[2]+0.5*ax[1]*rhoux[2]+0.5*rhouz[1]*az[2]+0.5*rhouy[1]*ay[2]+0.5*rhoux[1]*ax[2]; 

  const double avg = 0.5; 
  const double rho_a = avg*rho[0], Mx_a = avg*rhoux[0], My_a = avg*rhouy[0], Mz_a = avg*rhouz[0]; 
  const double E_a = avg*fluid[16]; 
  const double ax_a = avg*ax[0], ay_a = avg*ay[0], az_a = avg*az[0]; 
  const double Bx_a = avg*Bx[0], By_a = avg*By[0], Bz_a = avg*Bz[0]; 
  const double fx = rho_a*ax_a + My_a*Bz_a - Mz_a*By_a; 
  const double fy = rho_a*ay_a + Mz_a*Bx_a - Mx_a*Bz_a; 
  const double fz = rho_a*az_a + Mx_a*By_a - My_a*Bx_a; 
  const double e_a = E_a - 0.5*(Mx_a*Mx_a + My_a*My_a + Mz_a*Mz_a)/rho_a; 
  const double rate_freq = sqrt(omega_p_sq[0]) + sqrt(Bx_a*Bx_a + By_a*By_a + Bz_a*Bz_a); 
  const double rate_accel = sqrt(2.0*(fx*fx + fy*fy + fz*fz)/(rho_a*e_a)); 
  return fmax(rate_freq, rate_accel); 

} 
