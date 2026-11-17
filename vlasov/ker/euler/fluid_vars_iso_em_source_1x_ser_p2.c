#include <gkyl_iso_euler_kernels.h> 
#include <math.h> 
GKYL_CU_DH double fluid_vars_iso_em_source_1x_ser_p2(const double* qmem, const double* omega_p_sq, const double* fluid, double* GKYL_RESTRICT out) 
{ 
  // qmem:       [ax, ay, az, q/m Bx, q/m By, q/m Bz], acceleration (q/m E plus any applied
  //             acceleration) and q/m times the magnetic field.
  // omega_p_sq: Plasma frequency squared of the cell (one value; all the coupled species).
  // fluid:      [rho, rho ux, rho uy, rho uz], Fluid input state vector.
  // out:        Output increment
  // Returns the rate of the source for the time-step (from the cell averages).

  const double *ax = &qmem[0]; 
  const double *ay = &qmem[3]; 
  const double *az = &qmem[6]; 
  const double *Bx = &qmem[9]; 
  const double *By = &qmem[12]; 
  const double *Bz = &qmem[15]; 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[3]; 
  const double *rhouy = &fluid[6]; 
  const double *rhouz = &fluid[9]; 

  double *out_rhoux = &out[3]; 
  double *out_rhouy = &out[6]; 
  double *out_rhouz = &out[9]; 

  out_rhoux[0] += -(0.7071067811865475*By[2]*rhouz[2])+0.7071067811865475*Bz[2]*rhouy[2]+0.7071067811865475*ax[2]*rho[2]-0.7071067811865475*By[1]*rhouz[1]+0.7071067811865475*Bz[1]*rhouy[1]+0.7071067811865475*ax[1]*rho[1]-0.7071067811865475*By[0]*rhouz[0]+0.7071067811865475*Bz[0]*rhouy[0]+0.7071067811865475*ax[0]*rho[0]; 
  out_rhoux[1] += -(0.6324555320336759*By[1]*rhouz[2])+0.6324555320336759*Bz[1]*rhouy[2]+0.6324555320336759*ax[1]*rho[2]+0.6324555320336759*rho[1]*ax[2]+0.6324555320336759*rhouy[1]*Bz[2]-0.6324555320336759*rhouz[1]*By[2]-0.7071067811865475*By[0]*rhouz[1]+0.7071067811865475*Bz[0]*rhouy[1]+0.7071067811865475*ax[0]*rho[1]+0.7071067811865475*rho[0]*ax[1]+0.7071067811865475*rhouy[0]*Bz[1]-0.7071067811865475*rhouz[0]*By[1]; 
  out_rhoux[2] += -(0.45175395145262565*By[2]*rhouz[2])-0.7071067811865475*By[0]*rhouz[2]+0.45175395145262565*Bz[2]*rhouy[2]+0.7071067811865475*Bz[0]*rhouy[2]+0.45175395145262565*ax[2]*rho[2]+0.7071067811865475*ax[0]*rho[2]+0.7071067811865475*rho[0]*ax[2]+0.7071067811865475*rhouy[0]*Bz[2]-0.7071067811865475*rhouz[0]*By[2]-0.6324555320336759*By[1]*rhouz[1]+0.6324555320336759*Bz[1]*rhouy[1]+0.6324555320336759*ax[1]*rho[1]; 

  out_rhouy[0] += 0.7071067811865475*Bx[2]*rhouz[2]-0.7071067811865475*Bz[2]*rhoux[2]+0.7071067811865475*ay[2]*rho[2]+0.7071067811865475*Bx[1]*rhouz[1]-0.7071067811865475*Bz[1]*rhoux[1]+0.7071067811865475*ay[1]*rho[1]+0.7071067811865475*Bx[0]*rhouz[0]-0.7071067811865475*Bz[0]*rhoux[0]+0.7071067811865475*ay[0]*rho[0]; 
  out_rhouy[1] += 0.6324555320336759*Bx[1]*rhouz[2]-0.6324555320336759*Bz[1]*rhoux[2]+0.6324555320336759*ay[1]*rho[2]+0.6324555320336759*rho[1]*ay[2]-0.6324555320336759*rhoux[1]*Bz[2]+0.6324555320336759*rhouz[1]*Bx[2]+0.7071067811865475*Bx[0]*rhouz[1]-0.7071067811865475*Bz[0]*rhoux[1]+0.7071067811865475*ay[0]*rho[1]+0.7071067811865475*rho[0]*ay[1]-0.7071067811865475*rhoux[0]*Bz[1]+0.7071067811865475*rhouz[0]*Bx[1]; 
  out_rhouy[2] += 0.45175395145262565*Bx[2]*rhouz[2]+0.7071067811865475*Bx[0]*rhouz[2]-0.45175395145262565*Bz[2]*rhoux[2]-0.7071067811865475*Bz[0]*rhoux[2]+0.45175395145262565*ay[2]*rho[2]+0.7071067811865475*ay[0]*rho[2]+0.7071067811865475*rho[0]*ay[2]-0.7071067811865475*rhoux[0]*Bz[2]+0.7071067811865475*rhouz[0]*Bx[2]+0.6324555320336759*Bx[1]*rhouz[1]-0.6324555320336759*Bz[1]*rhoux[1]+0.6324555320336759*ay[1]*rho[1]; 

  out_rhouz[0] += -(0.7071067811865475*Bx[2]*rhouy[2])+0.7071067811865475*By[2]*rhoux[2]+0.7071067811865475*az[2]*rho[2]-0.7071067811865475*Bx[1]*rhouy[1]+0.7071067811865475*By[1]*rhoux[1]+0.7071067811865475*az[1]*rho[1]-0.7071067811865475*Bx[0]*rhouy[0]+0.7071067811865475*By[0]*rhoux[0]+0.7071067811865475*az[0]*rho[0]; 
  out_rhouz[1] += -(0.6324555320336759*Bx[1]*rhouy[2])+0.6324555320336759*By[1]*rhoux[2]+0.6324555320336759*az[1]*rho[2]+0.6324555320336759*rho[1]*az[2]+0.6324555320336759*rhoux[1]*By[2]-0.6324555320336759*rhouy[1]*Bx[2]-0.7071067811865475*Bx[0]*rhouy[1]+0.7071067811865475*By[0]*rhoux[1]+0.7071067811865475*az[0]*rho[1]+0.7071067811865475*rho[0]*az[1]+0.7071067811865475*rhoux[0]*By[1]-0.7071067811865475*rhouy[0]*Bx[1]; 
  out_rhouz[2] += -(0.45175395145262565*Bx[2]*rhouy[2])-0.7071067811865475*Bx[0]*rhouy[2]+0.45175395145262565*By[2]*rhoux[2]+0.7071067811865475*By[0]*rhoux[2]+0.45175395145262565*az[2]*rho[2]+0.7071067811865475*az[0]*rho[2]+0.7071067811865475*rho[0]*az[2]+0.7071067811865475*rhoux[0]*By[2]-0.7071067811865475*rhouy[0]*Bx[2]-0.6324555320336759*Bx[1]*rhouy[1]+0.6324555320336759*By[1]*rhoux[1]+0.6324555320336759*az[1]*rho[1]; 

  const double avg = 0.7071067811865475; 
  const double Bx_a = avg*Bx[0], By_a = avg*By[0], Bz_a = avg*Bz[0]; 
  return sqrt(omega_p_sq[0]) + sqrt(Bx_a*Bx_a + By_a*By_a + Bz_a*Bz_a); 

} 
