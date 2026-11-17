#include <gkyl_iso_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_iso_source_2x_ser_p1(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out) 
{ 
  // app_accel: External applied acceleration (external forces).
  // fluid:     [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  // out:       Output increment

  const double *app_accel_x = &app_accel[0]; 
  const double *app_accel_y = &app_accel[4]; 
  const double *app_accel_z = &app_accel[8]; 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[4]; 
  const double *rhouy = &fluid[8]; 
  const double *rhouz = &fluid[12]; 

  double *out_rhoux = &out[4]; 
  double *out_rhouy = &out[8]; 
  double *out_rhouz = &out[12]; 

  out_rhoux[0] += 0.5*app_accel_x[3]*rho[3]+0.5*app_accel_x[2]*rho[2]+0.5*app_accel_x[1]*rho[1]+0.5*app_accel_x[0]*rho[0]; 
  out_rhoux[1] += 0.5*app_accel_x[2]*rho[3]+0.5*rho[2]*app_accel_x[3]+0.5*app_accel_x[0]*rho[1]+0.5*rho[0]*app_accel_x[1]; 
  out_rhoux[2] += 0.5*app_accel_x[1]*rho[3]+0.5*rho[1]*app_accel_x[3]+0.5*app_accel_x[0]*rho[2]+0.5*rho[0]*app_accel_x[2]; 
  out_rhoux[3] += 0.5*app_accel_x[0]*rho[3]+0.5*rho[0]*app_accel_x[3]+0.5*app_accel_x[1]*rho[2]+0.5*rho[1]*app_accel_x[2]; 

  out_rhouy[0] += 0.5*app_accel_y[3]*rho[3]+0.5*app_accel_y[2]*rho[2]+0.5*app_accel_y[1]*rho[1]+0.5*app_accel_y[0]*rho[0]; 
  out_rhouy[1] += 0.5*app_accel_y[2]*rho[3]+0.5*rho[2]*app_accel_y[3]+0.5*app_accel_y[0]*rho[1]+0.5*rho[0]*app_accel_y[1]; 
  out_rhouy[2] += 0.5*app_accel_y[1]*rho[3]+0.5*rho[1]*app_accel_y[3]+0.5*app_accel_y[0]*rho[2]+0.5*rho[0]*app_accel_y[2]; 
  out_rhouy[3] += 0.5*app_accel_y[0]*rho[3]+0.5*rho[0]*app_accel_y[3]+0.5*app_accel_y[1]*rho[2]+0.5*rho[1]*app_accel_y[2]; 

  out_rhouz[0] += 0.5*app_accel_z[3]*rho[3]+0.5*app_accel_z[2]*rho[2]+0.5*app_accel_z[1]*rho[1]+0.5*app_accel_z[0]*rho[0]; 
  out_rhouz[1] += 0.5*app_accel_z[2]*rho[3]+0.5*rho[2]*app_accel_z[3]+0.5*app_accel_z[0]*rho[1]+0.5*rho[0]*app_accel_z[1]; 
  out_rhouz[2] += 0.5*app_accel_z[1]*rho[3]+0.5*rho[1]*app_accel_z[3]+0.5*app_accel_z[0]*rho[2]+0.5*rho[0]*app_accel_z[2]; 
  out_rhouz[3] += 0.5*app_accel_z[0]*rho[3]+0.5*rho[0]*app_accel_z[3]+0.5*app_accel_z[1]*rho[2]+0.5*rho[1]*app_accel_z[2]; 

} 
