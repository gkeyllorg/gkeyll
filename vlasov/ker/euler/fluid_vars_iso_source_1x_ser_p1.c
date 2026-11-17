#include <gkyl_iso_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_iso_source_1x_ser_p1(const double* app_accel, const double* fluid, double* GKYL_RESTRICT out) 
{ 
  // app_accel: External applied acceleration (external forces).
  // fluid:     [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  // out:       Output increment

  const double *app_accel_x = &app_accel[0]; 
  const double *app_accel_y = &app_accel[2]; 
  const double *app_accel_z = &app_accel[4]; 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[2]; 
  const double *rhouy = &fluid[4]; 
  const double *rhouz = &fluid[6]; 

  double *out_rhoux = &out[2]; 
  double *out_rhouy = &out[4]; 
  double *out_rhouz = &out[6]; 

  out_rhoux[0] += 0.7071067811865475*app_accel_x[1]*rho[1]+0.7071067811865475*app_accel_x[0]*rho[0]; 
  out_rhoux[1] += 0.7071067811865475*app_accel_x[0]*rho[1]+0.7071067811865475*rho[0]*app_accel_x[1]; 

  out_rhouy[0] += 0.7071067811865475*app_accel_y[1]*rho[1]+0.7071067811865475*app_accel_y[0]*rho[0]; 
  out_rhouy[1] += 0.7071067811865475*app_accel_y[0]*rho[1]+0.7071067811865475*rho[0]*app_accel_y[1]; 

  out_rhouz[0] += 0.7071067811865475*app_accel_z[1]*rho[1]+0.7071067811865475*app_accel_z[0]*rho[0]; 
  out_rhouz[1] += 0.7071067811865475*app_accel_z[0]*rho[1]+0.7071067811865475*rho[0]*app_accel_z[1]; 

} 
