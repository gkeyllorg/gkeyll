#include <gkyl_iso_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_iso_pressure_3x_ser_p1(double param, const double *fluid, const double *u, 
    double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf) 
{ 
  // param:  Thermal velocity vt of isothermal Euler.
  // fluid:  [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  // u:      Input volume expansion of flow velocity: [ux, uy, uz]. 
  // p:      Output volume expansion of pressure.
  // p_surf: Output surface expansion of pressure.
  //         [p_xl, p_xr, p_yl, p_yr, p_zl, p_zr] 

  const double *rho = &fluid[0]; 

  p[0] = param*param*rho[0]; 
  p[1] = param*param*rho[1]; 
  p[2] = param*param*rho[2]; 
  p[3] = param*param*rho[3]; 
  p[4] = param*param*rho[4]; 
  p[5] = param*param*rho[5]; 
  p[6] = param*param*rho[6]; 
  p[7] = param*param*rho[7]; 
  double *p_xl = &p_surf[0]; 
  double *p_xr = &p_surf[4]; 

  double *p_yl = &p_surf[8]; 
  double *p_yr = &p_surf[12]; 
 
  double *p_zl = &p_surf[16]; 
  double *p_zr = &p_surf[20]; 
 
  p_xl[0] = 0.7071067811865475*p[0]-1.224744871391589*p[1]; 
  p_xl[1] = 0.7071067811865475*p[2]-1.224744871391589*p[4]; 
  p_xl[2] = 0.7071067811865475*p[3]-1.224744871391589*p[5]; 
  p_xl[3] = 0.7071067811865475*p[6]-1.224744871391589*p[7]; 
 
  p_xr[0] = 1.224744871391589*p[1]+0.7071067811865475*p[0]; 
  p_xr[1] = 1.224744871391589*p[4]+0.7071067811865475*p[2]; 
  p_xr[2] = 1.224744871391589*p[5]+0.7071067811865475*p[3]; 
  p_xr[3] = 1.224744871391589*p[7]+0.7071067811865475*p[6]; 
 
  p_yl[0] = 0.7071067811865475*p[0]-1.224744871391589*p[2]; 
  p_yl[1] = 0.7071067811865475*p[1]-1.224744871391589*p[4]; 
  p_yl[2] = 0.7071067811865475*p[3]-1.224744871391589*p[6]; 
  p_yl[3] = 0.7071067811865475*p[5]-1.224744871391589*p[7]; 
 
  p_yr[0] = 1.224744871391589*p[2]+0.7071067811865475*p[0]; 
  p_yr[1] = 1.224744871391589*p[4]+0.7071067811865475*p[1]; 
  p_yr[2] = 1.224744871391589*p[6]+0.7071067811865475*p[3]; 
  p_yr[3] = 1.224744871391589*p[7]+0.7071067811865475*p[5]; 
 
  p_zl[0] = 0.7071067811865475*p[0]-1.224744871391589*p[3]; 
  p_zl[1] = 0.7071067811865475*p[1]-1.224744871391589*p[5]; 
  p_zl[2] = 0.7071067811865475*p[2]-1.224744871391589*p[6]; 
  p_zl[3] = 0.7071067811865475*p[4]-1.224744871391589*p[7]; 
 
  p_zr[0] = 1.224744871391589*p[3]+0.7071067811865475*p[0]; 
  p_zr[1] = 1.224744871391589*p[5]+0.7071067811865475*p[1]; 
  p_zr[2] = 1.224744871391589*p[6]+0.7071067811865475*p[2]; 
  p_zr[3] = 1.224744871391589*p[7]+0.7071067811865475*p[4]; 
 
} 
