#include <gkyl_iso_euler_kernels.h> 
GKYL_CU_DH void fluid_vars_iso_pressure_1x_ser_p2(double param, const double *fluid, const double *u, 
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
  double *p_xl = &p_surf[0]; 
  double *p_xr = &p_surf[1]; 

  p_xl[0] = 1.5811388300841895*p[2]-1.224744871391589*p[1]+0.7071067811865475*p[0]; 
  p_xr[0] = 1.5811388300841895*p[2]+1.224744871391589*p[1]+0.7071067811865475*p[0]; 
 
} 
