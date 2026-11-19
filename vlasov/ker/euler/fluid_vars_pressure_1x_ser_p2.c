#include <gkyl_euler_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void fluid_vars_pressure_1x_ser_p2(double param, const double *fluid, const double *u, 
    double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf) 
{ 
  // param:  Adiabatic index gas_gamma.
  // fluid:  [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  // u:      Input volume expansion of flow velocity: [ux, uy, uz]. 
  // p:      Output volume expansion of pressure.
  // p_surf: Output surface expansion of pressure.
  //         [p_xl, p_xr, p_yl, p_yr, p_zl, p_zr] 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[3]; 
  const double *rhouy = &fluid[6]; 
  const double *rhouz = &fluid[9]; 
  const double *energy = &fluid[12]; 

  const double *ux = &u[0]; 
  const double *uy = &u[3]; 
  const double *uz = &u[6]; 

  double rhoux2[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhoux, ux, rhoux2); 
 
  double rhouy2[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhouy, uy, rhouy2); 
 
  double rhouz2[3] = {0.0}; 
  binop_mul_1d_ser_p2(rhouz, uz, rhouz2); 
 
  p[0] = (param - 1.0)*(energy[0] - 0.5*(rhoux2[0] + rhouy2[0] + rhouz2[0])); 
  p[1] = (param - 1.0)*(energy[1] - 0.5*(rhoux2[1] + rhouy2[1] + rhouz2[1])); 
  p[2] = (param - 1.0)*(energy[2] - 0.5*(rhoux2[2] + rhouy2[2] + rhouz2[2])); 
  double *p_xl = &p_surf[0]; 
  double *p_xr = &p_surf[1]; 

  p_xl[0] = 1.58113883008419*p[2]-1.224744871391589*p[1]+0.7071067811865475*p[0]; 
  p_xr[0] = 1.58113883008419*p[2]+1.224744871391589*p[1]+0.7071067811865475*p[0]; 
 
} 
