#include <gkyl_euler_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
GKYL_CU_DH void fluid_vars_pressure_2x_tensor_p2(double param, const double *fluid, const double *u, 
    double* GKYL_RESTRICT p, double* GKYL_RESTRICT p_surf) 
{ 
  // param:  Adiabatic index gas_gamma.
  // fluid:  [rho, rho ux, rho uy, rho uz...], Fluid input state vector.
  // u:      Input volume expansion of flow velocity: [ux, uy, uz]. 
  // p:      Output volume expansion of pressure.
  // p_surf: Output surface expansion of pressure.
  //         [p_xl, p_xr, p_yl, p_yr, p_zl, p_zr] 

  const double *rho = &fluid[0]; 
  const double *rhoux = &fluid[9]; 
  const double *rhouy = &fluid[18]; 
  const double *rhouz = &fluid[27]; 
  const double *energy = &fluid[36]; 

  const double *ux = &u[0]; 
  const double *uy = &u[9]; 
  const double *uz = &u[18]; 

  double rhoux2[9] = {0.0}; 
  binop_mul_2d_tensor_p2(rhoux, ux, rhoux2); 
 
  double rhouy2[9] = {0.0}; 
  binop_mul_2d_tensor_p2(rhouy, uy, rhouy2); 
 
  double rhouz2[9] = {0.0}; 
  binop_mul_2d_tensor_p2(rhouz, uz, rhouz2); 
 
  p[0] = (param - 1.0)*(energy[0] - 0.5*(rhoux2[0] + rhouy2[0] + rhouz2[0])); 
  p[1] = (param - 1.0)*(energy[1] - 0.5*(rhoux2[1] + rhouy2[1] + rhouz2[1])); 
  p[2] = (param - 1.0)*(energy[2] - 0.5*(rhoux2[2] + rhouy2[2] + rhouz2[2])); 
  p[3] = (param - 1.0)*(energy[3] - 0.5*(rhoux2[3] + rhouy2[3] + rhouz2[3])); 
  p[4] = (param - 1.0)*(energy[4] - 0.5*(rhoux2[4] + rhouy2[4] + rhouz2[4])); 
  p[5] = (param - 1.0)*(energy[5] - 0.5*(rhoux2[5] + rhouy2[5] + rhouz2[5])); 
  p[6] = (param - 1.0)*(energy[6] - 0.5*(rhoux2[6] + rhouy2[6] + rhouz2[6])); 
  p[7] = (param - 1.0)*(energy[7] - 0.5*(rhoux2[7] + rhouy2[7] + rhouz2[7])); 
  p[8] = (param - 1.0)*(energy[8] - 0.5*(rhoux2[8] + rhouy2[8] + rhouz2[8])); 
  double *p_xl = &p_surf[0]; 
  double *p_xr = &p_surf[3]; 

  double *p_yl = &p_surf[6]; 
  double *p_yr = &p_surf[9]; 
 
  p_xl[0] = 1.58113883008419*p[4]-1.224744871391589*p[1]+0.7071067811865475*p[0]; 
  p_xl[1] = 1.58113883008419*p[6]-1.224744871391589*p[3]+0.7071067811865475*p[2]; 
  p_xl[2] = 1.58113883008419*p[8]-1.224744871391589*p[7]+0.7071067811865475*p[5]; 
 
  p_xr[0] = 1.58113883008419*p[4]+1.224744871391589*p[1]+0.7071067811865475*p[0]; 
  p_xr[1] = 1.58113883008419*p[6]+1.224744871391589*p[3]+0.7071067811865475*p[2]; 
  p_xr[2] = 1.58113883008419*p[8]+1.224744871391589*p[7]+0.7071067811865475*p[5]; 
 
  p_yl[0] = 1.58113883008419*p[5]-1.224744871391589*p[2]+0.7071067811865475*p[0]; 
  p_yl[1] = 1.58113883008419*p[7]-1.224744871391589*p[3]+0.7071067811865475*p[1]; 
  p_yl[2] = 1.58113883008419*p[8]-1.224744871391589*p[6]+0.7071067811865475*p[4]; 
 
  p_yr[0] = 1.58113883008419*p[5]+1.224744871391589*p[2]+0.7071067811865475*p[0]; 
  p_yr[1] = 1.58113883008419*p[7]+1.224744871391589*p[3]+0.7071067811865475*p[1]; 
  p_yr[2] = 1.58113883008419*p[8]+1.224744871391589*p[6]+0.7071067811865475*p[4]; 
 
} 
