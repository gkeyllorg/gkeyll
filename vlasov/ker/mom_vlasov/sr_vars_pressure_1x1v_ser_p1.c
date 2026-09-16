#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p3_exp_sq.h> 
GKYL_CU_DH void sr_vars_pressure_1x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
{ 
  // w:           Cell-center coordinates of velocity grid.
  // dxv:         Cell spacing of velocity grid.
  // vmap:        Momentum-space nonuniform mapping.
  // gamma:       Particle Lorentz boost factor sqrt(1 + p^2).
  // gamma_inv:   Inverse particle Lorentz boost factor 1/sqrt(1 + p^2).
  // u_i:         Spatial components of bulk four-velocity = GammaV*V_drift. 
  // u_i_sq:      Squared spatial components of bulk four-velocity = u_i^2. 
  // GammaV:      Bulk four-velocity Lorentz factor = sqrt(1 + |u_i|^2). 
  // GammaV_sq:   Squared bulk four-velocity Lorentz factor = 1 + |u_i|^2. 
  // f:           Input distribution function.
  // sr_pressure: Output relativistic pressure.
  const double volFact = dxv[1]/2; 
 
  const double *p0 = &vmap[0]; 
  const double *V_0 = &u_i[0]; 
  const double *V_0_sq = &u_i_sq[0]; 
 
  double temp[4] = {0.0}; 
  double temp_sq[4] = {0.0}; 
  double p_fac[4] = {0.0}; 
  double p0_sq[4] = {0.0}; 
  ser_1x_p3_exp_sq(p0, p0_sq); 
  temp[0] = V_0[0]*p0[0]; 
  temp[1] = p0[0]*V_0[1]; 
  temp[2] = V_0[0]*p0[1]; 
  temp[3] = V_0[1]*p0[1]; 

  temp_sq[0] = V_0_sq[0]*p0_sq[0]; 
  temp_sq[1] = p0_sq[0]*V_0_sq[1]; 
  temp_sq[2] = V_0_sq[0]*p0_sq[1]; 
  temp_sq[3] = V_0_sq[1]*p0_sq[1]; 

  p_fac[0] = 0.7071067811865475*gamma_inv[1]*temp_sq[2]-1.4142135623730951*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.7071067811865475*gamma_inv[0]*temp_sq[0]-1.4142135623730951*GammaV[0]*temp[0]-1.4142135623730951*gamma_inv[0]; 
  p_fac[1] = 0.7071067811865475*gamma_inv[1]*temp_sq[3]+0.7071067811865475*gamma_inv[0]*temp_sq[1]-1.4142135623730951*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.4142135623730951*temp[0]*GammaV[1]; 
  p_fac[2] = -(1.4142135623730951*GammaV[1]*temp[3])+0.7071067811865475*gamma_inv[0]*temp_sq[2]-1.4142135623730951*GammaV[0]*temp[2]+GammaV_sq[0]*gamma[1]+0.7071067811865475*temp_sq[0]*gamma_inv[1]-1.4142135623730951*gamma_inv[1]; 
  p_fac[3] = 0.7071067811865475*gamma_inv[0]*temp_sq[3]-1.4142135623730951*GammaV[0]*temp[3]-1.4142135623730951*GammaV[1]*temp[2]+GammaV_sq[1]*gamma[1]+0.7071067811865475*gamma_inv[1]*temp_sq[1]; 

  sr_pressure[0] += (volFact/1)*(0.7071067811865475*f[3]*p_fac[3]+0.7071067811865475*f[2]*p_fac[2]+0.7071067811865475*f[1]*p_fac[1]+0.7071067811865475*f[0]*p_fac[0]); 
  sr_pressure[1] += (volFact/1)*(0.7071067811865475*f[2]*p_fac[3]+0.7071067811865475*p_fac[2]*f[3]+0.7071067811865475*f[0]*p_fac[1]+0.7071067811865475*p_fac[0]*f[1]); 

} 
