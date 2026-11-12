#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p3_exp_sq.h> 
GKYL_CU_DH void sr_vars_pressure_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
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
  const double volFact = dxv[2]*dxv[3]/4; 
 
  const double *p0 = &vmap[0]; 
  const double *V_0 = &u_i[0]; 
  const double *V_0_sq = &u_i_sq[0]; 
 
  const double *p1 = &vmap[4]; 
  const double *V_1 = &u_i[4]; 
  const double *V_1_sq = &u_i_sq[4]; 
 
  double temp[16] = {0.0}; 
  double temp_sq[16] = {0.0}; 
  double p_fac[16] = {0.0}; 
  double p0_sq[4] = {0.0}; 
  ser_1x_p3_exp_sq(p0, p0_sq); 
  double p1_sq[4] = {0.0}; 
  ser_1x_p3_exp_sq(p1, p1_sq); 
  temp[0] = 1.4142135623730951*V_1[0]*p1[0]+1.4142135623730951*V_0[0]*p0[0]; 
  temp[1] = 1.4142135623730951*p1[0]*V_1[1]+1.4142135623730951*p0[0]*V_0[1]; 
  temp[2] = 1.4142135623730951*p1[0]*V_1[2]+1.4142135623730951*p0[0]*V_0[2]; 
  temp[3] = 1.4142135623730951*V_0[0]*p0[1]; 
  temp[4] = 1.4142135623730951*V_1[0]*p1[1]; 
  temp[5] = 1.4142135623730951*p1[0]*V_1[3]+1.4142135623730951*p0[0]*V_0[3]; 
  temp[6] = 1.4142135623730951*V_0[1]*p0[1]; 
  temp[7] = 1.4142135623730951*p0[1]*V_0[2]; 
  temp[8] = 1.4142135623730951*V_1[1]*p1[1]; 
  temp[9] = 1.4142135623730951*p1[1]*V_1[2]; 
  temp[11] = 1.4142135623730951*p0[1]*V_0[3]; 
  temp[12] = 1.4142135623730951*p1[1]*V_1[3]; 

  temp_sq[0] = p0[0]*p1[0]*V_0[3]*V_1[3]+p0[0]*p1[0]*V_0[2]*V_1[2]+p0[0]*p1[0]*V_0[1]*V_1[1]+1.4142135623730951*V_1_sq[0]*p1_sq[0]+V_0[0]*V_1[0]*p0[0]*p1[0]+1.4142135623730951*V_0_sq[0]*p0_sq[0]; 
  temp_sq[1] = p0[0]*p1[0]*V_0[2]*V_1[3]+p0[0]*p1[0]*V_1[2]*V_0[3]+1.4142135623730951*p1_sq[0]*V_1_sq[1]+V_0[0]*p0[0]*p1[0]*V_1[1]+1.4142135623730951*p0_sq[0]*V_0_sq[1]+V_1[0]*p0[0]*p1[0]*V_0[1]; 
  temp_sq[2] = p0[0]*p1[0]*V_0[1]*V_1[3]+p0[0]*p1[0]*V_1[1]*V_0[3]+1.4142135623730951*p1_sq[0]*V_1_sq[2]+V_0[0]*p0[0]*p1[0]*V_1[2]+1.4142135623730951*p0_sq[0]*V_0_sq[2]+V_1[0]*p0[0]*p1[0]*V_0[2]; 
  temp_sq[3] = p1[0]*p0[1]*V_0[3]*V_1[3]+p1[0]*p0[1]*V_0[2]*V_1[2]+1.4142135623730951*V_0_sq[0]*p0_sq[1]+p1[0]*V_0[1]*V_1[1]*p0[1]+V_0[0]*V_1[0]*p1[0]*p0[1]; 
  temp_sq[4] = p0[0]*p1[1]*V_0[3]*V_1[3]+p0[0]*p1[1]*V_0[2]*V_1[2]+1.4142135623730951*V_1_sq[0]*p1_sq[1]+p0[0]*V_0[1]*V_1[1]*p1[1]+V_0[0]*V_1[0]*p0[0]*p1[1]; 
  temp_sq[5] = 1.4142135623730951*p1_sq[0]*V_1_sq[3]+V_0[0]*p0[0]*p1[0]*V_1[3]+1.4142135623730951*p0_sq[0]*V_0_sq[3]+V_1[0]*p0[0]*p1[0]*V_0[3]+p0[0]*p1[0]*V_0[1]*V_1[2]+p0[0]*p1[0]*V_1[1]*V_0[2]; 
  temp_sq[6] = p1[0]*p0[1]*V_0[2]*V_1[3]+p1[0]*p0[1]*V_1[2]*V_0[3]+1.4142135623730951*V_0_sq[1]*p0_sq[1]+V_0[0]*p1[0]*V_1[1]*p0[1]+V_1[0]*p1[0]*V_0[1]*p0[1]; 
  temp_sq[7] = p1[0]*V_0[1]*p0[1]*V_1[3]+p1[0]*V_1[1]*p0[1]*V_0[3]+V_0[0]*p1[0]*p0[1]*V_1[2]+1.4142135623730951*p0_sq[1]*V_0_sq[2]+V_1[0]*p1[0]*p0[1]*V_0[2]; 
  temp_sq[8] = p0[0]*p1[1]*V_0[2]*V_1[3]+p0[0]*p1[1]*V_1[2]*V_0[3]+1.4142135623730951*V_1_sq[1]*p1_sq[1]+V_0[0]*p0[0]*V_1[1]*p1[1]+V_1[0]*p0[0]*V_0[1]*p1[1]; 
  temp_sq[9] = p0[0]*V_0[1]*p1[1]*V_1[3]+p0[0]*V_1[1]*p1[1]*V_0[3]+1.4142135623730951*p1_sq[1]*V_1_sq[2]+V_0[0]*p0[0]*p1[1]*V_1[2]+V_1[0]*p0[0]*p1[1]*V_0[2]; 
  temp_sq[10] = p0[1]*p1[1]*V_0[3]*V_1[3]+p0[1]*p1[1]*V_0[2]*V_1[2]+V_0[1]*V_1[1]*p0[1]*p1[1]+V_0[0]*V_1[0]*p0[1]*p1[1]; 
  temp_sq[11] = V_0[0]*p1[0]*p0[1]*V_1[3]+1.4142135623730951*p0_sq[1]*V_0_sq[3]+V_1[0]*p1[0]*p0[1]*V_0[3]+p1[0]*V_0[1]*p0[1]*V_1[2]+p1[0]*V_1[1]*p0[1]*V_0[2]; 
  temp_sq[12] = 1.4142135623730951*p1_sq[1]*V_1_sq[3]+V_0[0]*p0[0]*p1[1]*V_1[3]+V_1[0]*p0[0]*p1[1]*V_0[3]+p0[0]*V_0[1]*p1[1]*V_1[2]+p0[0]*V_1[1]*p1[1]*V_0[2]; 
  temp_sq[13] = p0[1]*p1[1]*V_0[2]*V_1[3]+p0[1]*p1[1]*V_1[2]*V_0[3]+V_0[0]*V_1[1]*p0[1]*p1[1]+V_1[0]*V_0[1]*p0[1]*p1[1]; 
  temp_sq[14] = V_0[1]*p0[1]*p1[1]*V_1[3]+V_1[1]*p0[1]*p1[1]*V_0[3]+V_0[0]*p0[1]*p1[1]*V_1[2]+V_1[0]*p0[1]*p1[1]*V_0[2]; 
  temp_sq[15] = V_0[0]*p0[1]*p1[1]*V_1[3]+V_1[0]*p0[1]*p1[1]*V_0[3]+V_0[1]*p0[1]*p1[1]*V_1[2]+V_1[1]*p0[1]*p1[1]*V_0[2]; 

  p_fac[0] = 0.5*gamma_inv[3]*temp_sq[10]-1.0*GammaV[3]*temp[5]+0.5*gamma_inv[2]*temp_sq[4]+0.5*gamma_inv[1]*temp_sq[3]-1.0*GammaV[2]*temp[2]-1.0*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.5*gamma_inv[0]*temp_sq[0]-1.0*GammaV[0]*temp[0]-2.0*gamma_inv[0]; 
  p_fac[1] = 0.5*gamma_inv[3]*temp_sq[13]+0.5*gamma_inv[2]*temp_sq[8]+0.5*gamma_inv[1]*temp_sq[6]-1.0*GammaV[2]*temp[5]-1.0*temp[2]*GammaV[3]+0.5*gamma_inv[0]*temp_sq[1]-1.0*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.0*temp[0]*GammaV[1]; 
  p_fac[2] = 0.5*gamma_inv[3]*temp_sq[14]+0.5*gamma_inv[2]*temp_sq[9]+0.5*gamma_inv[1]*temp_sq[7]-1.0*GammaV[1]*temp[5]-1.0*temp[1]*GammaV[3]+0.5*gamma_inv[0]*temp_sq[2]-1.0*GammaV[0]*temp[2]+gamma[0]*GammaV_sq[2]-1.0*temp[0]*GammaV[2]; 
  p_fac[3] = -(1.0*GammaV[3]*temp[11])+0.5*gamma_inv[2]*temp_sq[10]-1.0*GammaV[2]*temp[7]-1.0*GammaV[1]*temp[6]+0.5*gamma_inv[3]*temp_sq[4]+0.5*gamma_inv[0]*temp_sq[3]-1.0*GammaV[0]*temp[3]+GammaV_sq[0]*gamma[1]+0.5*temp_sq[0]*gamma_inv[1]-2.0*gamma_inv[1]; 
  p_fac[4] = -(1.0*GammaV[3]*temp[12])+0.5*gamma_inv[1]*temp_sq[10]-1.0*GammaV[2]*temp[9]-1.0*GammaV[1]*temp[8]+0.5*gamma_inv[0]*temp_sq[4]-1.0*GammaV[0]*temp[4]+0.5*gamma_inv[3]*temp_sq[3]+GammaV_sq[0]*gamma[2]+0.5*temp_sq[0]*gamma_inv[2]-2.0*gamma_inv[2]; 
  p_fac[5] = 0.5*gamma_inv[3]*temp_sq[15]+0.5*gamma_inv[2]*temp_sq[12]+0.5*gamma_inv[1]*temp_sq[11]+0.5*gamma_inv[0]*temp_sq[5]-1.0*GammaV[0]*temp[5]+gamma[0]*GammaV_sq[3]-1.0*temp[0]*GammaV[3]-1.0*GammaV[1]*temp[2]-1.0*temp[1]*GammaV[2]; 
  p_fac[6] = 0.5*gamma_inv[2]*temp_sq[13]-1.0*GammaV[2]*temp[11]+0.5*gamma_inv[3]*temp_sq[8]-1.0*GammaV[3]*temp[7]+0.5*gamma_inv[0]*temp_sq[6]-1.0*GammaV[0]*temp[6]-1.0*GammaV[1]*temp[3]+GammaV_sq[1]*gamma[1]+0.5*gamma_inv[1]*temp_sq[1]; 
  p_fac[7] = 0.5*gamma_inv[2]*temp_sq[14]-1.0*GammaV[1]*temp[11]+0.5*gamma_inv[3]*temp_sq[9]+0.5*gamma_inv[0]*temp_sq[7]-1.0*GammaV[0]*temp[7]-1.0*GammaV[3]*temp[6]-1.0*GammaV[2]*temp[3]+0.5*gamma_inv[1]*temp_sq[2]+gamma[1]*GammaV_sq[2]; 
  p_fac[8] = 0.5*gamma_inv[1]*temp_sq[13]-1.0*GammaV[2]*temp[12]-1.0*GammaV[3]*temp[9]+0.5*gamma_inv[0]*temp_sq[8]-1.0*GammaV[0]*temp[8]+0.5*gamma_inv[3]*temp_sq[6]-1.0*GammaV[1]*temp[4]+GammaV_sq[1]*gamma[2]+0.5*temp_sq[1]*gamma_inv[2]; 
  p_fac[9] = 0.5*gamma_inv[1]*temp_sq[14]-1.0*GammaV[1]*temp[12]+0.5*gamma_inv[0]*temp_sq[9]-1.0*GammaV[0]*temp[9]-1.0*GammaV[3]*temp[8]+0.5*gamma_inv[3]*temp_sq[7]-1.0*GammaV[2]*temp[4]+GammaV_sq[2]*gamma[2]+0.5*gamma_inv[2]*temp_sq[2]; 
  p_fac[10] = 0.5*gamma_inv[0]*temp_sq[10]+0.5*gamma_inv[1]*temp_sq[4]+GammaV_sq[0]*gamma[3]+0.5*gamma_inv[2]*temp_sq[3]+0.5*temp_sq[0]*gamma_inv[3]-2.0*gamma_inv[3]; 
  p_fac[11] = 0.5*gamma_inv[2]*temp_sq[15]+0.5*gamma_inv[3]*temp_sq[12]+0.5*gamma_inv[0]*temp_sq[11]-1.0*GammaV[0]*temp[11]-1.0*GammaV[1]*temp[7]-1.0*GammaV[2]*temp[6]+0.5*gamma_inv[1]*temp_sq[5]-1.0*GammaV[3]*temp[3]+gamma[1]*GammaV_sq[3]; 
  p_fac[12] = 0.5*gamma_inv[1]*temp_sq[15]+0.5*gamma_inv[0]*temp_sq[12]-1.0*GammaV[0]*temp[12]+0.5*gamma_inv[3]*temp_sq[11]-1.0*GammaV[1]*temp[9]-1.0*GammaV[2]*temp[8]+0.5*gamma_inv[2]*temp_sq[5]-1.0*GammaV[3]*temp[4]+gamma[2]*GammaV_sq[3]; 
  p_fac[13] = 0.5*gamma_inv[0]*temp_sq[13]+0.5*gamma_inv[1]*temp_sq[8]+0.5*gamma_inv[2]*temp_sq[6]+GammaV_sq[1]*gamma[3]+0.5*temp_sq[1]*gamma_inv[3]; 
  p_fac[14] = 0.5*gamma_inv[0]*temp_sq[14]+0.5*gamma_inv[1]*temp_sq[9]+0.5*gamma_inv[2]*temp_sq[7]+GammaV_sq[2]*gamma[3]+0.5*temp_sq[2]*gamma_inv[3]; 
  p_fac[15] = 0.5*gamma_inv[0]*temp_sq[15]+0.5*gamma_inv[1]*temp_sq[12]+0.5*gamma_inv[2]*temp_sq[11]+0.5*gamma_inv[3]*temp_sq[5]+GammaV_sq[3]*gamma[3]; 

  sr_pressure[0] += (volFact/2)*(0.5*f[15]*p_fac[15]+0.5*f[14]*p_fac[14]+0.5*f[13]*p_fac[13]+0.5*f[12]*p_fac[12]+0.5*f[11]*p_fac[11]+0.5*f[10]*p_fac[10]+0.5*f[9]*p_fac[9]+0.5*f[8]*p_fac[8]+0.5*f[7]*p_fac[7]+0.5*f[6]*p_fac[6]+0.5*f[5]*p_fac[5]+0.5*f[4]*p_fac[4]+0.5*f[3]*p_fac[3]+0.5*f[2]*p_fac[2]+0.5*f[1]*p_fac[1]+0.5*f[0]*p_fac[0]); 
  sr_pressure[1] += (volFact/2)*(0.5*f[14]*p_fac[15]+0.5*p_fac[14]*f[15]+0.5*f[10]*p_fac[13]+0.5*p_fac[10]*f[13]+0.5*f[9]*p_fac[12]+0.5*p_fac[9]*f[12]+0.5*f[7]*p_fac[11]+0.5*p_fac[7]*f[11]+0.5*f[4]*p_fac[8]+0.5*p_fac[4]*f[8]+0.5*f[3]*p_fac[6]+0.5*p_fac[3]*f[6]+0.5*f[2]*p_fac[5]+0.5*p_fac[2]*f[5]+0.5*f[0]*p_fac[1]+0.5*p_fac[0]*f[1]); 
  sr_pressure[2] += (volFact/2)*(0.5*f[13]*p_fac[15]+0.5*p_fac[13]*f[15]+0.5*f[10]*p_fac[14]+0.5*p_fac[10]*f[14]+0.5*f[8]*p_fac[12]+0.5*p_fac[8]*f[12]+0.5*f[6]*p_fac[11]+0.5*p_fac[6]*f[11]+0.5*f[4]*p_fac[9]+0.5*p_fac[4]*f[9]+0.5*f[3]*p_fac[7]+0.5*p_fac[3]*f[7]+0.5*f[1]*p_fac[5]+0.5*p_fac[1]*f[5]+0.5*f[0]*p_fac[2]+0.5*p_fac[0]*f[2]); 
  sr_pressure[3] += (volFact/2)*(0.5*f[10]*p_fac[15]+0.5*p_fac[10]*f[15]+0.5*f[13]*p_fac[14]+0.5*p_fac[13]*f[14]+0.5*f[4]*p_fac[12]+0.5*p_fac[4]*f[12]+0.5*f[3]*p_fac[11]+0.5*p_fac[3]*f[11]+0.5*f[8]*p_fac[9]+0.5*p_fac[8]*f[9]+0.5*f[6]*p_fac[7]+0.5*p_fac[6]*f[7]+0.5*f[0]*p_fac[5]+0.5*p_fac[0]*f[5]+0.5*f[1]*p_fac[2]+0.5*p_fac[1]*f[2]); 

} 
