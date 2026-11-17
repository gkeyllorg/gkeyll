#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p3_exp_sq.h> 
GKYL_CU_DH void sr_vars_pressure_1x2v_tensor_p1(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
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
  const double volFact = dxv[1]*dxv[2]/4; 
 
  const double *p0 = &vmap[0]; 
  const double *V_0 = &u_i[0]; 
  const double *V_0_sq = &u_i_sq[0]; 
 
  const double *p1 = &vmap[4]; 
  const double *V_1 = &u_i[2]; 
  const double *V_1_sq = &u_i_sq[2]; 
 
  double temp[18] = {0.0}; 
  double temp_sq[18] = {0.0}; 
  double p_fac[18] = {0.0}; 
  double p0_sq[4] = {0.0}; 
  ser_1x_p3_exp_sq(p0, p0_sq); 
  double p1_sq[4] = {0.0}; 
  ser_1x_p3_exp_sq(p1, p1_sq); 
  temp[0] = 1.4142135623730951*V_1[0]*p1[0]+1.4142135623730951*V_0[0]*p0[0]; 
  temp[1] = 1.4142135623730951*p1[0]*V_1[1]+1.4142135623730951*p0[0]*V_0[1]; 
  temp[2] = 1.4142135623730951*V_0[0]*p0[1]; 
  temp[3] = 1.4142135623730951*V_1[0]*p1[1]; 
  temp[4] = 1.4142135623730951*V_0[1]*p0[1]; 
  temp[5] = 1.4142135623730951*V_1[1]*p1[1]; 
  temp[7] = 1.4142135623730951*V_0[0]*p0[2]; 
  temp[8] = 1.4142135623730951*V_1[0]*p1[2]; 
  temp[10] = 1.4142135623730951*V_0[1]*p0[2]; 
  temp[12] = 1.4142135623730951*V_1[1]*p1[2]; 

  temp_sq[0] = 1.4142135623730951*p0[0]*p1[0]*V_0[1]*V_1[1]+1.4142135623730951*V_1_sq[0]*p1_sq[0]+1.4142135623730951*V_0[0]*V_1[0]*p0[0]*p1[0]+1.4142135623730951*V_0_sq[0]*p0_sq[0]; 
  temp_sq[1] = 1.4142135623730951*p1_sq[0]*V_1_sq[1]+1.4142135623730951*V_0[0]*p0[0]*p1[0]*V_1[1]+1.4142135623730951*p0_sq[0]*V_0_sq[1]+1.4142135623730951*V_1[0]*p0[0]*p1[0]*V_0[1]; 
  temp_sq[2] = 1.4142135623730951*V_0_sq[0]*p0_sq[1]+1.4142135623730951*p1[0]*V_0[1]*V_1[1]*p0[1]+1.4142135623730951*V_0[0]*V_1[0]*p1[0]*p0[1]; 
  temp_sq[3] = 1.4142135623730951*V_1_sq[0]*p1_sq[1]+1.4142135623730951*p0[0]*V_0[1]*V_1[1]*p1[1]+1.4142135623730951*V_0[0]*V_1[0]*p0[0]*p1[1]; 
  temp_sq[4] = 1.4142135623730951*V_0_sq[1]*p0_sq[1]+1.4142135623730951*V_0[0]*p1[0]*V_1[1]*p0[1]+1.4142135623730951*V_1[0]*p1[0]*V_0[1]*p0[1]; 
  temp_sq[5] = 1.4142135623730951*V_1_sq[1]*p1_sq[1]+1.4142135623730951*V_0[0]*p0[0]*V_1[1]*p1[1]+1.4142135623730951*V_1[0]*p0[0]*V_0[1]*p1[1]; 
  temp_sq[6] = 1.4142135623730951*V_0[1]*V_1[1]*p0[1]*p1[1]+1.4142135623730951*V_0[0]*V_1[0]*p0[1]*p1[1]; 
  temp_sq[7] = 1.4142135623730951*V_0_sq[0]*p0_sq[2]+1.4142135623730951*p1[0]*V_0[1]*V_1[1]*p0[2]+1.4142135623730951*V_0[0]*V_1[0]*p1[0]*p0[2]; 
  temp_sq[8] = 1.4142135623730951*V_1_sq[0]*p1_sq[2]+1.4142135623730951*p0[0]*V_0[1]*V_1[1]*p1[2]+1.4142135623730951*V_0[0]*V_1[0]*p0[0]*p1[2]; 
  temp_sq[9] = 1.4142135623730951*V_0[0]*V_1[1]*p0[1]*p1[1]+1.4142135623730951*V_1[0]*V_0[1]*p0[1]*p1[1]; 
  temp_sq[10] = 1.4142135623730951*V_0_sq[1]*p0_sq[2]+1.4142135623730951*V_0[0]*p1[0]*V_1[1]*p0[2]+1.4142135623730951*V_1[0]*p1[0]*V_0[1]*p0[2]; 
  temp_sq[11] = 1.4142135623730951*V_0[1]*V_1[1]*p1[1]*p0[2]+1.4142135623730951*V_0[0]*V_1[0]*p1[1]*p0[2]; 
  temp_sq[12] = 1.4142135623730951*V_1_sq[1]*p1_sq[2]+1.4142135623730951*V_0[0]*p0[0]*V_1[1]*p1[2]+1.4142135623730951*V_1[0]*p0[0]*V_0[1]*p1[2]; 
  temp_sq[13] = 1.4142135623730951*V_0[1]*V_1[1]*p0[1]*p1[2]+1.4142135623730951*V_0[0]*V_1[0]*p0[1]*p1[2]; 
  temp_sq[14] = 1.4142135623730951*V_0[0]*V_1[1]*p1[1]*p0[2]+1.4142135623730951*V_1[0]*V_0[1]*p1[1]*p0[2]; 
  temp_sq[15] = 1.4142135623730951*V_0[0]*V_1[1]*p0[1]*p1[2]+1.4142135623730951*V_1[0]*V_0[1]*p0[1]*p1[2]; 
  temp_sq[16] = 1.4142135623730951*V_0[1]*V_1[1]*p0[2]*p1[2]+1.4142135623730951*V_0[0]*V_1[0]*p0[2]*p1[2]; 
  temp_sq[17] = 1.4142135623730951*V_0[0]*V_1[1]*p0[2]*p1[2]+1.4142135623730951*V_1[0]*V_0[1]*p0[2]*p1[2]; 

  p_fac[0] = 0.5*gamma_inv[8]*temp_sq[16]+0.5*gamma_inv[7]*temp_sq[13]+0.5*gamma_inv[6]*temp_sq[11]+0.5*gamma_inv[5]*temp_sq[8]+0.5*gamma_inv[4]*temp_sq[7]+0.5*gamma_inv[3]*temp_sq[6]+0.5*gamma_inv[2]*temp_sq[3]+0.5*gamma_inv[1]*temp_sq[2]-1.4142135623730951*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.5*gamma_inv[0]*temp_sq[0]-1.4142135623730951*GammaV[0]*temp[0]-1.4142135623730951*gamma_inv[0]; 
  p_fac[1] = 0.5*gamma_inv[8]*temp_sq[17]+0.5000000000000001*gamma_inv[7]*temp_sq[15]+0.5000000000000001*gamma_inv[6]*temp_sq[14]+0.5000000000000001*gamma_inv[5]*temp_sq[12]+0.5000000000000001*gamma_inv[4]*temp_sq[10]+0.5*gamma_inv[3]*temp_sq[9]+0.5*gamma_inv[2]*temp_sq[5]+0.5*gamma_inv[1]*temp_sq[4]+0.5*gamma_inv[0]*temp_sq[1]-1.4142135623730951*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.4142135623730951*temp[0]*GammaV[1]; 
  p_fac[2] = 0.44721359549995804*gamma_inv[7]*temp_sq[16]+0.44721359549995804*gamma_inv[8]*temp_sq[13]+0.5000000000000001*gamma_inv[5]*temp_sq[13]+0.44721359549995804*gamma_inv[3]*temp_sq[11]+0.5000000000000001*gamma_inv[7]*temp_sq[8]+0.4472135954999579*gamma_inv[1]*temp_sq[7]+0.44721359549995804*gamma_inv[6]*temp_sq[6]+0.5*gamma_inv[2]*temp_sq[6]-1.4142135623730951*GammaV[1]*temp[4]+0.4472135954999579*temp_sq[2]*gamma_inv[4]+0.5*gamma_inv[3]*temp_sq[3]+0.5*gamma_inv[0]*temp_sq[2]-1.4142135623730951*GammaV[0]*temp[2]+GammaV_sq[0]*gamma[1]+0.5*temp_sq[0]*gamma_inv[1]-1.4142135623730951*gamma_inv[1]; 
  p_fac[3] = 0.44721359549995804*gamma_inv[6]*temp_sq[16]+0.44721359549995804*gamma_inv[3]*temp_sq[13]+0.44721359549995804*gamma_inv[8]*temp_sq[11]+0.5000000000000001*gamma_inv[4]*temp_sq[11]+0.4472135954999579*gamma_inv[2]*temp_sq[8]+0.5000000000000001*gamma_inv[6]*temp_sq[7]+0.44721359549995804*temp_sq[6]*gamma_inv[7]+0.5*gamma_inv[1]*temp_sq[6]-1.4142135623730951*GammaV[1]*temp[5]+0.4472135954999579*temp_sq[3]*gamma_inv[5]+0.5*gamma_inv[0]*temp_sq[3]-1.4142135623730951*GammaV[0]*temp[3]+0.5*temp_sq[2]*gamma_inv[3]+GammaV_sq[0]*gamma[2]+0.5*temp_sq[0]*gamma_inv[2]-1.4142135623730951*gamma_inv[2]; 
  p_fac[4] = 0.44721359549995804*gamma_inv[7]*temp_sq[17]+0.4472135954999579*gamma_inv[8]*temp_sq[15]+0.5*gamma_inv[5]*temp_sq[15]+0.4472135954999579*gamma_inv[3]*temp_sq[14]+0.5*gamma_inv[7]*temp_sq[12]+0.44721359549995804*gamma_inv[1]*temp_sq[10]+0.44721359549995804*gamma_inv[6]*temp_sq[9]+0.5*gamma_inv[2]*temp_sq[9]+0.5*gamma_inv[3]*temp_sq[5]+0.4472135954999579*gamma_inv[4]*temp_sq[4]+0.5*gamma_inv[0]*temp_sq[4]-1.4142135623730951*GammaV[0]*temp[4]-1.4142135623730951*GammaV[1]*temp[2]+GammaV_sq[1]*gamma[1]+0.5*gamma_inv[1]*temp_sq[1]; 
  p_fac[5] = 0.44721359549995804*gamma_inv[6]*temp_sq[17]+0.4472135954999579*gamma_inv[3]*temp_sq[15]+0.4472135954999579*gamma_inv[8]*temp_sq[14]+0.5*gamma_inv[4]*temp_sq[14]+0.44721359549995804*gamma_inv[2]*temp_sq[12]+0.5*gamma_inv[6]*temp_sq[10]+0.44721359549995804*gamma_inv[7]*temp_sq[9]+0.5*gamma_inv[1]*temp_sq[9]+0.4472135954999579*gamma_inv[5]*temp_sq[5]+0.5*gamma_inv[0]*temp_sq[5]-1.4142135623730951*GammaV[0]*temp[5]+0.5*gamma_inv[3]*temp_sq[4]-1.4142135623730951*GammaV[1]*temp[3]+GammaV_sq[1]*gamma[2]+0.5*temp_sq[1]*gamma_inv[2]; 
  p_fac[6] = 0.4*gamma_inv[3]*temp_sq[16]+0.4*gamma_inv[6]*temp_sq[13]+0.44721359549995804*gamma_inv[2]*temp_sq[13]+0.4*gamma_inv[7]*temp_sq[11]+0.44721359549995804*gamma_inv[1]*temp_sq[11]+0.4472135954999579*gamma_inv[3]*temp_sq[8]+0.4*temp_sq[6]*gamma_inv[8]+0.4472135954999579*gamma_inv[3]*temp_sq[7]+0.44721359549995804*temp_sq[3]*gamma_inv[7]+0.4472135954999579*gamma_inv[5]*temp_sq[6]+0.4472135954999579*gamma_inv[4]*temp_sq[6]+0.5*gamma_inv[0]*temp_sq[6]+0.44721359549995804*temp_sq[2]*gamma_inv[6]+GammaV_sq[0]*gamma[3]+0.5*gamma_inv[1]*temp_sq[3]+0.5*temp_sq[0]*gamma_inv[3]-1.4142135623730951*gamma_inv[3]+0.5*gamma_inv[2]*temp_sq[2]; 
  p_fac[7] = 0.31943828249996997*gamma_inv[8]*temp_sq[16]+0.5*gamma_inv[5]*temp_sq[16]+0.4472135954999579*gamma_inv[7]*temp_sq[13]+0.31943828249996997*gamma_inv[6]*temp_sq[11]+0.5000000000000001*gamma_inv[2]*temp_sq[11]-1.4142135623730951*GammaV[1]*temp[10]+0.5*gamma_inv[8]*temp_sq[8]+0.31943828249996997*gamma_inv[4]*temp_sq[7]+0.5*gamma_inv[0]*temp_sq[7]-1.4142135623730951*GammaV[0]*temp[7]+0.4472135954999579*gamma_inv[3]*temp_sq[6]+0.5000000000000001*temp_sq[3]*gamma_inv[6]+GammaV_sq[0]*gamma[4]+0.5*temp_sq[0]*gamma_inv[4]-1.4142135623730951*gamma_inv[4]+0.4472135954999579*gamma_inv[1]*temp_sq[2]; 
  p_fac[8] = 0.31943828249996997*gamma_inv[8]*temp_sq[16]+0.5*gamma_inv[4]*temp_sq[16]+0.31943828249996997*gamma_inv[7]*temp_sq[13]+0.5000000000000001*gamma_inv[1]*temp_sq[13]-1.4142135623730951*GammaV[1]*temp[12]+0.4472135954999579*gamma_inv[6]*temp_sq[11]+0.31943828249996997*gamma_inv[5]*temp_sq[8]+0.5*gamma_inv[0]*temp_sq[8]-1.4142135623730951*GammaV[0]*temp[8]+0.5*temp_sq[7]*gamma_inv[8]+0.5000000000000001*temp_sq[2]*gamma_inv[7]+0.4472135954999579*gamma_inv[3]*temp_sq[6]+GammaV_sq[0]*gamma[5]+0.5*temp_sq[0]*gamma_inv[5]-1.4142135623730951*gamma_inv[5]+0.4472135954999579*gamma_inv[2]*temp_sq[3]; 
  p_fac[9] = 0.4*gamma_inv[3]*temp_sq[17]+0.4*gamma_inv[6]*temp_sq[15]+0.4472135954999579*gamma_inv[2]*temp_sq[15]+0.4*gamma_inv[7]*temp_sq[14]+0.4472135954999579*gamma_inv[1]*temp_sq[14]+0.44721359549995804*gamma_inv[3]*temp_sq[12]+0.44721359549995804*gamma_inv[3]*temp_sq[10]+0.4*gamma_inv[8]*temp_sq[9]+0.4472135954999579*gamma_inv[5]*temp_sq[9]+0.4472135954999579*gamma_inv[4]*temp_sq[9]+0.5*gamma_inv[0]*temp_sq[9]+0.44721359549995804*temp_sq[5]*gamma_inv[7]+0.44721359549995804*temp_sq[4]*gamma_inv[6]+0.5*gamma_inv[1]*temp_sq[5]+0.5*gamma_inv[2]*temp_sq[4]+GammaV_sq[1]*gamma[3]+0.5*temp_sq[1]*gamma_inv[3]; 
  p_fac[10] = 0.31943828249997*gamma_inv[8]*temp_sq[17]+0.5000000000000001*gamma_inv[5]*temp_sq[17]+0.4472135954999579*gamma_inv[7]*temp_sq[15]+0.31943828249996997*gamma_inv[6]*temp_sq[14]+0.5000000000000001*gamma_inv[2]*temp_sq[14]+0.5*gamma_inv[8]*temp_sq[12]+0.31943828249996997*gamma_inv[4]*temp_sq[10]+0.5*gamma_inv[0]*temp_sq[10]-1.4142135623730951*GammaV[0]*temp[10]+0.44721359549995804*gamma_inv[3]*temp_sq[9]-1.4142135623730951*GammaV[1]*temp[7]+0.5*temp_sq[5]*gamma_inv[6]+1.0000000000000002*GammaV_sq[1]*gamma[4]+0.44721359549995804*gamma_inv[1]*temp_sq[4]+0.5000000000000001*temp_sq[1]*gamma_inv[4]; 
  p_fac[11] = 0.2857142857142857*gamma_inv[6]*temp_sq[16]+0.44721359549995804*gamma_inv[2]*temp_sq[16]+0.4*gamma_inv[3]*temp_sq[13]+0.2857142857142857*gamma_inv[8]*temp_sq[11]+0.4472135954999579*gamma_inv[5]*temp_sq[11]+0.31943828249996997*gamma_inv[4]*temp_sq[11]+0.5*gamma_inv[0]*temp_sq[11]+0.4472135954999579*gamma_inv[6]*temp_sq[8]+0.44721359549995804*temp_sq[3]*gamma_inv[8]+0.31943828249996997*gamma_inv[6]*temp_sq[7]+0.5000000000000001*gamma_inv[2]*temp_sq[7]+0.4*temp_sq[6]*gamma_inv[7]+GammaV_sq[0]*gamma[6]+0.44721359549995804*gamma_inv[1]*temp_sq[6]+0.5*temp_sq[0]*gamma_inv[6]-1.4142135623730951*gamma_inv[6]+0.5000000000000001*temp_sq[3]*gamma_inv[4]+0.44721359549995804*temp_sq[2]*gamma_inv[3]; 
  p_fac[12] = 0.31943828249997*gamma_inv[8]*temp_sq[17]+0.5000000000000001*gamma_inv[4]*temp_sq[17]+0.31943828249996997*gamma_inv[7]*temp_sq[15]+0.5000000000000001*gamma_inv[1]*temp_sq[15]+0.4472135954999579*gamma_inv[6]*temp_sq[14]+0.31943828249996997*gamma_inv[5]*temp_sq[12]+0.5*gamma_inv[0]*temp_sq[12]-1.4142135623730951*GammaV[0]*temp[12]+0.5*gamma_inv[8]*temp_sq[10]+0.44721359549995804*gamma_inv[3]*temp_sq[9]-1.4142135623730951*GammaV[1]*temp[8]+0.5*temp_sq[4]*gamma_inv[7]+1.0000000000000002*GammaV_sq[1]*gamma[5]+0.44721359549995804*gamma_inv[2]*temp_sq[5]+0.5000000000000001*temp_sq[1]*gamma_inv[5]; 
  p_fac[13] = 0.2857142857142857*gamma_inv[7]*temp_sq[16]+0.44721359549995804*gamma_inv[1]*temp_sq[16]+0.2857142857142857*gamma_inv[8]*temp_sq[13]+0.31943828249996997*gamma_inv[5]*temp_sq[13]+0.4472135954999579*gamma_inv[4]*temp_sq[13]+0.5*gamma_inv[0]*temp_sq[13]+0.4*gamma_inv[3]*temp_sq[11]+0.31943828249996997*gamma_inv[7]*temp_sq[8]+0.5000000000000001*gamma_inv[1]*temp_sq[8]+0.44721359549995804*temp_sq[2]*gamma_inv[8]+GammaV_sq[0]*gamma[7]+0.4472135954999579*gamma_inv[7]*temp_sq[7]+0.5*temp_sq[0]*gamma_inv[7]-1.4142135623730951*gamma_inv[7]+0.4*gamma_inv[6]*temp_sq[6]+0.44721359549995804*gamma_inv[2]*temp_sq[6]+0.5000000000000001*temp_sq[2]*gamma_inv[5]+0.44721359549995804*gamma_inv[3]*temp_sq[3]; 
  p_fac[14] = 0.28571428571428575*gamma_inv[6]*temp_sq[17]+0.4472135954999579*gamma_inv[2]*temp_sq[17]+0.4*gamma_inv[3]*temp_sq[15]+0.2857142857142857*gamma_inv[8]*temp_sq[14]+0.4472135954999579*gamma_inv[5]*temp_sq[14]+0.31943828249996997*gamma_inv[4]*temp_sq[14]+0.5*gamma_inv[0]*temp_sq[14]+0.4472135954999579*gamma_inv[6]*temp_sq[12]+0.31943828249996997*gamma_inv[6]*temp_sq[10]+0.5000000000000001*gamma_inv[2]*temp_sq[10]+0.4*gamma_inv[7]*temp_sq[9]+0.4472135954999579*gamma_inv[1]*temp_sq[9]+0.4472135954999579*temp_sq[5]*gamma_inv[8]+1.0000000000000002*GammaV_sq[1]*gamma[6]+0.5000000000000001*temp_sq[1]*gamma_inv[6]+0.5*gamma_inv[4]*temp_sq[5]+0.4472135954999579*gamma_inv[3]*temp_sq[4]; 
  p_fac[15] = 0.28571428571428575*gamma_inv[7]*temp_sq[17]+0.4472135954999579*gamma_inv[1]*temp_sq[17]+0.2857142857142857*gamma_inv[8]*temp_sq[15]+0.31943828249996997*gamma_inv[5]*temp_sq[15]+0.4472135954999579*gamma_inv[4]*temp_sq[15]+0.5*gamma_inv[0]*temp_sq[15]+0.4*gamma_inv[3]*temp_sq[14]+0.31943828249996997*gamma_inv[7]*temp_sq[12]+0.5000000000000001*gamma_inv[1]*temp_sq[12]+0.4472135954999579*gamma_inv[7]*temp_sq[10]+0.4*gamma_inv[6]*temp_sq[9]+0.4472135954999579*gamma_inv[2]*temp_sq[9]+0.4472135954999579*temp_sq[4]*gamma_inv[8]+1.0000000000000002*GammaV_sq[1]*gamma[7]+0.5000000000000001*temp_sq[1]*gamma_inv[7]+0.4472135954999579*gamma_inv[3]*temp_sq[5]+0.5*temp_sq[4]*gamma_inv[5]; 
  p_fac[16] = 0.20408163265306123*gamma_inv[8]*temp_sq[16]+0.31943828249996997*gamma_inv[5]*temp_sq[16]+0.31943828249996997*gamma_inv[4]*temp_sq[16]+0.5*gamma_inv[0]*temp_sq[16]+0.2857142857142857*gamma_inv[7]*temp_sq[13]+0.44721359549995804*gamma_inv[1]*temp_sq[13]+0.2857142857142857*gamma_inv[6]*temp_sq[11]+0.44721359549995804*gamma_inv[2]*temp_sq[11]+GammaV_sq[0]*gamma[8]+0.31943828249996997*gamma_inv[8]*temp_sq[8]+0.5*gamma_inv[4]*temp_sq[8]+0.31943828249996997*temp_sq[7]*gamma_inv[8]+0.5*temp_sq[0]*gamma_inv[8]-1.4142135623730951*gamma_inv[8]+0.5*gamma_inv[5]*temp_sq[7]+0.44721359549995804*temp_sq[2]*gamma_inv[7]+0.4*gamma_inv[3]*temp_sq[6]+0.44721359549995804*temp_sq[3]*gamma_inv[6]; 
  p_fac[17] = 0.20408163265306123*gamma_inv[8]*temp_sq[17]+0.31943828249996997*gamma_inv[5]*temp_sq[17]+0.31943828249996997*gamma_inv[4]*temp_sq[17]+0.5*gamma_inv[0]*temp_sq[17]+0.28571428571428575*gamma_inv[7]*temp_sq[15]+0.4472135954999579*gamma_inv[1]*temp_sq[15]+0.28571428571428575*gamma_inv[6]*temp_sq[14]+0.4472135954999579*gamma_inv[2]*temp_sq[14]+0.31943828249997*gamma_inv[8]*temp_sq[12]+0.5000000000000001*gamma_inv[4]*temp_sq[12]+0.31943828249997*gamma_inv[8]*temp_sq[10]+0.5000000000000001*gamma_inv[5]*temp_sq[10]+0.4*gamma_inv[3]*temp_sq[9]+GammaV_sq[1]*gamma[8]+0.5*temp_sq[1]*gamma_inv[8]+0.44721359549995804*temp_sq[4]*gamma_inv[7]+0.44721359549995804*temp_sq[5]*gamma_inv[6]; 

  sr_pressure[0] += (volFact/2)*(0.7071067811865475*f[17]*p_fac[17]+0.7071067811865475*f[16]*p_fac[16]+0.7071067811865475*f[15]*p_fac[15]+0.7071067811865475*f[14]*p_fac[14]+0.7071067811865475*f[13]*p_fac[13]+0.7071067811865475*f[12]*p_fac[12]+0.7071067811865475*f[11]*p_fac[11]+0.7071067811865475*f[10]*p_fac[10]+0.7071067811865475*f[9]*p_fac[9]+0.7071067811865475*f[8]*p_fac[8]+0.7071067811865475*f[7]*p_fac[7]+0.7071067811865475*f[6]*p_fac[6]+0.7071067811865475*f[5]*p_fac[5]+0.7071067811865475*f[4]*p_fac[4]+0.7071067811865475*f[3]*p_fac[3]+0.7071067811865475*f[2]*p_fac[2]+0.7071067811865475*f[1]*p_fac[1]+0.7071067811865475*f[0]*p_fac[0]); 
  sr_pressure[1] += (volFact/2)*(0.7071067811865475*f[16]*p_fac[17]+0.7071067811865475*p_fac[16]*f[17]+0.7071067811865475*f[13]*p_fac[15]+0.7071067811865475*p_fac[13]*f[15]+0.7071067811865475*f[11]*p_fac[14]+0.7071067811865475*p_fac[11]*f[14]+0.7071067811865475*f[8]*p_fac[12]+0.7071067811865475*p_fac[8]*f[12]+0.7071067811865475*f[7]*p_fac[10]+0.7071067811865475*p_fac[7]*f[10]+0.7071067811865475*f[6]*p_fac[9]+0.7071067811865475*p_fac[6]*f[9]+0.7071067811865475*f[3]*p_fac[5]+0.7071067811865475*p_fac[3]*f[5]+0.7071067811865475*f[2]*p_fac[4]+0.7071067811865475*p_fac[2]*f[4]+0.7071067811865475*f[0]*p_fac[1]+0.7071067811865475*p_fac[0]*f[1]); 

} 
