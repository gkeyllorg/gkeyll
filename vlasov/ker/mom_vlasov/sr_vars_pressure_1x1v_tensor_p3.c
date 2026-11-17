#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p3_exp_sq.h> 
GKYL_CU_DH void sr_vars_pressure_1x1v_tensor_p3(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
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
 
  double temp[16] = {0.0}; 
  double temp_sq[16] = {0.0}; 
  double p_fac[16] = {0.0}; 
  double p0_sq[4] = {0.0}; 
  ser_1x_p3_exp_sq(p0, p0_sq); 
  temp[0] = V_0[0]*p0[0]; 
  temp[1] = p0[0]*V_0[1]; 
  temp[2] = V_0[0]*p0[1]; 
  temp[3] = V_0[1]*p0[1]; 
  temp[4] = p0[0]*V_0[2]; 
  temp[5] = V_0[0]*p0[2]; 
  temp[6] = 1.0000000000000002*p0[1]*V_0[2]; 
  temp[7] = 1.0000000000000002*V_0[1]*p0[2]; 
  temp[8] = p0[0]*V_0[3]; 
  temp[9] = V_0[0]*p0[3]; 
  temp[10] = V_0[2]*p0[2]; 
  temp[11] = 1.0*p0[1]*V_0[3]; 
  temp[12] = 1.0*V_0[1]*p0[3]; 
  temp[13] = 1.0*p0[2]*V_0[3]; 
  temp[14] = 1.0*V_0[2]*p0[3]; 
  temp[15] = V_0[3]*p0[3]; 

  temp_sq[0] = V_0_sq[0]*p0_sq[0]; 
  temp_sq[1] = p0_sq[0]*V_0_sq[1]; 
  temp_sq[2] = V_0_sq[0]*p0_sq[1]; 
  temp_sq[3] = V_0_sq[1]*p0_sq[1]; 
  temp_sq[4] = p0_sq[0]*V_0_sq[2]; 
  temp_sq[5] = V_0_sq[0]*p0_sq[2]; 
  temp_sq[6] = 1.0000000000000002*p0_sq[1]*V_0_sq[2]; 
  temp_sq[7] = 1.0000000000000002*V_0_sq[1]*p0_sq[2]; 
  temp_sq[8] = p0_sq[0]*V_0_sq[3]; 
  temp_sq[9] = V_0_sq[0]*p0_sq[3]; 
  temp_sq[10] = V_0_sq[2]*p0_sq[2]; 
  temp_sq[11] = 1.0*p0_sq[1]*V_0_sq[3]; 
  temp_sq[12] = 1.0*V_0_sq[1]*p0_sq[3]; 
  temp_sq[13] = 1.0*p0_sq[2]*V_0_sq[3]; 
  temp_sq[14] = 1.0*V_0_sq[2]*p0_sq[3]; 
  temp_sq[15] = V_0_sq[3]*p0_sq[3]; 

  p_fac[0] = 0.7071067811865475*gamma_inv[3]*temp_sq[9]-1.4142135623730951*GammaV[3]*temp[8]+0.7071067811865475*gamma_inv[2]*temp_sq[5]-1.4142135623730951*GammaV[2]*temp[4]+0.7071067811865475*gamma_inv[1]*temp_sq[2]-1.4142135623730951*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.7071067811865475*gamma_inv[0]*temp_sq[0]-1.4142135623730951*GammaV[0]*temp[0]-1.4142135623730951*gamma_inv[0]; 
  p_fac[1] = 0.7071067811865474*gamma_inv[3]*temp_sq[12]-1.2421180068162374*GammaV[2]*temp[8]+0.7071067811865475*gamma_inv[2]*temp_sq[7]-1.2421180068162374*GammaV[3]*temp[4]-1.264911064067352*GammaV[1]*temp[4]+0.7071067811865475*gamma_inv[1]*temp_sq[3]-1.264911064067352*temp[1]*GammaV[2]+0.7071067811865475*gamma_inv[0]*temp_sq[1]-1.4142135623730951*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.4142135623730951*temp[0]*GammaV[1]; 
  p_fac[2] = -(1.414213562373095*GammaV[3]*temp[11])+0.6210590034081186*gamma_inv[2]*temp_sq[9]-1.4142135623730951*GammaV[2]*temp[6]+0.6210590034081186*gamma_inv[3]*temp_sq[5]+0.6324555320336759*gamma_inv[1]*temp_sq[5]-1.4142135623730951*GammaV[1]*temp[3]+0.6324555320336759*gamma_inv[2]*temp_sq[2]+0.7071067811865475*gamma_inv[0]*temp_sq[2]-1.4142135623730951*GammaV[0]*temp[2]+GammaV_sq[0]*gamma[1]+0.7071067811865475*temp_sq[0]*gamma_inv[1]-1.4142135623730951*gamma_inv[1]; 
  p_fac[3] = 0.6210590034081187*gamma_inv[2]*temp_sq[12]-1.2421180068162376*GammaV[2]*temp[11]+0.6210590034081187*gamma_inv[3]*temp_sq[7]+0.632455532033676*gamma_inv[1]*temp_sq[7]-1.2421180068162376*GammaV[3]*temp[6]-1.264911064067352*GammaV[1]*temp[6]+0.6324555320336759*gamma_inv[2]*temp_sq[3]+0.7071067811865475*gamma_inv[0]*temp_sq[3]-1.264911064067352*GammaV[2]*temp[3]-1.4142135623730951*GammaV[0]*temp[3]-1.4142135623730951*GammaV[1]*temp[2]+GammaV_sq[1]*gamma[1]+0.7071067811865475*gamma_inv[1]*temp_sq[1]; 
  p_fac[4] = 0.7071067811865474*gamma_inv[3]*temp_sq[14]+0.7071067811865475*gamma_inv[2]*temp_sq[10]-0.8432740427115681*GammaV[3]*temp[8]-1.2421180068162374*GammaV[1]*temp[8]+0.7071067811865475*gamma_inv[1]*temp_sq[6]+0.7071067811865475*gamma_inv[0]*temp_sq[4]-0.9035079029052515*GammaV[2]*temp[4]-1.4142135623730951*GammaV[0]*temp[4]-1.2421180068162374*temp[1]*GammaV[3]+gamma[0]*GammaV_sq[2]-1.4142135623730951*temp[0]*GammaV[2]-1.264911064067352*GammaV[1]*temp[1]; 
  p_fac[5] = -(1.4142135623730951*GammaV[3]*temp[13])-1.4142135623730951*GammaV[2]*temp[10]+0.421637021355784*gamma_inv[3]*temp_sq[9]+0.6210590034081186*gamma_inv[1]*temp_sq[9]-1.4142135623730951*GammaV[1]*temp[7]+0.45175395145262565*gamma_inv[2]*temp_sq[5]+0.7071067811865475*gamma_inv[0]*temp_sq[5]-1.4142135623730951*GammaV[0]*temp[5]+0.6210590034081186*temp_sq[2]*gamma_inv[3]+GammaV_sq[0]*gamma[2]+0.6324555320336759*gamma_inv[1]*temp_sq[2]+0.7071067811865475*temp_sq[0]*gamma_inv[2]-1.4142135623730951*gamma_inv[2]; 
  p_fac[6] = 0.6210590034081186*gamma_inv[2]*temp_sq[14]-0.8432740427115679*GammaV[3]*temp[11]-1.2421180068162374*GammaV[1]*temp[11]+0.6210590034081188*gamma_inv[3]*temp_sq[10]+0.632455532033676*gamma_inv[1]*temp_sq[10]+0.6324555320336759*gamma_inv[2]*temp_sq[6]+0.7071067811865475*gamma_inv[0]*temp_sq[6]-0.9035079029052515*GammaV[2]*temp[6]-1.4142135623730951*GammaV[0]*temp[6]+0.7071067811865475*gamma_inv[1]*temp_sq[4]-1.2421180068162379*GammaV[3]*temp[3]-1.264911064067352*GammaV[1]*temp[3]-1.4142135623730951*GammaV[2]*temp[2]+1.0000000000000002*gamma[1]*GammaV_sq[2]; 
  p_fac[7] = -(1.2421180068162376*GammaV[2]*temp[13])+0.4216370213557839*gamma_inv[3]*temp_sq[12]+0.6210590034081187*gamma_inv[1]*temp_sq[12]-1.2421180068162379*GammaV[3]*temp[10]-1.264911064067352*GammaV[1]*temp[10]+0.45175395145262565*gamma_inv[2]*temp_sq[7]+0.7071067811865475*gamma_inv[0]*temp_sq[7]-1.264911064067352*GammaV[2]*temp[7]-1.4142135623730951*GammaV[0]*temp[7]-1.4142135623730951*GammaV[1]*temp[5]+0.6210590034081188*gamma_inv[3]*temp_sq[3]+0.632455532033676*gamma_inv[1]*temp_sq[3]+1.0000000000000002*GammaV_sq[1]*gamma[2]+0.7071067811865475*temp_sq[1]*gamma_inv[2]; 
  p_fac[8] = 0.7071067811865475*gamma_inv[3]*temp_sq[15]+0.7071067811865474*gamma_inv[2]*temp_sq[13]+0.7071067811865474*gamma_inv[1]*temp_sq[11]+0.7071067811865475*gamma_inv[0]*temp_sq[8]-0.8432740427115681*GammaV[2]*temp[8]-1.4142135623730951*GammaV[0]*temp[8]-0.8432740427115681*GammaV[3]*temp[4]-1.2421180068162374*GammaV[1]*temp[4]+gamma[0]*GammaV_sq[3]-1.4142135623730951*temp[0]*GammaV[3]-1.2421180068162374*temp[1]*GammaV[2]; 
  p_fac[9] = -(1.4142135623730951*GammaV[3]*temp[15])-1.4142135623730951*GammaV[2]*temp[14]-1.414213562373095*GammaV[1]*temp[12]+0.421637021355784*gamma_inv[2]*temp_sq[9]+0.7071067811865475*gamma_inv[0]*temp_sq[9]-1.4142135623730951*GammaV[0]*temp[9]+0.421637021355784*gamma_inv[3]*temp_sq[5]+0.6210590034081186*gamma_inv[1]*temp_sq[5]+GammaV_sq[0]*gamma[3]+0.7071067811865475*temp_sq[0]*gamma_inv[3]-1.4142135623730951*gamma_inv[3]+0.6210590034081186*gamma_inv[2]*temp_sq[2]; 
  p_fac[10] = 0.42163702135578396*gamma_inv[3]*temp_sq[14]+0.6210590034081186*gamma_inv[1]*temp_sq[14]-0.8432740427115681*GammaV[3]*temp[13]-1.2421180068162374*GammaV[1]*temp[13]+0.45175395145262565*gamma_inv[2]*temp_sq[10]+0.7071067811865475*gamma_inv[0]*temp_sq[10]-0.9035079029052515*GammaV[2]*temp[10]-1.4142135623730951*GammaV[0]*temp[10]-1.2421180068162376*GammaV[3]*temp[7]-1.264911064067352*GammaV[1]*temp[7]+0.6210590034081187*gamma_inv[3]*temp_sq[6]+0.632455532033676*gamma_inv[1]*temp_sq[6]-1.4142135623730951*GammaV[2]*temp[5]+0.7071067811865475*gamma_inv[2]*temp_sq[4]+GammaV_sq[2]*gamma[2]; 
  p_fac[11] = 0.6210590034081186*gamma_inv[2]*temp_sq[15]+0.6210590034081186*gamma_inv[3]*temp_sq[13]+0.6324555320336759*gamma_inv[1]*temp_sq[13]+0.6324555320336759*gamma_inv[2]*temp_sq[11]+0.7071067811865475*gamma_inv[0]*temp_sq[11]-0.8432740427115681*GammaV[2]*temp[11]-1.4142135623730951*GammaV[0]*temp[11]+0.7071067811865474*gamma_inv[1]*temp_sq[8]-0.8432740427115679*GammaV[3]*temp[6]-1.2421180068162374*GammaV[1]*temp[6]-1.2421180068162374*GammaV[2]*temp[3]+1.0*gamma[1]*GammaV_sq[3]-1.414213562373095*temp[2]*GammaV[3]; 
  p_fac[12] = -(1.2421180068162374*GammaV[2]*temp[15])-1.2421180068162372*GammaV[3]*temp[14]-1.264911064067352*GammaV[1]*temp[14]+0.421637021355784*gamma_inv[2]*temp_sq[12]+0.7071067811865475*gamma_inv[0]*temp_sq[12]-1.264911064067352*GammaV[2]*temp[12]-1.4142135623730951*GammaV[0]*temp[12]-1.414213562373095*GammaV[1]*temp[9]+0.4216370213557839*gamma_inv[3]*temp_sq[7]+0.6210590034081187*gamma_inv[1]*temp_sq[7]+1.0*GammaV_sq[1]*gamma[3]+0.6210590034081186*gamma_inv[2]*temp_sq[3]+0.7071067811865474*temp_sq[1]*gamma_inv[3]; 
  p_fac[13] = 0.42163702135578396*gamma_inv[3]*temp_sq[15]+0.6210590034081186*gamma_inv[1]*temp_sq[15]+0.45175395145262565*gamma_inv[2]*temp_sq[13]+0.7071067811865475*gamma_inv[0]*temp_sq[13]-0.8432740427115681*GammaV[2]*temp[13]-1.4142135623730951*GammaV[0]*temp[13]+0.6210590034081186*gamma_inv[3]*temp_sq[11]+0.6324555320336759*gamma_inv[1]*temp_sq[11]-0.8432740427115681*GammaV[3]*temp[10]-1.2421180068162374*GammaV[1]*temp[10]+0.7071067811865474*gamma_inv[2]*temp_sq[8]-1.2421180068162376*GammaV[2]*temp[7]-1.4142135623730951*GammaV[3]*temp[5]+1.0*gamma[2]*GammaV_sq[3]; 
  p_fac[14] = -(0.8432740427115681*GammaV[3]*temp[15])-1.2421180068162374*GammaV[1]*temp[15]+0.421637021355784*gamma_inv[2]*temp_sq[14]+0.7071067811865475*gamma_inv[0]*temp_sq[14]-0.9035079029052515*GammaV[2]*temp[14]-1.4142135623730951*GammaV[0]*temp[14]-1.2421180068162372*GammaV[3]*temp[12]-1.264911064067352*GammaV[1]*temp[12]+0.42163702135578396*gamma_inv[3]*temp_sq[10]+0.6210590034081186*gamma_inv[1]*temp_sq[10]-1.4142135623730951*GammaV[2]*temp[9]+0.6210590034081186*gamma_inv[2]*temp_sq[6]+0.7071067811865474*gamma_inv[3]*temp_sq[4]+1.0*GammaV_sq[2]*gamma[3]; 
  p_fac[15] = 0.421637021355784*gamma_inv[2]*temp_sq[15]+0.7071067811865475*gamma_inv[0]*temp_sq[15]-0.8432740427115681*GammaV[2]*temp[15]-1.4142135623730951*GammaV[0]*temp[15]-0.8432740427115681*GammaV[3]*temp[14]-1.2421180068162374*GammaV[1]*temp[14]+0.42163702135578396*gamma_inv[3]*temp_sq[13]+0.6210590034081186*gamma_inv[1]*temp_sq[13]-1.2421180068162376*GammaV[2]*temp[12]+0.6210590034081187*gamma_inv[2]*temp_sq[11]-1.4142135623730951*GammaV[3]*temp[9]+0.7071067811865475*gamma_inv[3]*temp_sq[8]+GammaV_sq[3]*gamma[3]; 

  sr_pressure[0] += (volFact/1)*(0.7071067811865475*f[15]*p_fac[15]+0.7071067811865475*f[14]*p_fac[14]+0.7071067811865475*f[13]*p_fac[13]+0.7071067811865475*f[12]*p_fac[12]+0.7071067811865475*f[11]*p_fac[11]+0.7071067811865475*f[10]*p_fac[10]+0.7071067811865475*f[9]*p_fac[9]+0.7071067811865475*f[8]*p_fac[8]+0.7071067811865475*f[7]*p_fac[7]+0.7071067811865475*f[6]*p_fac[6]+0.7071067811865475*f[5]*p_fac[5]+0.7071067811865475*f[4]*p_fac[4]+0.7071067811865475*f[3]*p_fac[3]+0.7071067811865475*f[2]*p_fac[2]+0.7071067811865475*f[1]*p_fac[1]+0.7071067811865475*f[0]*p_fac[0]); 
  sr_pressure[1] += (volFact/1)*(0.6210590034081186*f[14]*p_fac[15]+0.6210590034081186*p_fac[14]*f[15]+0.6324555320336761*f[12]*p_fac[14]+0.6324555320336761*p_fac[12]*f[14]+0.6210590034081186*f[10]*p_fac[13]+0.6210590034081186*p_fac[10]*f[13]+0.7071067811865474*f[9]*p_fac[12]+0.7071067811865474*p_fac[9]*f[12]+0.6210590034081187*f[6]*p_fac[11]+0.6210590034081187*p_fac[6]*f[11]+0.6324555320336759*f[7]*p_fac[10]+0.6324555320336759*p_fac[7]*f[10]+0.6210590034081186*f[4]*p_fac[8]+0.6210590034081186*p_fac[4]*f[8]+0.7071067811865475*f[5]*p_fac[7]+0.7071067811865475*p_fac[5]*f[7]+0.6324555320336759*f[3]*p_fac[6]+0.6324555320336759*p_fac[3]*f[6]+0.6324555320336759*f[1]*p_fac[4]+0.6324555320336759*p_fac[1]*f[4]+0.7071067811865475*f[2]*p_fac[3]+0.7071067811865475*p_fac[2]*f[3]+0.7071067811865475*f[0]*p_fac[1]+0.7071067811865475*p_fac[0]*f[1]); 
  sr_pressure[2] += (volFact/1)*(0.421637021355784*f[15]*p_fac[15]+0.6210590034081187*f[12]*p_fac[15]+0.6210590034081187*p_fac[12]*f[15]+0.45175395145262565*f[14]*p_fac[14]+0.7071067811865474*f[9]*p_fac[14]+0.7071067811865474*p_fac[9]*f[14]+0.421637021355784*f[13]*p_fac[13]+0.6210590034081187*f[7]*p_fac[13]+0.6210590034081187*p_fac[7]*f[13]+0.6324555320336759*f[12]*p_fac[12]+0.421637021355784*f[11]*p_fac[11]+0.6210590034081187*f[3]*p_fac[11]+0.6210590034081187*p_fac[3]*f[11]+0.45175395145262565*f[10]*p_fac[10]+0.7071067811865475*f[5]*p_fac[10]+0.7071067811865475*p_fac[5]*f[10]+0.421637021355784*f[8]*p_fac[8]+0.6210590034081186*f[1]*p_fac[8]+0.6210590034081186*p_fac[1]*f[8]+0.6324555320336759*f[7]*p_fac[7]+0.45175395145262565*f[6]*p_fac[6]+0.7071067811865475*f[2]*p_fac[6]+0.7071067811865475*p_fac[2]*f[6]+0.45175395145262565*f[4]*p_fac[4]+0.7071067811865475*f[0]*p_fac[4]+0.7071067811865475*p_fac[0]*f[4]+0.6324555320336759*f[3]*p_fac[3]+0.6324555320336759*f[1]*p_fac[1]); 
  sr_pressure[3] += (volFact/1)*(0.421637021355784*f[14]*p_fac[15]+0.7071067811865475*f[9]*p_fac[15]+0.421637021355784*p_fac[14]*f[15]+0.7071067811865475*p_fac[9]*f[15]+0.6210590034081188*f[12]*p_fac[14]+0.6210590034081188*p_fac[12]*f[14]+0.421637021355784*f[10]*p_fac[13]+0.7071067811865474*f[5]*p_fac[13]+0.421637021355784*p_fac[10]*f[13]+0.7071067811865474*p_fac[5]*f[13]+0.42163702135578385*f[6]*p_fac[11]+0.7071067811865474*f[2]*p_fac[11]+0.42163702135578385*p_fac[6]*f[11]+0.7071067811865474*p_fac[2]*f[11]+0.6210590034081187*f[7]*p_fac[10]+0.6210590034081187*p_fac[7]*f[10]+0.421637021355784*f[4]*p_fac[8]+0.7071067811865475*f[0]*p_fac[8]+0.421637021355784*p_fac[4]*f[8]+0.7071067811865475*p_fac[0]*f[8]+0.6210590034081187*f[3]*p_fac[6]+0.6210590034081187*p_fac[3]*f[6]+0.6210590034081186*f[1]*p_fac[4]+0.6210590034081186*p_fac[1]*f[4]); 

} 
