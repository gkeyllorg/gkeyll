#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_1x_p3_exp_sq.h> 
GKYL_CU_DH void sr_vars_pressure_2x1v_ser_p2(const double *w, const double *dxv, const double *vmap, const double *gamma, const double *gamma_inv, const double *u_i, const double *u_i_sq, const double *GammaV, const double *GammaV_sq, const double *f, double* GKYL_RESTRICT sr_pressure) 
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
  const double volFact = dxv[2]/2; 
 
  const double *p0 = &vmap[0]; 
  const double *V_0 = &u_i[0]; 
  const double *V_0_sq = &u_i_sq[0]; 
 
  double temp[20] = {0.0}; 
  double temp_sq[20] = {0.0}; 
  double p_fac[20] = {0.0}; 
  double p0_sq[4] = {0.0}; 
  ser_1x_p3_exp_sq(p0, p0_sq); 
  temp[0] = V_0[0]*p0[0]; 
  temp[1] = p0[0]*V_0[1]; 
  temp[2] = p0[0]*V_0[2]; 
  temp[3] = V_0[0]*p0[1]; 
  temp[4] = p0[0]*V_0[3]; 
  temp[5] = V_0[1]*p0[1]; 
  temp[6] = p0[1]*V_0[2]; 
  temp[7] = p0[0]*V_0[4]; 
  temp[8] = p0[0]*V_0[5]; 
  temp[10] = p0[1]*V_0[3]; 
  temp[11] = p0[0]*V_0[6]; 
  temp[12] = p0[0]*V_0[7]; 
  temp[13] = 1.0000000000000002*p0[1]*V_0[4]; 
  temp[14] = 1.0000000000000002*p0[1]*V_0[5]; 
  temp[17] = 1.0000000000000002*p0[1]*V_0[6]; 
  temp[18] = 1.0000000000000002*p0[1]*V_0[7]; 

  temp_sq[0] = V_0_sq[0]*p0_sq[0]; 
  temp_sq[1] = p0_sq[0]*V_0_sq[1]; 
  temp_sq[2] = p0_sq[0]*V_0_sq[2]; 
  temp_sq[3] = V_0_sq[0]*p0_sq[1]; 
  temp_sq[4] = p0_sq[0]*V_0_sq[3]; 
  temp_sq[5] = V_0_sq[1]*p0_sq[1]; 
  temp_sq[6] = p0_sq[1]*V_0_sq[2]; 
  temp_sq[7] = p0_sq[0]*V_0_sq[4]; 
  temp_sq[8] = p0_sq[0]*V_0_sq[5]; 
  temp_sq[9] = V_0_sq[0]*p0_sq[2]; 
  temp_sq[10] = p0_sq[1]*V_0_sq[3]; 
  temp_sq[11] = p0_sq[0]*V_0_sq[6]; 
  temp_sq[12] = p0_sq[0]*V_0_sq[7]; 
  temp_sq[13] = 1.0000000000000002*p0_sq[1]*V_0_sq[4]; 
  temp_sq[14] = 1.0000000000000002*p0_sq[1]*V_0_sq[5]; 
  temp_sq[15] = 1.0000000000000002*V_0_sq[1]*p0_sq[2]; 
  temp_sq[16] = 1.0000000000000002*V_0_sq[2]*p0_sq[2]; 
  temp_sq[17] = 1.0000000000000002*p0_sq[1]*V_0_sq[6]; 
  temp_sq[18] = 1.0000000000000002*p0_sq[1]*V_0_sq[7]; 
  temp_sq[19] = p0_sq[2]*V_0_sq[3]; 

  p_fac[0] = -(1.0*GammaV[7]*temp[12])-1.0*GammaV[6]*temp[11]+0.7071067811865475*gamma_inv[2]*temp_sq[9]-1.0*GammaV[5]*temp[8]-1.0*GammaV[4]*temp[7]-1.0*GammaV[3]*temp[4]+0.7071067811865475*gamma_inv[1]*temp_sq[3]-1.0*GammaV[2]*temp[2]-1.0*GammaV[1]*temp[1]+GammaV_sq[0]*gamma[0]+0.7071067811865475*gamma_inv[0]*temp_sq[0]-1.0*GammaV[0]*temp[0]-2.0*gamma_inv[0]; 
  p_fac[1] = 0.7071067811865475*gamma_inv[2]*temp_sq[15]-1.0000000000000002*GammaV[5]*temp[12]-0.8944271909999161*GammaV[3]*temp[11]-1.0000000000000002*GammaV[7]*temp[8]-0.8944271909999159*GammaV[1]*temp[7]-0.8944271909999161*temp[4]*GammaV[6]+0.7071067811865475*gamma_inv[1]*temp_sq[5]-1.0*GammaV[2]*temp[4]-0.8944271909999159*temp[1]*GammaV[4]-1.0*temp[2]*GammaV[3]+0.7071067811865475*gamma_inv[0]*temp_sq[1]-1.0*GammaV[0]*temp[1]+gamma[0]*GammaV_sq[1]-1.0*temp[0]*GammaV[1]; 
  p_fac[2] = 0.7071067811865475*gamma_inv[2]*temp_sq[16]-0.8944271909999161*GammaV[3]*temp[12]-1.0000000000000002*GammaV[4]*temp[11]-0.8944271909999159*GammaV[2]*temp[8]-1.0000000000000002*GammaV[6]*temp[7]-0.8944271909999161*temp[4]*GammaV[7]+0.7071067811865475*gamma_inv[1]*temp_sq[6]-0.8944271909999159*temp[2]*GammaV[5]-1.0*GammaV[1]*temp[4]-1.0*temp[1]*GammaV[3]+0.7071067811865475*gamma_inv[0]*temp_sq[2]-1.0*GammaV[0]*temp[2]+gamma[0]*GammaV_sq[2]-1.0*temp[0]*GammaV[2]; 
  p_fac[3] = -(1.0000000000000002*GammaV[7]*temp[18])-1.0000000000000002*GammaV[6]*temp[17]-1.0000000000000002*GammaV[5]*temp[14]-1.0000000000000002*GammaV[4]*temp[13]-1.0*GammaV[3]*temp[10]+0.6324555320336759*gamma_inv[1]*temp_sq[9]-1.0*GammaV[2]*temp[6]-1.0*GammaV[1]*temp[5]+0.6324555320336759*gamma_inv[2]*temp_sq[3]+0.7071067811865475*gamma_inv[0]*temp_sq[3]-1.0*GammaV[0]*temp[3]+GammaV_sq[0]*gamma[1]+0.7071067811865475*temp_sq[0]*gamma_inv[1]-2.0*gamma_inv[1]; 
  p_fac[4] = 0.7071067811865475*gamma_inv[2]*temp_sq[19]-0.8*GammaV[6]*temp[12]-0.8944271909999161*GammaV[2]*temp[12]-0.8*GammaV[7]*temp[11]-0.8944271909999161*GammaV[1]*temp[11]+0.7071067811865475*gamma_inv[1]*temp_sq[10]-0.8944271909999159*GammaV[3]*temp[8]-0.8944271909999159*GammaV[3]*temp[7]-0.8944271909999161*temp[2]*GammaV[7]-0.8944271909999161*temp[1]*GammaV[6]-0.8944271909999159*temp[4]*GammaV[5]+0.7071067811865475*gamma_inv[0]*temp_sq[4]-0.8944271909999159*GammaV[4]*temp[4]-1.0*GammaV[0]*temp[4]+gamma[0]*GammaV_sq[3]-1.0*temp[0]*GammaV[3]-1.0*GammaV[1]*temp[2]-1.0*temp[1]*GammaV[2]; 
  p_fac[5] = -(1.0*GammaV[5]*temp[18])-0.8944271909999159*GammaV[3]*temp[17]+0.632455532033676*gamma_inv[1]*temp_sq[15]-1.0*GammaV[7]*temp[14]-0.8944271909999161*GammaV[1]*temp[13]-0.8944271909999161*GammaV[6]*temp[10]-1.0*GammaV[2]*temp[10]-1.0*GammaV[3]*temp[6]+0.6324555320336759*gamma_inv[2]*temp_sq[5]+0.7071067811865475*gamma_inv[0]*temp_sq[5]-0.8944271909999159*GammaV[4]*temp[5]-1.0*GammaV[0]*temp[5]-1.0*GammaV[1]*temp[3]+GammaV_sq[1]*gamma[1]+0.7071067811865475*gamma_inv[1]*temp_sq[1]; 
  p_fac[6] = -(0.8944271909999159*GammaV[3]*temp[18])-1.0*GammaV[4]*temp[17]+0.632455532033676*gamma_inv[1]*temp_sq[16]-0.8944271909999161*GammaV[2]*temp[14]-1.0*GammaV[6]*temp[13]-0.8944271909999161*GammaV[7]*temp[10]-1.0*GammaV[1]*temp[10]+0.6324555320336759*gamma_inv[2]*temp_sq[6]+0.7071067811865475*gamma_inv[0]*temp_sq[6]-0.8944271909999159*GammaV[5]*temp[6]-1.0*GammaV[0]*temp[6]-1.0*GammaV[3]*temp[5]-1.0*GammaV[2]*temp[3]+0.7071067811865475*gamma_inv[1]*temp_sq[2]+gamma[1]*GammaV_sq[2]; 
  p_fac[7] = 0.7071067811865475*gamma_inv[1]*temp_sq[13]-0.8944271909999159*GammaV[7]*temp[12]-0.6388765649999399*GammaV[6]*temp[11]-1.0000000000000002*GammaV[2]*temp[11]+0.7071067811865475*gamma_inv[0]*temp_sq[7]-0.6388765649999399*GammaV[4]*temp[7]-1.0*GammaV[0]*temp[7]-1.0000000000000002*temp[2]*GammaV[6]-0.8944271909999159*GammaV[3]*temp[4]+gamma[0]*GammaV_sq[4]-1.0*temp[0]*GammaV[4]-0.8944271909999159*GammaV[1]*temp[1]; 
  p_fac[8] = 0.7071067811865475*gamma_inv[1]*temp_sq[14]-0.6388765649999399*GammaV[7]*temp[12]-1.0000000000000002*GammaV[1]*temp[12]-0.8944271909999159*GammaV[6]*temp[11]+0.7071067811865475*gamma_inv[0]*temp_sq[8]-0.6388765649999399*GammaV[5]*temp[8]-1.0*GammaV[0]*temp[8]-1.0000000000000002*temp[1]*GammaV[7]+gamma[0]*GammaV_sq[5]-1.0*temp[0]*GammaV[5]-0.8944271909999159*GammaV[3]*temp[4]-0.8944271909999159*GammaV[2]*temp[2]; 
  p_fac[9] = 0.45175395145262565*gamma_inv[2]*temp_sq[9]+0.7071067811865475*gamma_inv[0]*temp_sq[9]+0.6324555320336759*gamma_inv[1]*temp_sq[3]+GammaV_sq[0]*gamma[2]+0.7071067811865475*temp_sq[0]*gamma_inv[2]-2.0*gamma_inv[2]; 
  p_fac[10] = 0.6324555320336759*gamma_inv[1]*temp_sq[19]-0.8*GammaV[6]*temp[18]-0.8944271909999159*GammaV[2]*temp[18]-0.8*GammaV[7]*temp[17]-0.8944271909999159*GammaV[1]*temp[17]-0.8944271909999161*GammaV[3]*temp[14]-0.8944271909999161*GammaV[3]*temp[13]+0.6324555320336759*gamma_inv[2]*temp_sq[10]+0.7071067811865475*gamma_inv[0]*temp_sq[10]-0.8944271909999159*GammaV[5]*temp[10]-0.8944271909999159*GammaV[4]*temp[10]-1.0*GammaV[0]*temp[10]-0.8944271909999161*temp[6]*GammaV[7]-1.0*GammaV[1]*temp[6]-0.8944271909999161*temp[5]*GammaV[6]-1.0*GammaV[2]*temp[5]+0.7071067811865475*gamma_inv[1]*temp_sq[4]-1.0*GammaV[3]*temp[3]+gamma[1]*GammaV_sq[3]; 
  p_fac[11] = 0.7071067811865475*gamma_inv[1]*temp_sq[17]-0.8*GammaV[3]*temp[12]+0.7071067811865475*gamma_inv[0]*temp_sq[11]-0.8944271909999159*GammaV[5]*temp[11]-0.6388765649999399*GammaV[4]*temp[11]-1.0*GammaV[0]*temp[11]-0.8944271909999159*GammaV[6]*temp[8]-0.6388765649999399*GammaV[6]*temp[7]-1.0000000000000002*GammaV[2]*temp[7]-0.8*temp[4]*GammaV[7]+gamma[0]*GammaV_sq[6]-1.0*temp[0]*GammaV[6]-0.8944271909999161*GammaV[1]*temp[4]-1.0000000000000002*temp[2]*GammaV[4]-0.8944271909999161*temp[1]*GammaV[3]; 
  p_fac[12] = 0.7071067811865475*gamma_inv[1]*temp_sq[18]+0.7071067811865475*gamma_inv[0]*temp_sq[12]-0.6388765649999399*GammaV[5]*temp[12]-0.8944271909999159*GammaV[4]*temp[12]-1.0*GammaV[0]*temp[12]-0.8*GammaV[3]*temp[11]-0.6388765649999399*GammaV[7]*temp[8]-1.0000000000000002*GammaV[1]*temp[8]-0.8944271909999159*GammaV[7]*temp[7]+gamma[0]*GammaV_sq[7]-1.0*temp[0]*GammaV[7]-0.8*temp[4]*GammaV[6]-1.0000000000000002*temp[1]*GammaV[5]-0.8944271909999161*GammaV[2]*temp[4]-0.8944271909999161*temp[2]*GammaV[3]; 
  p_fac[13] = -(0.8944271909999159*GammaV[7]*temp[18])-0.6388765649999399*GammaV[6]*temp[17]-1.0000000000000002*GammaV[2]*temp[17]+0.6324555320336759*gamma_inv[2]*temp_sq[13]+0.7071067811865475*gamma_inv[0]*temp_sq[13]-0.6388765649999399*GammaV[4]*temp[13]-1.0*GammaV[0]*temp[13]-0.8944271909999161*GammaV[3]*temp[10]+0.7071067811865475*gamma_inv[1]*temp_sq[7]-1.0*GammaV[6]*temp[6]-0.8944271909999161*GammaV[1]*temp[5]+1.0000000000000002*gamma[1]*GammaV_sq[4]-1.0000000000000002*temp[3]*GammaV[4]; 
  p_fac[14] = -(0.6388765649999399*GammaV[7]*temp[18])-1.0000000000000002*GammaV[1]*temp[18]-0.8944271909999159*GammaV[6]*temp[17]+0.6324555320336759*gamma_inv[2]*temp_sq[14]+0.7071067811865475*gamma_inv[0]*temp_sq[14]-0.6388765649999399*GammaV[5]*temp[14]-1.0*GammaV[0]*temp[14]-0.8944271909999161*GammaV[3]*temp[10]+0.7071067811865475*gamma_inv[1]*temp_sq[8]-1.0*temp[5]*GammaV[7]-0.8944271909999161*GammaV[2]*temp[6]+1.0000000000000002*gamma[1]*GammaV_sq[5]-1.0000000000000002*temp[3]*GammaV[5]; 
  p_fac[15] = 0.45175395145262565*gamma_inv[2]*temp_sq[15]+0.7071067811865475*gamma_inv[0]*temp_sq[15]+0.632455532033676*gamma_inv[1]*temp_sq[5]+1.0000000000000002*GammaV_sq[1]*gamma[2]+0.7071067811865475*temp_sq[1]*gamma_inv[2]; 
  p_fac[16] = 0.45175395145262565*gamma_inv[2]*temp_sq[16]+0.7071067811865475*gamma_inv[0]*temp_sq[16]+0.632455532033676*gamma_inv[1]*temp_sq[6]+1.0000000000000002*GammaV_sq[2]*gamma[2]+0.7071067811865475*gamma_inv[2]*temp_sq[2]; 
  p_fac[17] = -(0.8*GammaV[3]*temp[18])+0.6324555320336759*gamma_inv[2]*temp_sq[17]+0.7071067811865475*gamma_inv[0]*temp_sq[17]-0.8944271909999159*GammaV[5]*temp[17]-0.6388765649999399*GammaV[4]*temp[17]-1.0*GammaV[0]*temp[17]-0.8944271909999159*GammaV[6]*temp[14]-0.6388765649999399*GammaV[6]*temp[13]-1.0000000000000002*GammaV[2]*temp[13]+0.7071067811865475*gamma_inv[1]*temp_sq[11]-0.8*GammaV[7]*temp[10]-0.8944271909999159*GammaV[1]*temp[10]-1.0*GammaV[4]*temp[6]+1.0000000000000002*gamma[1]*GammaV_sq[6]-1.0000000000000002*temp[3]*GammaV[6]-0.8944271909999159*GammaV[3]*temp[5]; 
  p_fac[18] = 0.6324555320336759*gamma_inv[2]*temp_sq[18]+0.7071067811865475*gamma_inv[0]*temp_sq[18]-0.6388765649999399*GammaV[5]*temp[18]-0.8944271909999159*GammaV[4]*temp[18]-1.0*GammaV[0]*temp[18]-0.8*GammaV[3]*temp[17]-0.6388765649999399*GammaV[7]*temp[14]-1.0000000000000002*GammaV[1]*temp[14]-0.8944271909999159*GammaV[7]*temp[13]+0.7071067811865475*gamma_inv[1]*temp_sq[12]-0.8*GammaV[6]*temp[10]-0.8944271909999159*GammaV[2]*temp[10]+1.0000000000000002*gamma[1]*GammaV_sq[7]-1.0000000000000002*temp[3]*GammaV[7]-0.8944271909999159*GammaV[3]*temp[6]-1.0*GammaV[5]*temp[5]; 
  p_fac[19] = 0.45175395145262565*gamma_inv[2]*temp_sq[19]+0.7071067811865475*gamma_inv[0]*temp_sq[19]+0.6324555320336759*gamma_inv[1]*temp_sq[10]+0.7071067811865475*gamma_inv[2]*temp_sq[4]+gamma[2]*GammaV_sq[3]; 

  sr_pressure[0] += (volFact/1)*(0.5*f[19]*p_fac[19]+0.5*f[18]*p_fac[18]+0.5*f[17]*p_fac[17]+0.5*f[16]*p_fac[16]+0.5*f[15]*p_fac[15]+0.5*f[14]*p_fac[14]+0.5*f[13]*p_fac[13]+0.5*f[12]*p_fac[12]+0.5*f[11]*p_fac[11]+0.5*f[10]*p_fac[10]+0.5*f[9]*p_fac[9]+0.5*f[8]*p_fac[8]+0.5*f[7]*p_fac[7]+0.5*f[6]*p_fac[6]+0.5*f[5]*p_fac[5]+0.5*f[4]*p_fac[4]+0.5*f[3]*p_fac[3]+0.5*f[2]*p_fac[2]+0.5*f[1]*p_fac[1]+0.5*f[0]*p_fac[0]); 
  sr_pressure[1] += (volFact/1)*(0.5000000000000001*f[16]*p_fac[19]+0.5000000000000001*p_fac[16]*f[19]+0.5000000000000001*f[14]*p_fac[18]+0.5000000000000001*p_fac[14]*f[18]+0.4472135954999579*f[10]*p_fac[17]+0.4472135954999579*p_fac[10]*f[17]+0.5000000000000001*f[9]*p_fac[15]+0.5000000000000001*p_fac[9]*f[15]+0.44721359549995787*f[5]*p_fac[13]+0.44721359549995787*p_fac[5]*f[13]+0.5000000000000001*f[8]*p_fac[12]+0.5000000000000001*p_fac[8]*f[12]+0.44721359549995787*f[4]*p_fac[11]+0.44721359549995787*p_fac[4]*f[11]+0.5*f[6]*p_fac[10]+0.5*p_fac[6]*f[10]+0.4472135954999579*f[1]*p_fac[7]+0.4472135954999579*p_fac[1]*f[7]+0.5*f[3]*p_fac[5]+0.5*p_fac[3]*f[5]+0.5*f[2]*p_fac[4]+0.5*p_fac[2]*f[4]+0.5*f[0]*p_fac[1]+0.5*p_fac[0]*f[1]); 
  sr_pressure[2] += (volFact/1)*(0.5000000000000001*f[15]*p_fac[19]+0.5000000000000001*p_fac[15]*f[19]+0.4472135954999579*f[10]*p_fac[18]+0.4472135954999579*p_fac[10]*f[18]+0.5000000000000001*f[13]*p_fac[17]+0.5000000000000001*p_fac[13]*f[17]+0.5000000000000001*f[9]*p_fac[16]+0.5000000000000001*p_fac[9]*f[16]+0.44721359549995787*f[6]*p_fac[14]+0.44721359549995787*p_fac[6]*f[14]+0.44721359549995787*f[4]*p_fac[12]+0.44721359549995787*p_fac[4]*f[12]+0.5000000000000001*f[7]*p_fac[11]+0.5000000000000001*p_fac[7]*f[11]+0.5*f[5]*p_fac[10]+0.5*p_fac[5]*f[10]+0.4472135954999579*f[2]*p_fac[8]+0.4472135954999579*p_fac[2]*f[8]+0.5*f[3]*p_fac[6]+0.5*p_fac[3]*f[6]+0.5*f[1]*p_fac[4]+0.5*p_fac[1]*f[4]+0.5*f[0]*p_fac[2]+0.5*p_fac[0]*f[2]); 
  sr_pressure[3] += (volFact/1)*(0.5*f[9]*p_fac[19]+0.5*p_fac[9]*f[19]+0.4*f[17]*p_fac[18]+0.4472135954999579*f[6]*p_fac[18]+0.4*p_fac[17]*f[18]+0.4472135954999579*p_fac[6]*f[18]+0.4472135954999579*f[5]*p_fac[17]+0.4472135954999579*p_fac[5]*f[17]+0.5*f[15]*p_fac[16]+0.5*p_fac[15]*f[16]+0.44721359549995804*f[10]*p_fac[14]+0.44721359549995804*p_fac[10]*f[14]+0.44721359549995804*f[10]*p_fac[13]+0.44721359549995804*p_fac[10]*f[13]+0.4*f[11]*p_fac[12]+0.44721359549995804*f[2]*p_fac[12]+0.4*p_fac[11]*f[12]+0.44721359549995804*p_fac[2]*f[12]+0.44721359549995804*f[1]*p_fac[11]+0.44721359549995804*p_fac[1]*f[11]+0.5*f[3]*p_fac[10]+0.5*p_fac[3]*f[10]+0.4472135954999579*f[4]*p_fac[8]+0.4472135954999579*p_fac[4]*f[8]+0.4472135954999579*f[4]*p_fac[7]+0.4472135954999579*p_fac[4]*f[7]+0.5*f[5]*p_fac[6]+0.5*p_fac[5]*f[6]+0.5*f[0]*p_fac[4]+0.5*p_fac[0]*f[4]+0.5*f[1]*p_fac[2]+0.5*p_fac[1]*f[2]); 
  sr_pressure[4] += (volFact/1)*(0.4472135954999579*f[19]*p_fac[19]+0.4472135954999579*f[18]*p_fac[18]+0.31943828249996997*f[17]*p_fac[17]+0.5*f[6]*p_fac[17]+0.5*p_fac[6]*f[17]+0.4472135954999579*f[15]*p_fac[15]+0.31943828249996997*f[13]*p_fac[13]+0.5000000000000001*f[3]*p_fac[13]+0.5000000000000001*p_fac[3]*f[13]+0.4472135954999579*f[12]*p_fac[12]+0.31943828249996997*f[11]*p_fac[11]+0.5000000000000001*f[2]*p_fac[11]+0.5000000000000001*p_fac[2]*f[11]+0.4472135954999579*f[10]*p_fac[10]+0.31943828249996997*f[7]*p_fac[7]+0.5*f[0]*p_fac[7]+0.5*p_fac[0]*f[7]+0.4472135954999579*f[5]*p_fac[5]+0.4472135954999579*f[4]*p_fac[4]+0.4472135954999579*f[1]*p_fac[1]); 
  sr_pressure[5] += (volFact/1)*(0.4472135954999579*f[19]*p_fac[19]+0.31943828249996997*f[18]*p_fac[18]+0.5*f[5]*p_fac[18]+0.5*p_fac[5]*f[18]+0.4472135954999579*f[17]*p_fac[17]+0.4472135954999579*f[16]*p_fac[16]+0.31943828249996997*f[14]*p_fac[14]+0.5000000000000001*f[3]*p_fac[14]+0.5000000000000001*p_fac[3]*f[14]+0.31943828249996997*f[12]*p_fac[12]+0.5000000000000001*f[1]*p_fac[12]+0.5000000000000001*p_fac[1]*f[12]+0.4472135954999579*f[11]*p_fac[11]+0.4472135954999579*f[10]*p_fac[10]+0.31943828249996997*f[8]*p_fac[8]+0.5*f[0]*p_fac[8]+0.5*p_fac[0]*f[8]+0.4472135954999579*f[6]*p_fac[6]+0.4472135954999579*f[4]*p_fac[4]+0.4472135954999579*f[2]*p_fac[2]); 
  sr_pressure[6] += (volFact/1)*(0.4472135954999579*f[15]*p_fac[19]+0.4472135954999579*p_fac[15]*f[19]+0.4*f[10]*p_fac[18]+0.4*p_fac[10]*f[18]+0.4472135954999579*f[14]*p_fac[17]+0.31943828249996997*f[13]*p_fac[17]+0.5000000000000001*f[3]*p_fac[17]+0.4472135954999579*p_fac[14]*f[17]+0.31943828249996997*p_fac[13]*f[17]+0.5000000000000001*p_fac[3]*f[17]+0.5*f[6]*p_fac[13]+0.5*p_fac[6]*f[13]+0.4*f[4]*p_fac[12]+0.4*p_fac[4]*f[12]+0.4472135954999579*f[8]*p_fac[11]+0.31943828249996997*f[7]*p_fac[11]+0.5*f[0]*p_fac[11]+0.4472135954999579*p_fac[8]*f[11]+0.31943828249996997*p_fac[7]*f[11]+0.5*p_fac[0]*f[11]+0.44721359549995804*f[5]*p_fac[10]+0.44721359549995804*p_fac[5]*f[10]+0.5000000000000001*f[2]*p_fac[7]+0.5000000000000001*p_fac[2]*f[7]+0.44721359549995804*f[1]*p_fac[4]+0.44721359549995804*p_fac[1]*f[4]); 
  sr_pressure[7] += (volFact/1)*(0.4472135954999579*f[16]*p_fac[19]+0.4472135954999579*p_fac[16]*f[19]+0.31943828249996997*f[14]*p_fac[18]+0.4472135954999579*f[13]*p_fac[18]+0.5000000000000001*f[3]*p_fac[18]+0.31943828249996997*p_fac[14]*f[18]+0.4472135954999579*p_fac[13]*f[18]+0.5000000000000001*p_fac[3]*f[18]+0.4*f[10]*p_fac[17]+0.4*p_fac[10]*f[17]+0.5*f[5]*p_fac[14]+0.5*p_fac[5]*f[14]+0.31943828249996997*f[8]*p_fac[12]+0.4472135954999579*f[7]*p_fac[12]+0.5*f[0]*p_fac[12]+0.31943828249996997*p_fac[8]*f[12]+0.4472135954999579*p_fac[7]*f[12]+0.5*p_fac[0]*f[12]+0.4*f[4]*p_fac[11]+0.4*p_fac[4]*f[11]+0.44721359549995804*f[6]*p_fac[10]+0.44721359549995804*p_fac[6]*f[10]+0.5000000000000001*f[1]*p_fac[8]+0.5000000000000001*p_fac[1]*f[8]+0.44721359549995804*f[2]*p_fac[4]+0.44721359549995804*p_fac[2]*f[4]); 

} 
