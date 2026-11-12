#include <gkyl_vlasov_cross_prim_moms_bgk_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
 
GKYL_CU_DH void vlasov_cross_prim_moms_bgk_1x1v_ser_p1(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross) 
{ 
  // delta_sr: Free parameter in relation between nu_sr and alpha_E. 
  // betaGreenep1: Greene's beta + 1. 
  // m_self: mass of this species. 
  // prim_moms_self: primitive moments of this species. 
  // m_other: mass of other species. 
  // prim_moms_other: primitive moments of other species. 
  // prim_moms_cross: cross primitive moments. 
 
  const double m_s = m_self; 
  const double m_r = m_other; 
  const double *n_s = &prim_moms_self[0]; 
  const double *u_s = &prim_moms_self[2]; 
  const double *vtsq_s = &prim_moms_self[4]; 
  const double *n_r = &prim_moms_other[0]; 
  const double *u_r = &prim_moms_other[2]; 
  const double *vtsq_r = &prim_moms_other[4]; 
 
  double *n_sr = &prim_moms_cross[0]; 
  double *u_sr = &prim_moms_cross[2]; 
  double *vtsq_sr = &prim_moms_cross[4]; 

  n_sr[0] = n_s[0]; 
  n_sr[1] = n_s[1]; 
 
  const double deltaTbetap1 = delta_sr*betaGreenep1;
 
  double dU[2]; 
  double dUSq[2] = {0.0}; 
  double dUSq_sum[2] = {0.0}; 
  dU[0] = u_r[0] - u_s[0]; 
  dU[1] = u_r[1] - u_s[1]; 
  binop_mul_1d_ser_p1(dU, dU, dUSq); 
  dUSq_sum[0] += dUSq[0]; 
  dUSq_sum[1] += dUSq[1]; 
  u_sr[0] = u_s[0] + 0.5*deltaTbetap1*dU[0]; 
  u_sr[1] = u_s[1] + 0.5*deltaTbetap1*dU[1]; 
 
  const double vdimPhys = 1; 
  const double mSum = m_self + m_other;
 
  double dTemp[2] = {0.0}; 
  dTemp[0] = vdimPhys * (m_r*vtsq_r[0]-m_s*vtsq_s[0]); 
  dTemp[1] = vdimPhys * (m_r*vtsq_r[1]-m_s*vtsq_s[1]); 
  const double dUSq_fac = (m_r-0.25*deltaTbetap1*mSum/vdimPhys);
  vtsq_sr[0] = vtsq_s[0] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[0]+dUSq_fac*dUSq[0]); 
  vtsq_sr[1] = vtsq_s[1] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[1]+dUSq_fac*dUSq[1]); 
 
  // If vtsq_sr is negative at a corner, turn off collisions.
  bool negative_cross_temp = false; 
  if (-(0.5*(2.4494897427831783*vtsq_sr[1]-1.4142135623730951*vtsq_sr[0])) < 0.0) negative_cross_temp = true; 
  if (0.5*(2.4494897427831783*vtsq_sr[1]+1.4142135623730951*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (negative_cross_temp) { 
    u_sr[0] = u_s[0]; 
    u_sr[1] = u_s[1]; 
    vtsq_sr[0] = vtsq_s[0]; 
    vtsq_sr[1] = vtsq_s[1]; 
  } 
} 
