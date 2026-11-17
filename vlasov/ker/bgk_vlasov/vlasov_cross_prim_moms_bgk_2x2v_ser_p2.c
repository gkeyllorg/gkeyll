#include <gkyl_vlasov_cross_prim_moms_bgk_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
 
GKYL_CU_DH void vlasov_cross_prim_moms_bgk_2x2v_ser_p2(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross) 
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
  const double *u_s = &prim_moms_self[8]; 
  const double *vtsq_s = &prim_moms_self[24]; 
  const double *n_r = &prim_moms_other[0]; 
  const double *u_r = &prim_moms_other[8]; 
  const double *vtsq_r = &prim_moms_other[24]; 
 
  double *n_sr = &prim_moms_cross[0]; 
  double *u_sr = &prim_moms_cross[8]; 
  double *vtsq_sr = &prim_moms_cross[24]; 

  n_sr[0] = n_s[0]; 
  n_sr[1] = n_s[1]; 
  n_sr[2] = n_s[2]; 
  n_sr[3] = n_s[3]; 
  n_sr[4] = n_s[4]; 
  n_sr[5] = n_s[5]; 
  n_sr[6] = n_s[6]; 
  n_sr[7] = n_s[7]; 
 
  const double deltaTbetap1 = delta_sr*betaGreenep1;
 
  double dU[8]; 
  double dUSq[8] = {0.0}; 
  double dUSq_sum[8] = {0.0}; 
  dU[0] = u_r[0] - u_s[0]; 
  dU[1] = u_r[1] - u_s[1]; 
  dU[2] = u_r[2] - u_s[2]; 
  dU[3] = u_r[3] - u_s[3]; 
  dU[4] = u_r[4] - u_s[4]; 
  dU[5] = u_r[5] - u_s[5]; 
  dU[6] = u_r[6] - u_s[6]; 
  dU[7] = u_r[7] - u_s[7]; 
  binop_mul_2d_ser_p2(dU, dU, dUSq); 
  dUSq_sum[0] += dUSq[0]; 
  dUSq_sum[1] += dUSq[1]; 
  dUSq_sum[2] += dUSq[2]; 
  dUSq_sum[3] += dUSq[3]; 
  dUSq_sum[4] += dUSq[4]; 
  dUSq_sum[5] += dUSq[5]; 
  dUSq_sum[6] += dUSq[6]; 
  dUSq_sum[7] += dUSq[7]; 
  u_sr[0] = u_s[0] + 0.5*deltaTbetap1*dU[0]; 
  u_sr[1] = u_s[1] + 0.5*deltaTbetap1*dU[1]; 
  u_sr[2] = u_s[2] + 0.5*deltaTbetap1*dU[2]; 
  u_sr[3] = u_s[3] + 0.5*deltaTbetap1*dU[3]; 
  u_sr[4] = u_s[4] + 0.5*deltaTbetap1*dU[4]; 
  u_sr[5] = u_s[5] + 0.5*deltaTbetap1*dU[5]; 
  u_sr[6] = u_s[6] + 0.5*deltaTbetap1*dU[6]; 
  u_sr[7] = u_s[7] + 0.5*deltaTbetap1*dU[7]; 
 
  dU[0] = u_r[8] - u_s[8]; 
  dU[1] = u_r[9] - u_s[9]; 
  dU[2] = u_r[10] - u_s[10]; 
  dU[3] = u_r[11] - u_s[11]; 
  dU[4] = u_r[12] - u_s[12]; 
  dU[5] = u_r[13] - u_s[13]; 
  dU[6] = u_r[14] - u_s[14]; 
  dU[7] = u_r[15] - u_s[15]; 
  binop_mul_2d_ser_p2(dU, dU, dUSq); 
  dUSq_sum[0] += dUSq[0]; 
  dUSq_sum[1] += dUSq[1]; 
  dUSq_sum[2] += dUSq[2]; 
  dUSq_sum[3] += dUSq[3]; 
  dUSq_sum[4] += dUSq[4]; 
  dUSq_sum[5] += dUSq[5]; 
  dUSq_sum[6] += dUSq[6]; 
  dUSq_sum[7] += dUSq[7]; 
  u_sr[8] = u_s[8] + 0.5*deltaTbetap1*dU[0]; 
  u_sr[9] = u_s[9] + 0.5*deltaTbetap1*dU[1]; 
  u_sr[10] = u_s[10] + 0.5*deltaTbetap1*dU[2]; 
  u_sr[11] = u_s[11] + 0.5*deltaTbetap1*dU[3]; 
  u_sr[12] = u_s[12] + 0.5*deltaTbetap1*dU[4]; 
  u_sr[13] = u_s[13] + 0.5*deltaTbetap1*dU[5]; 
  u_sr[14] = u_s[14] + 0.5*deltaTbetap1*dU[6]; 
  u_sr[15] = u_s[15] + 0.5*deltaTbetap1*dU[7]; 
 
  const double vdimPhys = 2; 
  const double mSum = m_self + m_other;
 
  double dTemp[8] = {0.0}; 
  dTemp[0] = vdimPhys * (m_r*vtsq_r[0]-m_s*vtsq_s[0]); 
  dTemp[1] = vdimPhys * (m_r*vtsq_r[1]-m_s*vtsq_s[1]); 
  dTemp[2] = vdimPhys * (m_r*vtsq_r[2]-m_s*vtsq_s[2]); 
  dTemp[3] = vdimPhys * (m_r*vtsq_r[3]-m_s*vtsq_s[3]); 
  dTemp[4] = vdimPhys * (m_r*vtsq_r[4]-m_s*vtsq_s[4]); 
  dTemp[5] = vdimPhys * (m_r*vtsq_r[5]-m_s*vtsq_s[5]); 
  dTemp[6] = vdimPhys * (m_r*vtsq_r[6]-m_s*vtsq_s[6]); 
  dTemp[7] = vdimPhys * (m_r*vtsq_r[7]-m_s*vtsq_s[7]); 
  const double dUSq_fac = (m_r-0.25*deltaTbetap1*mSum/vdimPhys);
  vtsq_sr[0] = vtsq_s[0] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[0]+dUSq_fac*dUSq[0]); 
  vtsq_sr[1] = vtsq_s[1] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[1]+dUSq_fac*dUSq[1]); 
  vtsq_sr[2] = vtsq_s[2] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[2]+dUSq_fac*dUSq[2]); 
  vtsq_sr[3] = vtsq_s[3] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[3]+dUSq_fac*dUSq[3]); 
  vtsq_sr[4] = vtsq_s[4] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[4]+dUSq_fac*dUSq[4]); 
  vtsq_sr[5] = vtsq_s[5] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[5]+dUSq_fac*dUSq[5]); 
  vtsq_sr[6] = vtsq_s[6] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[6]+dUSq_fac*dUSq[6]); 
  vtsq_sr[7] = vtsq_s[7] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[7]+dUSq_fac*dUSq[7]); 
 
  // If vtsq_sr is negative at a corner, turn off collisions.
  bool negative_cross_temp = false; 
  if (-(0.5*(3.872983346207417*(vtsq_sr[7]+vtsq_sr[6])-2.23606797749979*(vtsq_sr[5]+vtsq_sr[4])-3.0*vtsq_sr[3]+1.7320508075688772*(vtsq_sr[2]+vtsq_sr[1])-1.0*vtsq_sr[0])) < 0.0) negative_cross_temp = true; 
  if (0.5*(3.872983346207417*vtsq_sr[7]-3.872983346207417*vtsq_sr[6]+2.23606797749979*(vtsq_sr[5]+vtsq_sr[4])-3.0*vtsq_sr[3]-1.7320508075688772*vtsq_sr[2]+1.7320508075688772*vtsq_sr[1]+vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (-(0.5*(3.872983346207417*vtsq_sr[7]-3.872983346207417*vtsq_sr[6]-2.23606797749979*(vtsq_sr[5]+vtsq_sr[4])+3.0*vtsq_sr[3]-1.7320508075688772*vtsq_sr[2]+1.7320508075688772*vtsq_sr[1]-1.0*vtsq_sr[0])) < 0.0) negative_cross_temp = true; 
  if (0.5*(3.872983346207417*(vtsq_sr[7]+vtsq_sr[6])+2.23606797749979*(vtsq_sr[5]+vtsq_sr[4])+3.0*vtsq_sr[3]+1.7320508075688772*(vtsq_sr[2]+vtsq_sr[1])+vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (negative_cross_temp) { 
    u_sr[0] = u_s[0]; 
    u_sr[1] = u_s[1]; 
    u_sr[2] = u_s[2]; 
    u_sr[3] = u_s[3]; 
    u_sr[4] = u_s[4]; 
    u_sr[5] = u_s[5]; 
    u_sr[6] = u_s[6]; 
    u_sr[7] = u_s[7]; 
    u_sr[8] = u_s[8]; 
    u_sr[9] = u_s[9]; 
    u_sr[10] = u_s[10]; 
    u_sr[11] = u_s[11]; 
    u_sr[12] = u_s[12]; 
    u_sr[13] = u_s[13]; 
    u_sr[14] = u_s[14]; 
    u_sr[15] = u_s[15]; 
    vtsq_sr[0] = vtsq_s[0]; 
    vtsq_sr[1] = vtsq_s[1]; 
    vtsq_sr[2] = vtsq_s[2]; 
    vtsq_sr[3] = vtsq_s[3]; 
    vtsq_sr[4] = vtsq_s[4]; 
    vtsq_sr[5] = vtsq_s[5]; 
    vtsq_sr[6] = vtsq_s[6]; 
    vtsq_sr[7] = vtsq_s[7]; 
  } 
} 
