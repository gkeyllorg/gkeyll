#include <gkyl_gyrokinetic_cross_prim_moms_bgk_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
 
GKYL_CU_DH void gyrokinetic_cross_prim_moms_bgk_3x2v_ser_p1(const double delta_sr, const double betaGreenep1, const double m_self, const double *prim_moms_self, const double m_other, const double *prim_moms_other, double *prim_moms_cross) 
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
  const double *upar_s = &prim_moms_self[8]; 
  const double *vtsq_s = &prim_moms_self[16]; 
  const double *n_r = &prim_moms_other[0]; 
  const double *upar_r = &prim_moms_other[8]; 
  const double *vtsq_r = &prim_moms_other[16]; 
 
  double *n_sr = &prim_moms_cross[0]; 
  double *upar_sr = &prim_moms_cross[8]; 
  double *vtsq_sr = &prim_moms_cross[16]; 

  n_sr[0] = n_s[0]; 
  n_sr[1] = n_s[1]; 
  n_sr[2] = n_s[2]; 
  n_sr[3] = n_s[3]; 
  n_sr[4] = n_s[4]; 
  n_sr[5] = n_s[5]; 
  n_sr[6] = n_s[6]; 
  n_sr[7] = n_s[7]; 
 
  const double deltaTbetap1 = delta_sr*betaGreenep1;
 
  double dUpar[8]; 
  dUpar[0] = upar_r[0] - upar_s[0]; 
  dUpar[1] = upar_r[1] - upar_s[1]; 
  dUpar[2] = upar_r[2] - upar_s[2]; 
  dUpar[3] = upar_r[3] - upar_s[3]; 
  dUpar[4] = upar_r[4] - upar_s[4]; 
  dUpar[5] = upar_r[5] - upar_s[5]; 
  dUpar[6] = upar_r[6] - upar_s[6]; 
  dUpar[7] = upar_r[7] - upar_s[7]; 
  upar_sr[0] = upar_s[0] + 0.5*deltaTbetap1*dUpar[0]; 
  upar_sr[1] = upar_s[1] + 0.5*deltaTbetap1*dUpar[1]; 
  upar_sr[2] = upar_s[2] + 0.5*deltaTbetap1*dUpar[2]; 
  upar_sr[3] = upar_s[3] + 0.5*deltaTbetap1*dUpar[3]; 
  upar_sr[4] = upar_s[4] + 0.5*deltaTbetap1*dUpar[4]; 
  upar_sr[5] = upar_s[5] + 0.5*deltaTbetap1*dUpar[5]; 
  upar_sr[6] = upar_s[6] + 0.5*deltaTbetap1*dUpar[6]; 
  upar_sr[7] = upar_s[7] + 0.5*deltaTbetap1*dUpar[7]; 
 
  const double vdimPhys = 3; 
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
  double dUparSq[8] = {0.0}; 
  binop_mul_3d_ser_p1(dUpar, dUpar, dUparSq); 
  const double dUparSq_fac = (m_r-0.25*deltaTbetap1*mSum);
  vtsq_sr[0] = vtsq_s[0] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[0]+dUparSq_fac*dUparSq[0]); 
  vtsq_sr[1] = vtsq_s[1] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[1]+dUparSq_fac*dUparSq[1]); 
  vtsq_sr[2] = vtsq_s[2] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[2]+dUparSq_fac*dUparSq[2]); 
  vtsq_sr[3] = vtsq_s[3] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[3]+dUparSq_fac*dUparSq[3]); 
  vtsq_sr[4] = vtsq_s[4] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[4]+dUparSq_fac*dUparSq[4]); 
  vtsq_sr[5] = vtsq_s[5] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[5]+dUparSq_fac*dUparSq[5]); 
  vtsq_sr[6] = vtsq_s[6] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[6]+dUparSq_fac*dUparSq[6]); 
  vtsq_sr[7] = vtsq_s[7] + deltaTbetap1/(vdimPhys*mSum)*(dTemp[7]+dUparSq_fac*dUparSq[7]); 
 
  // If vtsq_sr is negative at a corner, turn off collisions.
  bool negative_cross_temp = false; 
  if (-(0.25*(7.348469228349534*vtsq_sr[7]-4.242640687119286*(vtsq_sr[6]+vtsq_sr[5]+vtsq_sr[4])+2.4494897427831783*(vtsq_sr[3]+vtsq_sr[2]+vtsq_sr[1])-1.4142135623730951*vtsq_sr[0])) < 0.0) negative_cross_temp = true; 
  if (0.25*(7.348469228349534*vtsq_sr[7]+4.242640687119286*vtsq_sr[6]-4.242640687119286*(vtsq_sr[5]+vtsq_sr[4])-2.4494897427831783*(vtsq_sr[3]+vtsq_sr[2])+2.4494897427831783*vtsq_sr[1]+1.4142135623730951*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (0.25*(7.348469228349534*vtsq_sr[7]-4.242640687119286*vtsq_sr[6]+4.242640687119286*vtsq_sr[5]-4.242640687119286*vtsq_sr[4]-2.4494897427831783*vtsq_sr[3]+2.4494897427831783*vtsq_sr[2]-2.4494897427831783*vtsq_sr[1]+1.4142135623730951*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (-(0.25*(7.348469228349534*vtsq_sr[7]+4.242640687119286*(vtsq_sr[6]+vtsq_sr[5])-4.242640687119286*vtsq_sr[4]+2.4494897427831783*vtsq_sr[3]-2.4494897427831783*(vtsq_sr[2]+vtsq_sr[1])-1.4142135623730951*vtsq_sr[0])) < 0.0) negative_cross_temp = true; 
  if (0.25*(7.348469228349534*vtsq_sr[7]-4.242640687119286*(vtsq_sr[6]+vtsq_sr[5])+4.242640687119286*vtsq_sr[4]+2.4494897427831783*vtsq_sr[3]-2.4494897427831783*(vtsq_sr[2]+vtsq_sr[1])+1.4142135623730951*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (-(0.25*(7.348469228349534*vtsq_sr[7]+4.242640687119286*vtsq_sr[6]-4.242640687119286*vtsq_sr[5]+4.242640687119286*vtsq_sr[4]-2.4494897427831783*vtsq_sr[3]+2.4494897427831783*vtsq_sr[2]-2.4494897427831783*vtsq_sr[1]-1.4142135623730951*vtsq_sr[0])) < 0.0) negative_cross_temp = true; 
  if (-(0.25*(7.348469228349534*vtsq_sr[7]-4.242640687119286*vtsq_sr[6]+4.242640687119286*(vtsq_sr[5]+vtsq_sr[4])-2.4494897427831783*(vtsq_sr[3]+vtsq_sr[2])+2.4494897427831783*vtsq_sr[1]-1.4142135623730951*vtsq_sr[0])) < 0.0) negative_cross_temp = true; 
  if (0.25*(7.348469228349534*vtsq_sr[7]+4.242640687119286*(vtsq_sr[6]+vtsq_sr[5]+vtsq_sr[4])+2.4494897427831783*(vtsq_sr[3]+vtsq_sr[2]+vtsq_sr[1])+1.4142135623730951*vtsq_sr[0]) < 0.0) negative_cross_temp = true; 
  if (negative_cross_temp) { 
    upar_sr[0] = upar_s[0]; 
    vtsq_sr[0] = vtsq_s[0]; 
    upar_sr[1] = upar_s[1]; 
    vtsq_sr[1] = vtsq_s[1]; 
    upar_sr[2] = upar_s[2]; 
    vtsq_sr[2] = vtsq_s[2]; 
    upar_sr[3] = upar_s[3]; 
    vtsq_sr[3] = vtsq_s[3]; 
    upar_sr[4] = upar_s[4]; 
    vtsq_sr[4] = vtsq_s[4]; 
    upar_sr[5] = upar_s[5]; 
    vtsq_sr[5] = vtsq_s[5]; 
    upar_sr[6] = upar_s[6]; 
    vtsq_sr[6] = vtsq_s[6]; 
    upar_sr[7] = upar_s[7]; 
    vtsq_sr[7] = vtsq_s[7]; 
  } 
} 
