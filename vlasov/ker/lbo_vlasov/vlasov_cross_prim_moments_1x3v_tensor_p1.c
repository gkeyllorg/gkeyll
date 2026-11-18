#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_cross_prim_moments_1x3v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *alphaE, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu) 
{ 
  // A:                    Left-side matrix. 
  // rhs:                  Right-side vector. 
  // alphaE:               Morse's alpha_E factor. 
  // m_self:               self mass. 
  // moms_self:            self moments of the distribution function. 
  // prim_mom_self:        self primitive moments: mean flow velocity and thermal speed squared. 
  // m_other:              mass of other species. 
  // moms_other:           moments of the other distribution function. 
  // prim_mom_other:       cross primitive moments: mean flow velocity and thermal speed squared. 
  // boundary_corrections: corrections to momentum and energy conservation due to finite velocity space. 
  // nu:                   Cross-species collision frequency. 
 
  const double m_sumDms = (m_self+m_other)/m_self;
  const double *u_self = &prim_mom_self[0];
  const double *vtsq_self = &prim_mom_self[6];
  const double *u_other = &prim_mom_other[0];
  const double *vtsq_other = &prim_mom_other[6];
 
  double m0r[2] = {0.0}; 
  double m1r[6] = {0.0}; 
  double m2r[2] = {0.0}; 
  double cMr[6] = {0.0}; 
  double cEr[2] = {0.0}; 
  double u_selfr[6] = {0.0}; 
  double u_otherr[6] = {0.0}; 
 
  if (nu[0] > 0.0 && moms_self[8] > 0.0) { 
  
  // If a corner value is below zero, use cell average m0.
  bool notCellAvg = true;
  if (notCellAvg && (-(0.5*(2.4494897427831783*moms_self[1]-1.4142135623730951*moms_self[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*moms_self[1]+1.4142135623730951*moms_self[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(2.4494897427831783*moms_self[9]-1.4142135623730951*moms_self[8])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*moms_self[9]+1.4142135623730951*moms_self[8]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(2.4494897427831783*vtsq_self[1]-1.4142135623730951*vtsq_self[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*vtsq_self[1]+1.4142135623730951*vtsq_self[0]) < 0)) notCellAvg = false; 
 
  if (notCellAvg && (-(0.5*(2.4494897427831783*moms_other[1]-1.4142135623730951*moms_other[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*moms_other[1]+1.4142135623730951*moms_other[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(2.4494897427831783*moms_other[9]-1.4142135623730951*moms_other[8])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*moms_other[9]+1.4142135623730951*moms_other[8]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(2.4494897427831783*vtsq_other[1]-1.4142135623730951*vtsq_other[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*vtsq_other[1]+1.4142135623730951*vtsq_other[0]) < 0)) notCellAvg = false; 
 
  if (notCellAvg && (-(0.5*(2.4494897427831783*alphaE[1]-1.4142135623730951*alphaE[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*alphaE[1]+1.4142135623730951*alphaE[0]) < 0)) notCellAvg = false; 
  if (notCellAvg) { 
    m0r[0] = moms_self[0]; 
    m0r[1] = moms_self[1]; 
    m1r[0] = moms_self[2]; 
    m1r[1] = moms_self[3]; 
    m1r[2] = moms_self[4]; 
    m1r[3] = moms_self[5]; 
    m1r[4] = moms_self[6]; 
    m1r[5] = moms_self[7]; 
    m2r[0] = moms_self[8]; 
    m2r[1] = moms_self[9]; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = boundary_corrections[5]; 
    cEr[0] = boundary_corrections[6]; 
    cEr[1] = boundary_corrections[7]; 
    u_selfr[0] = u_self[0]; 
    u_selfr[1] = u_self[1]; 
    u_selfr[2] = u_self[2]; 
    u_selfr[3] = u_self[3]; 
    u_selfr[4] = u_self[4]; 
    u_selfr[5] = u_self[5]; 
    u_otherr[0] = u_other[0]; 
    u_otherr[1] = u_other[1]; 
    u_otherr[2] = u_other[2]; 
    u_otherr[3] = u_other[3]; 
    u_otherr[4] = u_other[4]; 
    u_otherr[5] = u_other[5]; 
  } else { 
    m0r[0] = moms_self[0]; 
    m0r[1] = 0.0; 
    m1r[0] = moms_self[2]; 
    m1r[1] = 0.0; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    m1r[2] = moms_self[4]; 
    m1r[3] = 0.0; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = 0.0; 
    m1r[4] = moms_self[6]; 
    m1r[5] = 0.0; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = 0.0; 
    m2r[0] = moms_self[8]; 
    m2r[1] = 0.0; 
    cEr[0] = boundary_corrections[6]; 
    cEr[1] = 0.0; 
    u_selfr[0] = u_self[0]; 
    u_selfr[1] = 0.0; 
    u_selfr[2] = 0.0; 
    u_selfr[3] = 0.0; 
    u_selfr[4] = 0.0; 
    u_selfr[5] = 0.0; 
    u_otherr[0] = u_other[0]; 
    u_otherr[1] = 0.0; 
    u_otherr[2] = 0.0; 
    u_otherr[3] = 0.0; 
    u_otherr[4] = 0.0; 
    u_otherr[5] = 0.0; 
  } 
 
  } else { 
  
    m0r[0] = 1.0; 
    m0r[1] = 0.0; 
    m1r[0] = 1.0; 
    m1r[1] = 0.0; 
    cMr[0] = 0.0; 
    cMr[1] = 0.0; 
    m1r[2] = 1.0; 
    m1r[3] = 0.0; 
    cMr[2] = 0.0; 
    cMr[3] = 0.0; 
    m1r[4] = 1.0; 
    m1r[5] = 0.0; 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    m2r[0] = 1.0; 
    m2r[1] = 0.0; 
    cEr[0] = 0.0; 
    cEr[1] = 0.0; 
    u_selfr[0] = 1.0; 
    u_selfr[1] = 0.0; 
    u_selfr[2] = 0.0; 
    u_selfr[3] = 0.0; 
    u_selfr[4] = 0.0; 
    u_selfr[5] = 0.0; 
    u_otherr[0] = 1.0; 
    u_otherr[1] = 0.0; 
    u_otherr[2] = 0.0; 
    u_otherr[3] = 0.0; 
    u_otherr[4] = 0.0; 
    u_otherr[5] = 0.0; 
  
  }
  
  double momRHS[6] = {0.0}; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossX ... // 
  gkyl_mat_set(A,0,0,1.4142135623730951*m0r[0]); 
  gkyl_mat_set(A,0,1,1.4142135623730951*m0r[1]); 
  gkyl_mat_set(A,1,0,1.4142135623730951*m0r[1]); 
  gkyl_mat_set(A,1,1,1.4142135623730951*m0r[0]); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,0,6,-(1.4142135623730951*cMr[0])); 
  gkyl_mat_set(A,0,7,-(1.4142135623730951*cMr[1])); 
  gkyl_mat_set(A,1,6,-(1.4142135623730951*cMr[1])); 
  gkyl_mat_set(A,1,7,-(1.4142135623730951*cMr[0])); 
 
  // ... Block from weak multiply of m_self, nu, m1X and uCrossX ... // 
  gkyl_mat_set(A,6,0,-(0.5*m0r[1]*u_selfr[1])-0.5*m0r[1]*u_otherr[1]-0.5*m0r[0]*u_selfr[0]-0.5*m0r[0]*u_otherr[0]+1.4142135623730951*m1r[0]); 
  gkyl_mat_set(A,6,1,-(0.5*m0r[0]*u_selfr[1])-0.5*m0r[0]*u_otherr[1]+1.4142135623730951*m1r[1]-0.5*u_selfr[0]*m0r[1]-0.5*u_otherr[0]*m0r[1]); 
  gkyl_mat_set(A,7,0,-(0.5*m0r[0]*u_selfr[1])-0.5*m0r[0]*u_otherr[1]+1.4142135623730951*m1r[1]-0.5*u_selfr[0]*m0r[1]-0.5*u_otherr[0]*m0r[1]); 
  gkyl_mat_set(A,7,1,-(0.9*m0r[1]*u_selfr[1])-0.9*m0r[1]*u_otherr[1]-0.5*m0r[0]*u_selfr[0]-0.5*m0r[0]*u_otherr[0]+1.4142135623730951*m1r[0]); 
 
  momRHS[0] += -(0.7071067811865475*((alphaE[1]*u_selfr[1]-1.0*alphaE[1]*u_otherr[1]+alphaE[0]*u_selfr[0]-1.0*alphaE[0]*u_otherr[0])*m_sumDms-2.8284271247461907*m1r[0])); 
  momRHS[1] += -(0.7071067811865475*((alphaE[0]*u_selfr[1]-1.0*alphaE[0]*u_otherr[1]+(u_selfr[0]-1.0*u_otherr[0])*alphaE[1])*m_sumDms-2.8284271247461907*m1r[1])); 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossY ... // 
  gkyl_mat_set(A,2,2,1.4142135623730951*m0r[0]); 
  gkyl_mat_set(A,2,3,1.4142135623730951*m0r[1]); 
  gkyl_mat_set(A,3,2,1.4142135623730951*m0r[1]); 
  gkyl_mat_set(A,3,3,1.4142135623730951*m0r[0]); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,2,6,-(1.4142135623730951*cMr[2])); 
  gkyl_mat_set(A,2,7,-(1.4142135623730951*cMr[3])); 
  gkyl_mat_set(A,3,6,-(1.4142135623730951*cMr[3])); 
  gkyl_mat_set(A,3,7,-(1.4142135623730951*cMr[2])); 
 
  // ... Block from weak multiply of m_self, nu, m1Y and uCrossY ... // 
  gkyl_mat_set(A,6,2,-(0.5*m0r[1]*u_selfr[3])-0.5*m0r[1]*u_otherr[3]-0.5*m0r[0]*u_selfr[2]-0.5*m0r[0]*u_otherr[2]+1.4142135623730951*m1r[2]); 
  gkyl_mat_set(A,6,3,-(0.5*m0r[0]*u_selfr[3])-0.5*m0r[0]*u_otherr[3]+1.4142135623730951*m1r[3]-0.5*m0r[1]*u_selfr[2]-0.5*m0r[1]*u_otherr[2]); 
  gkyl_mat_set(A,7,2,-(0.5*m0r[0]*u_selfr[3])-0.5*m0r[0]*u_otherr[3]+1.4142135623730951*m1r[3]-0.5*m0r[1]*u_selfr[2]-0.5*m0r[1]*u_otherr[2]); 
  gkyl_mat_set(A,7,3,-(0.9*m0r[1]*u_selfr[3])-0.9*m0r[1]*u_otherr[3]-0.5*m0r[0]*u_selfr[2]-0.5*m0r[0]*u_otherr[2]+1.4142135623730951*m1r[2]); 
 
  momRHS[2] += -(0.7071067811865475*((alphaE[1]*u_selfr[3]-1.0*alphaE[1]*u_otherr[3]+alphaE[0]*u_selfr[2]-1.0*alphaE[0]*u_otherr[2])*m_sumDms-2.8284271247461907*m1r[2])); 
  momRHS[3] += -(0.7071067811865475*((alphaE[0]*u_selfr[3]-1.0*alphaE[0]*u_otherr[3]+alphaE[1]*u_selfr[2]-1.0*alphaE[1]*u_otherr[2])*m_sumDms-2.8284271247461907*m1r[3])); 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossZ ... // 
  gkyl_mat_set(A,4,4,1.4142135623730951*m0r[0]); 
  gkyl_mat_set(A,4,5,1.4142135623730951*m0r[1]); 
  gkyl_mat_set(A,5,4,1.4142135623730951*m0r[1]); 
  gkyl_mat_set(A,5,5,1.4142135623730951*m0r[0]); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,4,6,-(1.4142135623730951*cMr[4])); 
  gkyl_mat_set(A,4,7,-(1.4142135623730951*cMr[5])); 
  gkyl_mat_set(A,5,6,-(1.4142135623730951*cMr[5])); 
  gkyl_mat_set(A,5,7,-(1.4142135623730951*cMr[4])); 
 
  // ... Block from weak multiply of m_self, nu, m1Z and uCrossZ ... // 
  gkyl_mat_set(A,6,4,-(0.5*m0r[1]*u_selfr[5])-0.5*m0r[1]*u_otherr[5]-0.5*m0r[0]*u_selfr[4]-0.5*m0r[0]*u_otherr[4]+1.4142135623730951*m1r[4]); 
  gkyl_mat_set(A,6,5,-(0.5*m0r[0]*u_selfr[5])-0.5*m0r[0]*u_otherr[5]+1.4142135623730951*m1r[5]-0.5*m0r[1]*u_selfr[4]-0.5*m0r[1]*u_otherr[4]); 
  gkyl_mat_set(A,7,4,-(0.5*m0r[0]*u_selfr[5])-0.5*m0r[0]*u_otherr[5]+1.4142135623730951*m1r[5]-0.5*m0r[1]*u_selfr[4]-0.5*m0r[1]*u_otherr[4]); 
  gkyl_mat_set(A,7,5,-(0.9*m0r[1]*u_selfr[5])-0.9*m0r[1]*u_otherr[5]-0.5*m0r[0]*u_selfr[4]-0.5*m0r[0]*u_otherr[4]+1.4142135623730951*m1r[4]); 
 
  momRHS[4] += -(0.7071067811865475*((alphaE[1]*u_selfr[5]-1.0*alphaE[1]*u_otherr[5]+alphaE[0]*u_selfr[4]-1.0*alphaE[0]*u_otherr[4])*m_sumDms-2.8284271247461907*m1r[4])); 
  momRHS[5] += -(0.7071067811865475*((alphaE[0]*u_selfr[5]-1.0*alphaE[0]*u_otherr[5]+alphaE[1]*u_selfr[4]-1.0*alphaE[1]*u_otherr[4])*m_sumDms-2.8284271247461907*m1r[5])); 
 
  double ucMSelf[2] = {0.0}; 
  double ucMOther[2] = {0.0}; 
  for (int vd=0; vd<3; vd++) 
  { 
    int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.7071067811865475*cMr[a0+1]*u_selfr[a0+1]+0.7071067811865475*cMr[a0]*u_selfr[a0]; 
    ucMSelf[1] += 0.7071067811865475*cMr[a0]*u_selfr[a0+1]+0.7071067811865475*u_selfr[a0]*cMr[a0+1]; 
    ucMOther[0] += 0.7071067811865475*cMr[a0+1]*u_otherr[a0+1]+0.7071067811865475*cMr[a0]*u_otherr[a0]; 
    ucMOther[1] += 0.7071067811865475*cMr[a0]*u_otherr[a0+1]+0.7071067811865475*u_otherr[a0]*cMr[a0+1]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  gkyl_mat_set(A,6,6,0.7071067811865475*ucMSelf[0]+0.7071067811865475*ucMOther[0]+4.242640687119286*m0r[0]-1.4142135623730951*cEr[0]); 
  gkyl_mat_set(A,6,7,0.7071067811865475*ucMSelf[1]+0.7071067811865475*ucMOther[1]+4.242640687119286*m0r[1]-1.4142135623730951*cEr[1]); 
  gkyl_mat_set(A,7,6,0.7071067811865475*ucMSelf[1]+0.7071067811865475*ucMOther[1]+4.242640687119286*m0r[1]-1.4142135623730951*cEr[1]); 
  gkyl_mat_set(A,7,7,0.7071067811865475*ucMSelf[0]+0.7071067811865475*ucMOther[0]+4.242640687119286*m0r[0]-1.4142135623730951*cEr[0]); 
 
  double uM1Self[2] = {0.0}; 
  double uM1Other[2] = {0.0}; 
  double uSumSq[2] = {0.0}; 
  for (int vd=0; vd<3; vd++) 
  { 
    int a0 = 2*vd; 
    // Dot product terms in energy equation RHS. 
    uM1Self[0] += 0.7071067811865475*m1r[a0+1]*u_selfr[a0+1]+0.7071067811865475*m1r[a0]*u_selfr[a0]; 
    uM1Self[1] += 0.7071067811865475*m1r[a0]*u_selfr[a0+1]+0.7071067811865475*u_selfr[a0]*m1r[a0+1]; 
    uM1Other[0] += 0.7071067811865475*m1r[a0+1]*u_otherr[a0+1]+0.7071067811865475*m1r[a0]*u_otherr[a0]; 
    uM1Other[1] += 0.7071067811865475*m1r[a0]*u_otherr[a0+1]+0.7071067811865475*u_otherr[a0]*m1r[a0+1]; 
  const double u_selfr0R2 = pow(u_selfr[a0],2);
  const double u_selfr1R2 = pow(u_selfr[a0+1],2);
  const double u_otherr0R2 = pow(u_otherr[a0],2);
  const double u_otherr1R2 = pow(u_otherr[a0+1],2);

  uSumSq[0] += 0.7071067811865475*u_selfr1R2-1.4142135623730951*u_otherr[a0+1]*u_selfr[a0+1]+0.7071067811865475*u_otherr1R2+0.7071067811865475*u_selfr0R2-1.4142135623730951*u_otherr[a0]*u_selfr[a0]+0.7071067811865475*u_otherr0R2; 
  uSumSq[1] += 1.4142135623730951*u_selfr[a0]*u_selfr[a0+1]-1.4142135623730951*u_otherr[a0]*u_selfr[a0+1]-1.4142135623730951*u_selfr[a0]*u_otherr[a0+1]+1.4142135623730951*u_otherr[a0]*u_otherr[a0+1]; 
  } 
 
  double m_diff = m_other-m_self;
  double enRHS[2] = {0.0}; 
  enRHS[0] = (2.1213203435596424*alphaE[1]*vtsq_other[1]*m_other)/m_self+(2.1213203435596424*alphaE[0]*vtsq_other[0]*m_other)/m_self+(0.3535533905932737*alphaE[1]*uSumSq[1]*m_diff)/m_self+(0.3535533905932737*alphaE[0]*uSumSq[0]*m_diff)/m_self-2.1213203435596424*alphaE[1]*vtsq_self[1]-2.1213203435596424*alphaE[0]*vtsq_self[0]-1.0*uM1Self[0]-1.0*uM1Other[0]+4.0*m2r[0]; 
  enRHS[1] = (2.1213203435596424*alphaE[0]*vtsq_other[1]*m_other)/m_self+(2.1213203435596424*vtsq_other[0]*alphaE[1]*m_other)/m_self+(0.3535533905932737*alphaE[0]*uSumSq[1]*m_diff)/m_self+(0.3535533905932737*uSumSq[0]*alphaE[1]*m_diff)/m_self-2.1213203435596424*alphaE[0]*vtsq_self[1]-1.0*uM1Self[1]-1.0*uM1Other[1]+4.0*m2r[1]-2.1213203435596424*vtsq_self[0]*alphaE[1]; 
 
  gkyl_mat_set(rhs,0,0,momRHS[0]); 
  gkyl_mat_set(rhs,1,0,momRHS[1]); 
  gkyl_mat_set(rhs,2,0,momRHS[2]); 
  gkyl_mat_set(rhs,3,0,momRHS[3]); 
  gkyl_mat_set(rhs,4,0,momRHS[4]); 
  gkyl_mat_set(rhs,5,0,momRHS[5]); 
  gkyl_mat_set(rhs,6,0,enRHS[0]); 
  gkyl_mat_set(rhs,7,0,enRHS[1]); 
} 
 
