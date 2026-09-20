#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *alphaE, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu) 
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
  const double *vtsq_self = &prim_mom_self[2];
  const double *u_other = &prim_mom_other[0];
  const double *vtsq_other = &prim_mom_other[2];
 
  double m0r[2] = {0.0}; 
  double m1r[2] = {0.0}; 
  double m2r[2] = {0.0}; 
  double cMr[2] = {0.0}; 
  double cEr[2] = {0.0}; 
  double u_selfr[2] = {0.0}; 
  double u_otherr[2] = {0.0}; 
 
  if (nu[0] > 0.0 && moms_self[4] > 0.0) { 
  
  // If a corner value is below zero, use cell average m0.
  bool notCellAvg = true;
  if (notCellAvg && (-(0.5*(2.4494897427831783*moms_self[1]-1.4142135623730951*moms_self[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*moms_self[1]+1.4142135623730951*moms_self[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(2.4494897427831783*moms_self[5]-1.4142135623730951*moms_self[4])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*moms_self[5]+1.4142135623730951*moms_self[4]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(2.4494897427831783*vtsq_self[1]-1.4142135623730951*vtsq_self[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*vtsq_self[1]+1.4142135623730951*vtsq_self[0]) < 0)) notCellAvg = false; 
 
  if (notCellAvg && (-(0.5*(2.4494897427831783*moms_other[1]-1.4142135623730951*moms_other[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*moms_other[1]+1.4142135623730951*moms_other[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(2.4494897427831783*moms_other[5]-1.4142135623730951*moms_other[4])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*moms_other[5]+1.4142135623730951*moms_other[4]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(2.4494897427831783*vtsq_other[1]-1.4142135623730951*vtsq_other[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*vtsq_other[1]+1.4142135623730951*vtsq_other[0]) < 0)) notCellAvg = false; 
 
  if (notCellAvg && (-(0.5*(2.4494897427831783*alphaE[1]-1.4142135623730951*alphaE[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*alphaE[1]+1.4142135623730951*alphaE[0]) < 0)) notCellAvg = false; 
  if (notCellAvg) { 
    m0r[0] = moms_self[0]; 
    m0r[1] = moms_self[1]; 
    m1r[0] = moms_self[2]; 
    m1r[1] = moms_self[3]; 
    m2r[0] = moms_self[4]; 
    m2r[1] = moms_self[5]; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cEr[0] = boundary_corrections[2]; 
    cEr[1] = boundary_corrections[3]; 
    u_selfr[0] = u_self[0]; 
    u_selfr[1] = u_self[1]; 
    u_otherr[0] = u_other[0]; 
    u_otherr[1] = u_other[1]; 
  } else { 
    m0r[0] = moms_self[0]; 
    m0r[1] = 0.0; 
    m1r[0] = moms_self[2]; 
    m1r[1] = 0.0; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    m2r[0] = moms_self[4]; 
    m2r[1] = 0.0; 
    cEr[0] = boundary_corrections[2]; 
    cEr[1] = 0.0; 
    u_selfr[0] = u_self[0]; 
    u_selfr[1] = 0.0; 
    u_otherr[0] = u_other[0]; 
    u_otherr[1] = 0.0; 
  } 
 
  } else { 
  
    m0r[0] = 1.0; 
    m0r[1] = 0.0; 
    m1r[0] = 1.0; 
    m1r[1] = 0.0; 
    cMr[0] = 0.0; 
    cMr[1] = 0.0; 
    m2r[0] = 1.0; 
    m2r[1] = 0.0; 
    cEr[0] = 0.0; 
    cEr[1] = 0.0; 
    u_selfr[0] = 1.0; 
    u_selfr[1] = 0.0; 
    u_otherr[0] = 1.0; 
    u_otherr[1] = 0.0; 
  
  }
  
  double momRHS[2] = {0.0}; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossX ... // 
  gkyl_mat_set(A,0,0,1.4142135623730951*m0r[0]); 
  gkyl_mat_set(A,0,1,1.4142135623730951*m0r[1]); 
  gkyl_mat_set(A,1,0,1.4142135623730951*m0r[1]); 
  gkyl_mat_set(A,1,1,1.4142135623730951*m0r[0]); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,0,2,-(1.4142135623730951*cMr[0])); 
  gkyl_mat_set(A,0,3,-(1.4142135623730951*cMr[1])); 
  gkyl_mat_set(A,1,2,-(1.4142135623730951*cMr[1])); 
  gkyl_mat_set(A,1,3,-(1.4142135623730951*cMr[0])); 
 
  // ... Block from weak multiply of m_self, nu, m1X and uCrossX ... // 
  gkyl_mat_set(A,2,0,-(0.5*m0r[1]*u_selfr[1])-0.5*m0r[1]*u_otherr[1]-0.5*m0r[0]*u_selfr[0]-0.5*m0r[0]*u_otherr[0]+1.4142135623730951*m1r[0]); 
  gkyl_mat_set(A,2,1,-(0.5*m0r[0]*u_selfr[1])-0.5*m0r[0]*u_otherr[1]+1.4142135623730951*m1r[1]-0.5*u_selfr[0]*m0r[1]-0.5*u_otherr[0]*m0r[1]); 
  gkyl_mat_set(A,3,0,-(0.5*m0r[0]*u_selfr[1])-0.5*m0r[0]*u_otherr[1]+1.4142135623730951*m1r[1]-0.5*u_selfr[0]*m0r[1]-0.5*u_otherr[0]*m0r[1]); 
  gkyl_mat_set(A,3,1,-(0.9*m0r[1]*u_selfr[1])-0.9*m0r[1]*u_otherr[1]-0.5*m0r[0]*u_selfr[0]-0.5*m0r[0]*u_otherr[0]+1.4142135623730951*m1r[0]); 
 
  momRHS[0] += -(0.7071067811865475*((alphaE[1]*u_selfr[1]-1.0*alphaE[1]*u_otherr[1]+alphaE[0]*u_selfr[0]-1.0*alphaE[0]*u_otherr[0])*m_sumDms-2.8284271247461907*m1r[0])); 
  momRHS[1] += -(0.7071067811865475*((alphaE[0]*u_selfr[1]-1.0*alphaE[0]*u_otherr[1]+(u_selfr[0]-1.0*u_otherr[0])*alphaE[1])*m_sumDms-2.8284271247461907*m1r[1])); 
 
  double ucMSelf[2] = {0.0}; 
  double ucMOther[2] = {0.0}; 
  for (int vd=0; vd<1; vd++) 
  { 
    int a0 = 2*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.7071067811865475*cMr[a0+1]*u_selfr[a0+1]+0.7071067811865475*cMr[a0]*u_selfr[a0]; 
    ucMSelf[1] += 0.7071067811865475*cMr[a0]*u_selfr[a0+1]+0.7071067811865475*u_selfr[a0]*cMr[a0+1]; 
    ucMOther[0] += 0.7071067811865475*cMr[a0+1]*u_otherr[a0+1]+0.7071067811865475*cMr[a0]*u_otherr[a0]; 
    ucMOther[1] += 0.7071067811865475*cMr[a0]*u_otherr[a0+1]+0.7071067811865475*u_otherr[a0]*cMr[a0+1]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  gkyl_mat_set(A,2,2,0.7071067811865475*ucMSelf[0]+0.7071067811865475*ucMOther[0]+1.4142135623730951*m0r[0]-1.4142135623730951*cEr[0]); 
  gkyl_mat_set(A,2,3,0.7071067811865475*ucMSelf[1]+0.7071067811865475*ucMOther[1]+1.4142135623730951*m0r[1]-1.4142135623730951*cEr[1]); 
  gkyl_mat_set(A,3,2,0.7071067811865475*ucMSelf[1]+0.7071067811865475*ucMOther[1]+1.4142135623730951*m0r[1]-1.4142135623730951*cEr[1]); 
  gkyl_mat_set(A,3,3,0.7071067811865475*ucMSelf[0]+0.7071067811865475*ucMOther[0]+1.4142135623730951*m0r[0]-1.4142135623730951*cEr[0]); 
 
  double uM1Self[2] = {0.0}; 
  double uM1Other[2] = {0.0}; 
  double uSumSq[2] = {0.0}; 
  for (int vd=0; vd<1; vd++) 
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
  enRHS[0] = (0.7071067811865475*alphaE[1]*vtsq_other[1]*m_other)/m_self+(0.7071067811865475*alphaE[0]*vtsq_other[0]*m_other)/m_self+(0.3535533905932737*alphaE[1]*uSumSq[1]*m_diff)/m_self+(0.3535533905932737*alphaE[0]*uSumSq[0]*m_diff)/m_self-0.7071067811865475*alphaE[1]*vtsq_self[1]-0.7071067811865475*alphaE[0]*vtsq_self[0]-1.0*uM1Self[0]-1.0*uM1Other[0]+4.0*m2r[0]; 
  enRHS[1] = (0.7071067811865475*alphaE[0]*vtsq_other[1]*m_other)/m_self+(0.7071067811865475*vtsq_other[0]*alphaE[1]*m_other)/m_self+(0.3535533905932737*alphaE[0]*uSumSq[1]*m_diff)/m_self+(0.3535533905932737*uSumSq[0]*alphaE[1]*m_diff)/m_self-0.7071067811865475*alphaE[0]*vtsq_self[1]-1.0*uM1Self[1]-1.0*uM1Other[1]+4.0*m2r[1]-0.7071067811865475*vtsq_self[0]*alphaE[1]; 
 
  gkyl_mat_set(rhs,0,0,momRHS[0]); 
  gkyl_mat_set(rhs,1,0,momRHS[1]); 
  gkyl_mat_set(rhs,2,0,enRHS[0]); 
  gkyl_mat_set(rhs,3,0,enRHS[1]); 
} 
 
