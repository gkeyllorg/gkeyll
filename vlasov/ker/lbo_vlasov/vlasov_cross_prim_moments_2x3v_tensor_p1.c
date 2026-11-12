#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_cross_prim_moments_2x3v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *alphaE, const double m_self, const double *moms_self, const double *prim_mom_self, const double m_other, const double *moms_other, const double *prim_mom_other, const double *boundary_corrections, const double *nu) 
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
  const double *vtsq_self = &prim_mom_self[12];
  const double *u_other = &prim_mom_other[0];
  const double *vtsq_other = &prim_mom_other[12];
 
  double m0r[4] = {0.0}; 
  double m1r[12] = {0.0}; 
  double m2r[4] = {0.0}; 
  double cMr[12] = {0.0}; 
  double cEr[4] = {0.0}; 
  double u_selfr[12] = {0.0}; 
  double u_otherr[12] = {0.0}; 
 
  if (nu[0] > 0.0 && moms_self[16] > 0.0) { 
  
  // If a corner value is below zero, use cell average m0.
  bool notCellAvg = true;
  if (notCellAvg && (0.5*(3.0*moms_self[3]-1.7320508075688772*(moms_self[2]+moms_self[1])+moms_self[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms_self[3]+1.7320508075688772*moms_self[2]-1.7320508075688772*moms_self[1]-1.0*moms_self[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms_self[3]-1.7320508075688772*moms_self[2]+1.7320508075688772*moms_self[1]-1.0*moms_self[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*moms_self[3]+1.7320508075688772*(moms_self[2]+moms_self[1])+moms_self[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*moms_self[19]-1.7320508075688772*(moms_self[18]+moms_self[17])+moms_self[16]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms_self[19]+1.7320508075688772*moms_self[18]-1.7320508075688772*moms_self[17]-1.0*moms_self[16])) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms_self[19]-1.7320508075688772*moms_self[18]+1.7320508075688772*moms_self[17]-1.0*moms_self[16])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*moms_self[19]+1.7320508075688772*(moms_self[18]+moms_self[17])+moms_self[16]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*vtsq_self[3]-1.7320508075688772*(vtsq_self[2]+vtsq_self[1])+vtsq_self[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*vtsq_self[3]+1.7320508075688772*vtsq_self[2]-1.7320508075688772*vtsq_self[1]-1.0*vtsq_self[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*vtsq_self[3]-1.7320508075688772*vtsq_self[2]+1.7320508075688772*vtsq_self[1]-1.0*vtsq_self[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*vtsq_self[3]+1.7320508075688772*(vtsq_self[2]+vtsq_self[1])+vtsq_self[0]) < 0)) notCellAvg = false; 
 
  if (notCellAvg && (0.5*(3.0*moms_other[3]-1.7320508075688772*(moms_other[2]+moms_other[1])+moms_other[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms_other[3]+1.7320508075688772*moms_other[2]-1.7320508075688772*moms_other[1]-1.0*moms_other[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms_other[3]-1.7320508075688772*moms_other[2]+1.7320508075688772*moms_other[1]-1.0*moms_other[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*moms_other[3]+1.7320508075688772*(moms_other[2]+moms_other[1])+moms_other[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*moms_other[19]-1.7320508075688772*(moms_other[18]+moms_other[17])+moms_other[16]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms_other[19]+1.7320508075688772*moms_other[18]-1.7320508075688772*moms_other[17]-1.0*moms_other[16])) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms_other[19]-1.7320508075688772*moms_other[18]+1.7320508075688772*moms_other[17]-1.0*moms_other[16])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*moms_other[19]+1.7320508075688772*(moms_other[18]+moms_other[17])+moms_other[16]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*vtsq_other[3]-1.7320508075688772*(vtsq_other[2]+vtsq_other[1])+vtsq_other[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*vtsq_other[3]+1.7320508075688772*vtsq_other[2]-1.7320508075688772*vtsq_other[1]-1.0*vtsq_other[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*vtsq_other[3]-1.7320508075688772*vtsq_other[2]+1.7320508075688772*vtsq_other[1]-1.0*vtsq_other[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*vtsq_other[3]+1.7320508075688772*(vtsq_other[2]+vtsq_other[1])+vtsq_other[0]) < 0)) notCellAvg = false; 
 
  if (notCellAvg && (0.5*(3.0*alphaE[3]-1.7320508075688772*(alphaE[2]+alphaE[1])+alphaE[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*alphaE[3]+1.7320508075688772*alphaE[2]-1.7320508075688772*alphaE[1]-1.0*alphaE[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*alphaE[3]-1.7320508075688772*alphaE[2]+1.7320508075688772*alphaE[1]-1.0*alphaE[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*alphaE[3]+1.7320508075688772*(alphaE[2]+alphaE[1])+alphaE[0]) < 0)) notCellAvg = false; 
  if (notCellAvg) { 
    m0r[0] = moms_self[0]; 
    m0r[1] = moms_self[1]; 
    m0r[2] = moms_self[2]; 
    m0r[3] = moms_self[3]; 
    m1r[0] = moms_self[4]; 
    m1r[1] = moms_self[5]; 
    m1r[2] = moms_self[6]; 
    m1r[3] = moms_self[7]; 
    m1r[4] = moms_self[8]; 
    m1r[5] = moms_self[9]; 
    m1r[6] = moms_self[10]; 
    m1r[7] = moms_self[11]; 
    m1r[8] = moms_self[12]; 
    m1r[9] = moms_self[13]; 
    m1r[10] = moms_self[14]; 
    m1r[11] = moms_self[15]; 
    m2r[0] = moms_self[16]; 
    m2r[1] = moms_self[17]; 
    m2r[2] = moms_self[18]; 
    m2r[3] = moms_self[19]; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = boundary_corrections[5]; 
    cMr[6] = boundary_corrections[6]; 
    cMr[7] = boundary_corrections[7]; 
    cMr[8] = boundary_corrections[8]; 
    cMr[9] = boundary_corrections[9]; 
    cMr[10] = boundary_corrections[10]; 
    cMr[11] = boundary_corrections[11]; 
    cEr[0] = boundary_corrections[12]; 
    cEr[1] = boundary_corrections[13]; 
    cEr[2] = boundary_corrections[14]; 
    cEr[3] = boundary_corrections[15]; 
    u_selfr[0] = u_self[0]; 
    u_selfr[1] = u_self[1]; 
    u_selfr[2] = u_self[2]; 
    u_selfr[3] = u_self[3]; 
    u_selfr[4] = u_self[4]; 
    u_selfr[5] = u_self[5]; 
    u_selfr[6] = u_self[6]; 
    u_selfr[7] = u_self[7]; 
    u_selfr[8] = u_self[8]; 
    u_selfr[9] = u_self[9]; 
    u_selfr[10] = u_self[10]; 
    u_selfr[11] = u_self[11]; 
    u_otherr[0] = u_other[0]; 
    u_otherr[1] = u_other[1]; 
    u_otherr[2] = u_other[2]; 
    u_otherr[3] = u_other[3]; 
    u_otherr[4] = u_other[4]; 
    u_otherr[5] = u_other[5]; 
    u_otherr[6] = u_other[6]; 
    u_otherr[7] = u_other[7]; 
    u_otherr[8] = u_other[8]; 
    u_otherr[9] = u_other[9]; 
    u_otherr[10] = u_other[10]; 
    u_otherr[11] = u_other[11]; 
  } else { 
    m0r[0] = moms_self[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m1r[0] = moms_self[4]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    cMr[3] = 0.0; 
    m1r[4] = moms_self[8]; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = 0.0; 
    cMr[6] = 0.0; 
    cMr[7] = 0.0; 
    m1r[8] = moms_self[12]; 
    m1r[9] = 0.0; 
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    cMr[8] = boundary_corrections[8]; 
    cMr[9] = 0.0; 
    cMr[10] = 0.0; 
    cMr[11] = 0.0; 
    m2r[0] = moms_self[16]; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    cEr[0] = boundary_corrections[12]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    u_selfr[0] = u_self[0]; 
    u_selfr[1] = 0.0; 
    u_selfr[2] = 0.0; 
    u_selfr[3] = 0.0; 
    u_selfr[4] = 0.0; 
    u_selfr[5] = 0.0; 
    u_selfr[6] = 0.0; 
    u_selfr[7] = 0.0; 
    u_selfr[8] = 0.0; 
    u_selfr[9] = 0.0; 
    u_selfr[10] = 0.0; 
    u_selfr[11] = 0.0; 
    u_otherr[0] = u_other[0]; 
    u_otherr[1] = 0.0; 
    u_otherr[2] = 0.0; 
    u_otherr[3] = 0.0; 
    u_otherr[4] = 0.0; 
    u_otherr[5] = 0.0; 
    u_otherr[6] = 0.0; 
    u_otherr[7] = 0.0; 
    u_otherr[8] = 0.0; 
    u_otherr[9] = 0.0; 
    u_otherr[10] = 0.0; 
    u_otherr[11] = 0.0; 
  } 
 
  } else { 
  
    m0r[0] = 1.0; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m1r[0] = 1.0; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    cMr[0] = 0.0; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    cMr[3] = 0.0; 
    m1r[4] = 1.0; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    cMr[6] = 0.0; 
    cMr[7] = 0.0; 
    m1r[8] = 1.0; 
    m1r[9] = 0.0; 
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    cMr[8] = 0.0; 
    cMr[9] = 0.0; 
    cMr[10] = 0.0; 
    cMr[11] = 0.0; 
    m2r[0] = 1.0; 
    m2r[1] = 0.0; 
    m2r[2] = 0.0; 
    m2r[3] = 0.0; 
    cEr[0] = 0.0; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    u_selfr[0] = 1.0; 
    u_selfr[1] = 0.0; 
    u_selfr[2] = 0.0; 
    u_selfr[3] = 0.0; 
    u_selfr[4] = 0.0; 
    u_selfr[5] = 0.0; 
    u_selfr[6] = 0.0; 
    u_selfr[7] = 0.0; 
    u_selfr[8] = 0.0; 
    u_selfr[9] = 0.0; 
    u_selfr[10] = 0.0; 
    u_selfr[11] = 0.0; 
    u_otherr[0] = 1.0; 
    u_otherr[1] = 0.0; 
    u_otherr[2] = 0.0; 
    u_otherr[3] = 0.0; 
    u_otherr[4] = 0.0; 
    u_otherr[5] = 0.0; 
    u_otherr[6] = 0.0; 
    u_otherr[7] = 0.0; 
    u_otherr[8] = 0.0; 
    u_otherr[9] = 0.0; 
    u_otherr[10] = 0.0; 
    u_otherr[11] = 0.0; 
  
  }
  
  double momRHS[12] = {0.0}; 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossX ... // 
  gkyl_mat_set(A,0,0,m0r[0]); 
  gkyl_mat_set(A,0,1,m0r[1]); 
  gkyl_mat_set(A,0,2,m0r[2]); 
  gkyl_mat_set(A,0,3,m0r[3]); 
  gkyl_mat_set(A,1,0,m0r[1]); 
  gkyl_mat_set(A,1,1,m0r[0]); 
  gkyl_mat_set(A,1,2,m0r[3]); 
  gkyl_mat_set(A,1,3,m0r[2]); 
  gkyl_mat_set(A,2,0,m0r[2]); 
  gkyl_mat_set(A,2,1,m0r[3]); 
  gkyl_mat_set(A,2,2,m0r[0]); 
  gkyl_mat_set(A,2,3,m0r[1]); 
  gkyl_mat_set(A,3,0,m0r[3]); 
  gkyl_mat_set(A,3,1,m0r[2]); 
  gkyl_mat_set(A,3,2,m0r[1]); 
  gkyl_mat_set(A,3,3,m0r[0]); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,0,12,-(1.0*cMr[0])); 
  gkyl_mat_set(A,0,13,-(1.0*cMr[1])); 
  gkyl_mat_set(A,0,14,-(1.0*cMr[2])); 
  gkyl_mat_set(A,0,15,-(1.0*cMr[3])); 
  gkyl_mat_set(A,1,12,-(1.0*cMr[1])); 
  gkyl_mat_set(A,1,13,-(1.0*cMr[0])); 
  gkyl_mat_set(A,1,14,-(1.0*cMr[3])); 
  gkyl_mat_set(A,1,15,-(1.0*cMr[2])); 
  gkyl_mat_set(A,2,12,-(1.0*cMr[2])); 
  gkyl_mat_set(A,2,13,-(1.0*cMr[3])); 
  gkyl_mat_set(A,2,14,-(1.0*cMr[0])); 
  gkyl_mat_set(A,2,15,-(1.0*cMr[1])); 
  gkyl_mat_set(A,3,12,-(1.0*cMr[3])); 
  gkyl_mat_set(A,3,13,-(1.0*cMr[2])); 
  gkyl_mat_set(A,3,14,-(1.0*cMr[1])); 
  gkyl_mat_set(A,3,15,-(1.0*cMr[0])); 
 
  // ... Block from weak multiply of m_self, nu, m1X and uCrossX ... // 
  gkyl_mat_set(A,12,0,-(0.25*m0r[3]*u_selfr[3])-0.25*m0r[3]*u_otherr[3]-0.25*m0r[2]*u_selfr[2]-0.25*m0r[2]*u_otherr[2]-0.25*m0r[1]*u_selfr[1]-0.25*m0r[1]*u_otherr[1]-0.25*m0r[0]*u_selfr[0]-0.25*m0r[0]*u_otherr[0]+m1r[0]); 
  gkyl_mat_set(A,12,1,-(0.25*m0r[2]*u_selfr[3])-0.25*m0r[2]*u_otherr[3]-0.25*u_selfr[2]*m0r[3]-0.25*u_otherr[2]*m0r[3]-0.25*m0r[0]*u_selfr[1]-0.25*m0r[0]*u_otherr[1]+m1r[1]-0.25*u_selfr[0]*m0r[1]-0.25*u_otherr[0]*m0r[1]); 
  gkyl_mat_set(A,12,2,-(0.25*m0r[1]*u_selfr[3])-0.25*m0r[1]*u_otherr[3]-0.25*u_selfr[1]*m0r[3]-0.25*u_otherr[1]*m0r[3]-0.25*m0r[0]*u_selfr[2]-0.25*m0r[0]*u_otherr[2]+m1r[2]-0.25*u_selfr[0]*m0r[2]-0.25*u_otherr[0]*m0r[2]); 
  gkyl_mat_set(A,12,3,-(0.25*m0r[0]*u_selfr[3])-0.25*m0r[0]*u_otherr[3]+m1r[3]-0.25*u_selfr[0]*m0r[3]-0.25*u_otherr[0]*m0r[3]-0.25*m0r[1]*u_selfr[2]-0.25*m0r[1]*u_otherr[2]-0.25*u_selfr[1]*m0r[2]-0.25*u_otherr[1]*m0r[2]); 
  gkyl_mat_set(A,13,0,-(0.25*m0r[2]*u_selfr[3])-0.25*m0r[2]*u_otherr[3]-0.25*u_selfr[2]*m0r[3]-0.25*u_otherr[2]*m0r[3]-0.25*m0r[0]*u_selfr[1]-0.25*m0r[0]*u_otherr[1]+m1r[1]-0.25*u_selfr[0]*m0r[1]-0.25*u_otherr[0]*m0r[1]); 
  gkyl_mat_set(A,13,1,-(0.45*m0r[3]*u_selfr[3])-0.45*m0r[3]*u_otherr[3]-0.25*m0r[2]*u_selfr[2]-0.25*m0r[2]*u_otherr[2]-0.45*m0r[1]*u_selfr[1]-0.45*m0r[1]*u_otherr[1]-0.25*m0r[0]*u_selfr[0]-0.25*m0r[0]*u_otherr[0]+m1r[0]); 
  gkyl_mat_set(A,13,2,-(0.25*m0r[0]*u_selfr[3])-0.25*m0r[0]*u_otherr[3]+m1r[3]-0.25*u_selfr[0]*m0r[3]-0.25*u_otherr[0]*m0r[3]-0.25*m0r[1]*u_selfr[2]-0.25*m0r[1]*u_otherr[2]-0.25*u_selfr[1]*m0r[2]-0.25*u_otherr[1]*m0r[2]); 
  gkyl_mat_set(A,13,3,-(0.45*m0r[1]*u_selfr[3])-0.45*m0r[1]*u_otherr[3]-0.45*u_selfr[1]*m0r[3]-0.45*u_otherr[1]*m0r[3]-0.25*m0r[0]*u_selfr[2]-0.25*m0r[0]*u_otherr[2]+m1r[2]-0.25*u_selfr[0]*m0r[2]-0.25*u_otherr[0]*m0r[2]); 
  gkyl_mat_set(A,14,0,-(0.25*m0r[1]*u_selfr[3])-0.25*m0r[1]*u_otherr[3]-0.25*u_selfr[1]*m0r[3]-0.25*u_otherr[1]*m0r[3]-0.25*m0r[0]*u_selfr[2]-0.25*m0r[0]*u_otherr[2]+m1r[2]-0.25*u_selfr[0]*m0r[2]-0.25*u_otherr[0]*m0r[2]); 
  gkyl_mat_set(A,14,1,-(0.25*m0r[0]*u_selfr[3])-0.25*m0r[0]*u_otherr[3]+m1r[3]-0.25*u_selfr[0]*m0r[3]-0.25*u_otherr[0]*m0r[3]-0.25*m0r[1]*u_selfr[2]-0.25*m0r[1]*u_otherr[2]-0.25*u_selfr[1]*m0r[2]-0.25*u_otherr[1]*m0r[2]); 
  gkyl_mat_set(A,14,2,-(0.45*m0r[3]*u_selfr[3])-0.45*m0r[3]*u_otherr[3]-0.45*m0r[2]*u_selfr[2]-0.45*m0r[2]*u_otherr[2]-0.25*m0r[1]*u_selfr[1]-0.25*m0r[1]*u_otherr[1]-0.25*m0r[0]*u_selfr[0]-0.25*m0r[0]*u_otherr[0]+m1r[0]); 
  gkyl_mat_set(A,14,3,-(0.45*m0r[2]*u_selfr[3])-0.45*m0r[2]*u_otherr[3]-0.45*u_selfr[2]*m0r[3]-0.45*u_otherr[2]*m0r[3]-0.25*m0r[0]*u_selfr[1]-0.25*m0r[0]*u_otherr[1]+m1r[1]-0.25*u_selfr[0]*m0r[1]-0.25*u_otherr[0]*m0r[1]); 
  gkyl_mat_set(A,15,0,-(0.25*m0r[0]*u_selfr[3])-0.25*m0r[0]*u_otherr[3]+m1r[3]-0.25*u_selfr[0]*m0r[3]-0.25*u_otherr[0]*m0r[3]-0.25*m0r[1]*u_selfr[2]-0.25*m0r[1]*u_otherr[2]-0.25*u_selfr[1]*m0r[2]-0.25*u_otherr[1]*m0r[2]); 
  gkyl_mat_set(A,15,1,-(0.45*m0r[1]*u_selfr[3])-0.45*m0r[1]*u_otherr[3]-0.45*u_selfr[1]*m0r[3]-0.45*u_otherr[1]*m0r[3]-0.25*m0r[0]*u_selfr[2]-0.25*m0r[0]*u_otherr[2]+m1r[2]-0.25*u_selfr[0]*m0r[2]-0.25*u_otherr[0]*m0r[2]); 
  gkyl_mat_set(A,15,2,-(0.45*m0r[2]*u_selfr[3])-0.45*m0r[2]*u_otherr[3]-0.45*u_selfr[2]*m0r[3]-0.45*u_otherr[2]*m0r[3]-0.25*m0r[0]*u_selfr[1]-0.25*m0r[0]*u_otherr[1]+m1r[1]-0.25*u_selfr[0]*m0r[1]-0.25*u_otherr[0]*m0r[1]); 
  gkyl_mat_set(A,15,3,-(0.81*m0r[3]*u_selfr[3])-0.81*m0r[3]*u_otherr[3]-0.45*m0r[2]*u_selfr[2]-0.45*m0r[2]*u_otherr[2]-0.45*m0r[1]*u_selfr[1]-0.45*m0r[1]*u_otherr[1]-0.25*m0r[0]*u_selfr[0]-0.25*m0r[0]*u_otherr[0]+m1r[0]); 
 
  momRHS[0] += -(0.5*((alphaE[3]*u_selfr[3]-1.0*alphaE[3]*u_otherr[3]+alphaE[2]*u_selfr[2]-1.0*alphaE[2]*u_otherr[2]+alphaE[1]*u_selfr[1]-1.0*alphaE[1]*u_otherr[1]+alphaE[0]*u_selfr[0]-1.0*alphaE[0]*u_otherr[0])*m_sumDms-4.0*m1r[0])); 
  momRHS[1] += -(0.5*((alphaE[2]*u_selfr[3]-1.0*alphaE[2]*u_otherr[3]+(u_selfr[2]-1.0*u_otherr[2])*alphaE[3]+alphaE[0]*u_selfr[1]-1.0*alphaE[0]*u_otherr[1]+(u_selfr[0]-1.0*u_otherr[0])*alphaE[1])*m_sumDms-4.0*m1r[1])); 
  momRHS[2] += -(0.5*((alphaE[1]*u_selfr[3]-1.0*alphaE[1]*u_otherr[3]+(u_selfr[1]-1.0*u_otherr[1])*alphaE[3]+alphaE[0]*u_selfr[2]-1.0*alphaE[0]*u_otherr[2]+(u_selfr[0]-1.0*u_otherr[0])*alphaE[2])*m_sumDms-4.0*m1r[2])); 
  momRHS[3] += -(0.5*((alphaE[0]*u_selfr[3]-1.0*alphaE[0]*u_otherr[3]+(u_selfr[0]-1.0*u_otherr[0])*alphaE[3]+alphaE[1]*u_selfr[2]-1.0*alphaE[1]*u_otherr[2]+(u_selfr[1]-1.0*u_otherr[1])*alphaE[2])*m_sumDms-4.0*m1r[3])); 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossY ... // 
  gkyl_mat_set(A,4,4,m0r[0]); 
  gkyl_mat_set(A,4,5,m0r[1]); 
  gkyl_mat_set(A,4,6,m0r[2]); 
  gkyl_mat_set(A,4,7,m0r[3]); 
  gkyl_mat_set(A,5,4,m0r[1]); 
  gkyl_mat_set(A,5,5,m0r[0]); 
  gkyl_mat_set(A,5,6,m0r[3]); 
  gkyl_mat_set(A,5,7,m0r[2]); 
  gkyl_mat_set(A,6,4,m0r[2]); 
  gkyl_mat_set(A,6,5,m0r[3]); 
  gkyl_mat_set(A,6,6,m0r[0]); 
  gkyl_mat_set(A,6,7,m0r[1]); 
  gkyl_mat_set(A,7,4,m0r[3]); 
  gkyl_mat_set(A,7,5,m0r[2]); 
  gkyl_mat_set(A,7,6,m0r[1]); 
  gkyl_mat_set(A,7,7,m0r[0]); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,4,12,-(1.0*cMr[4])); 
  gkyl_mat_set(A,4,13,-(1.0*cMr[5])); 
  gkyl_mat_set(A,4,14,-(1.0*cMr[6])); 
  gkyl_mat_set(A,4,15,-(1.0*cMr[7])); 
  gkyl_mat_set(A,5,12,-(1.0*cMr[5])); 
  gkyl_mat_set(A,5,13,-(1.0*cMr[4])); 
  gkyl_mat_set(A,5,14,-(1.0*cMr[7])); 
  gkyl_mat_set(A,5,15,-(1.0*cMr[6])); 
  gkyl_mat_set(A,6,12,-(1.0*cMr[6])); 
  gkyl_mat_set(A,6,13,-(1.0*cMr[7])); 
  gkyl_mat_set(A,6,14,-(1.0*cMr[4])); 
  gkyl_mat_set(A,6,15,-(1.0*cMr[5])); 
  gkyl_mat_set(A,7,12,-(1.0*cMr[7])); 
  gkyl_mat_set(A,7,13,-(1.0*cMr[6])); 
  gkyl_mat_set(A,7,14,-(1.0*cMr[5])); 
  gkyl_mat_set(A,7,15,-(1.0*cMr[4])); 
 
  // ... Block from weak multiply of m_self, nu, m1Y and uCrossY ... // 
  gkyl_mat_set(A,12,4,-(0.25*m0r[3]*u_selfr[7])-0.25*m0r[3]*u_otherr[7]-0.25*m0r[2]*u_selfr[6]-0.25*m0r[2]*u_otherr[6]-0.25*m0r[1]*u_selfr[5]-0.25*m0r[1]*u_otherr[5]-0.25*m0r[0]*u_selfr[4]-0.25*m0r[0]*u_otherr[4]+m1r[4]); 
  gkyl_mat_set(A,12,5,-(0.25*m0r[2]*u_selfr[7])-0.25*m0r[2]*u_otherr[7]-0.25*m0r[3]*u_selfr[6]-0.25*m0r[3]*u_otherr[6]-0.25*m0r[0]*u_selfr[5]-0.25*m0r[0]*u_otherr[5]+m1r[5]-0.25*m0r[1]*u_selfr[4]-0.25*m0r[1]*u_otherr[4]); 
  gkyl_mat_set(A,12,6,-(0.25*m0r[1]*u_selfr[7])-0.25*m0r[1]*u_otherr[7]-0.25*m0r[0]*u_selfr[6]-0.25*m0r[0]*u_otherr[6]+m1r[6]-0.25*m0r[3]*u_selfr[5]-0.25*m0r[3]*u_otherr[5]-0.25*m0r[2]*u_selfr[4]-0.25*m0r[2]*u_otherr[4]); 
  gkyl_mat_set(A,12,7,-(0.25*m0r[0]*u_selfr[7])-0.25*m0r[0]*u_otherr[7]+m1r[7]-0.25*m0r[1]*u_selfr[6]-0.25*m0r[1]*u_otherr[6]-0.25*m0r[2]*u_selfr[5]-0.25*m0r[2]*u_otherr[5]-0.25*m0r[3]*u_selfr[4]-0.25*m0r[3]*u_otherr[4]); 
  gkyl_mat_set(A,13,4,-(0.25*m0r[2]*u_selfr[7])-0.25*m0r[2]*u_otherr[7]-0.25*m0r[3]*u_selfr[6]-0.25*m0r[3]*u_otherr[6]-0.25*m0r[0]*u_selfr[5]-0.25*m0r[0]*u_otherr[5]+m1r[5]-0.25*m0r[1]*u_selfr[4]-0.25*m0r[1]*u_otherr[4]); 
  gkyl_mat_set(A,13,5,-(0.45*m0r[3]*u_selfr[7])-0.45*m0r[3]*u_otherr[7]-0.25*m0r[2]*u_selfr[6]-0.25*m0r[2]*u_otherr[6]-0.45*m0r[1]*u_selfr[5]-0.45*m0r[1]*u_otherr[5]-0.25*m0r[0]*u_selfr[4]-0.25*m0r[0]*u_otherr[4]+m1r[4]); 
  gkyl_mat_set(A,13,6,-(0.25*m0r[0]*u_selfr[7])-0.25*m0r[0]*u_otherr[7]+m1r[7]-0.25*m0r[1]*u_selfr[6]-0.25*m0r[1]*u_otherr[6]-0.25*m0r[2]*u_selfr[5]-0.25*m0r[2]*u_otherr[5]-0.25*m0r[3]*u_selfr[4]-0.25*m0r[3]*u_otherr[4]); 
  gkyl_mat_set(A,13,7,-(0.45*m0r[1]*u_selfr[7])-0.45*m0r[1]*u_otherr[7]-0.25*m0r[0]*u_selfr[6]-0.25*m0r[0]*u_otherr[6]+m1r[6]-0.45*m0r[3]*u_selfr[5]-0.45*m0r[3]*u_otherr[5]-0.25*m0r[2]*u_selfr[4]-0.25*m0r[2]*u_otherr[4]); 
  gkyl_mat_set(A,14,4,-(0.25*m0r[1]*u_selfr[7])-0.25*m0r[1]*u_otherr[7]-0.25*m0r[0]*u_selfr[6]-0.25*m0r[0]*u_otherr[6]+m1r[6]-0.25*m0r[3]*u_selfr[5]-0.25*m0r[3]*u_otherr[5]-0.25*m0r[2]*u_selfr[4]-0.25*m0r[2]*u_otherr[4]); 
  gkyl_mat_set(A,14,5,-(0.25*m0r[0]*u_selfr[7])-0.25*m0r[0]*u_otherr[7]+m1r[7]-0.25*m0r[1]*u_selfr[6]-0.25*m0r[1]*u_otherr[6]-0.25*m0r[2]*u_selfr[5]-0.25*m0r[2]*u_otherr[5]-0.25*m0r[3]*u_selfr[4]-0.25*m0r[3]*u_otherr[4]); 
  gkyl_mat_set(A,14,6,-(0.45*m0r[3]*u_selfr[7])-0.45*m0r[3]*u_otherr[7]-0.45*m0r[2]*u_selfr[6]-0.45*m0r[2]*u_otherr[6]-0.25*m0r[1]*u_selfr[5]-0.25*m0r[1]*u_otherr[5]-0.25*m0r[0]*u_selfr[4]-0.25*m0r[0]*u_otherr[4]+m1r[4]); 
  gkyl_mat_set(A,14,7,-(0.45*m0r[2]*u_selfr[7])-0.45*m0r[2]*u_otherr[7]-0.45*m0r[3]*u_selfr[6]-0.45*m0r[3]*u_otherr[6]-0.25*m0r[0]*u_selfr[5]-0.25*m0r[0]*u_otherr[5]+m1r[5]-0.25*m0r[1]*u_selfr[4]-0.25*m0r[1]*u_otherr[4]); 
  gkyl_mat_set(A,15,4,-(0.25*m0r[0]*u_selfr[7])-0.25*m0r[0]*u_otherr[7]+m1r[7]-0.25*m0r[1]*u_selfr[6]-0.25*m0r[1]*u_otherr[6]-0.25*m0r[2]*u_selfr[5]-0.25*m0r[2]*u_otherr[5]-0.25*m0r[3]*u_selfr[4]-0.25*m0r[3]*u_otherr[4]); 
  gkyl_mat_set(A,15,5,-(0.45*m0r[1]*u_selfr[7])-0.45*m0r[1]*u_otherr[7]-0.25*m0r[0]*u_selfr[6]-0.25*m0r[0]*u_otherr[6]+m1r[6]-0.45*m0r[3]*u_selfr[5]-0.45*m0r[3]*u_otherr[5]-0.25*m0r[2]*u_selfr[4]-0.25*m0r[2]*u_otherr[4]); 
  gkyl_mat_set(A,15,6,-(0.45*m0r[2]*u_selfr[7])-0.45*m0r[2]*u_otherr[7]-0.45*m0r[3]*u_selfr[6]-0.45*m0r[3]*u_otherr[6]-0.25*m0r[0]*u_selfr[5]-0.25*m0r[0]*u_otherr[5]+m1r[5]-0.25*m0r[1]*u_selfr[4]-0.25*m0r[1]*u_otherr[4]); 
  gkyl_mat_set(A,15,7,-(0.81*m0r[3]*u_selfr[7])-0.81*m0r[3]*u_otherr[7]-0.45*m0r[2]*u_selfr[6]-0.45*m0r[2]*u_otherr[6]-0.45*m0r[1]*u_selfr[5]-0.45*m0r[1]*u_otherr[5]-0.25*m0r[0]*u_selfr[4]-0.25*m0r[0]*u_otherr[4]+m1r[4]); 
 
  momRHS[4] += -(0.5*((alphaE[3]*u_selfr[7]-1.0*alphaE[3]*u_otherr[7]+alphaE[2]*u_selfr[6]-1.0*alphaE[2]*u_otherr[6]+alphaE[1]*u_selfr[5]-1.0*alphaE[1]*u_otherr[5]+alphaE[0]*u_selfr[4]-1.0*alphaE[0]*u_otherr[4])*m_sumDms-4.0*m1r[4])); 
  momRHS[5] += -(0.5*((alphaE[2]*u_selfr[7]-1.0*alphaE[2]*u_otherr[7]+alphaE[3]*u_selfr[6]-1.0*alphaE[3]*u_otherr[6]+alphaE[0]*u_selfr[5]-1.0*alphaE[0]*u_otherr[5]+alphaE[1]*u_selfr[4]-1.0*alphaE[1]*u_otherr[4])*m_sumDms-4.0*m1r[5])); 
  momRHS[6] += -(0.5*((alphaE[1]*u_selfr[7]-1.0*alphaE[1]*u_otherr[7]+alphaE[0]*u_selfr[6]-1.0*alphaE[0]*u_otherr[6]+alphaE[3]*u_selfr[5]-1.0*alphaE[3]*u_otherr[5]+alphaE[2]*u_selfr[4]-1.0*alphaE[2]*u_otherr[4])*m_sumDms-4.0*m1r[6])); 
  momRHS[7] += -(0.5*((alphaE[0]*u_selfr[7]-1.0*alphaE[0]*u_otherr[7]+alphaE[1]*u_selfr[6]-1.0*alphaE[1]*u_otherr[6]+alphaE[2]*u_selfr[5]-1.0*alphaE[2]*u_otherr[5]+alphaE[3]*u_selfr[4]-1.0*alphaE[3]*u_otherr[4])*m_sumDms-4.0*m1r[7])); 
 
  // ... Block from weak multiply of m_self, nu, M0 and uCrossZ ... // 
  gkyl_mat_set(A,8,8,m0r[0]); 
  gkyl_mat_set(A,8,9,m0r[1]); 
  gkyl_mat_set(A,8,10,m0r[2]); 
  gkyl_mat_set(A,8,11,m0r[3]); 
  gkyl_mat_set(A,9,8,m0r[1]); 
  gkyl_mat_set(A,9,9,m0r[0]); 
  gkyl_mat_set(A,9,10,m0r[3]); 
  gkyl_mat_set(A,9,11,m0r[2]); 
  gkyl_mat_set(A,10,8,m0r[2]); 
  gkyl_mat_set(A,10,9,m0r[3]); 
  gkyl_mat_set(A,10,10,m0r[0]); 
  gkyl_mat_set(A,10,11,m0r[1]); 
  gkyl_mat_set(A,11,8,m0r[3]); 
  gkyl_mat_set(A,11,9,m0r[2]); 
  gkyl_mat_set(A,11,10,m0r[1]); 
  gkyl_mat_set(A,11,11,m0r[0]); 
 
  // ... Block from correction to momentum conservation (self) ... // 
  gkyl_mat_set(A,8,12,-(1.0*cMr[8])); 
  gkyl_mat_set(A,8,13,-(1.0*cMr[9])); 
  gkyl_mat_set(A,8,14,-(1.0*cMr[10])); 
  gkyl_mat_set(A,8,15,-(1.0*cMr[11])); 
  gkyl_mat_set(A,9,12,-(1.0*cMr[9])); 
  gkyl_mat_set(A,9,13,-(1.0*cMr[8])); 
  gkyl_mat_set(A,9,14,-(1.0*cMr[11])); 
  gkyl_mat_set(A,9,15,-(1.0*cMr[10])); 
  gkyl_mat_set(A,10,12,-(1.0*cMr[10])); 
  gkyl_mat_set(A,10,13,-(1.0*cMr[11])); 
  gkyl_mat_set(A,10,14,-(1.0*cMr[8])); 
  gkyl_mat_set(A,10,15,-(1.0*cMr[9])); 
  gkyl_mat_set(A,11,12,-(1.0*cMr[11])); 
  gkyl_mat_set(A,11,13,-(1.0*cMr[10])); 
  gkyl_mat_set(A,11,14,-(1.0*cMr[9])); 
  gkyl_mat_set(A,11,15,-(1.0*cMr[8])); 
 
  // ... Block from weak multiply of m_self, nu, m1Z and uCrossZ ... // 
  gkyl_mat_set(A,12,8,-(0.25*m0r[3]*u_selfr[11])-0.25*m0r[3]*u_otherr[11]-0.25*m0r[2]*u_selfr[10]-0.25*m0r[2]*u_otherr[10]-0.25*m0r[1]*u_selfr[9]-0.25*m0r[1]*u_otherr[9]-0.25*m0r[0]*u_selfr[8]-0.25*m0r[0]*u_otherr[8]+m1r[8]); 
  gkyl_mat_set(A,12,9,-(0.25*m0r[2]*u_selfr[11])-0.25*m0r[2]*u_otherr[11]-0.25*m0r[3]*u_selfr[10]-0.25*m0r[3]*u_otherr[10]-0.25*m0r[0]*u_selfr[9]-0.25*m0r[0]*u_otherr[9]+m1r[9]-0.25*m0r[1]*u_selfr[8]-0.25*m0r[1]*u_otherr[8]); 
  gkyl_mat_set(A,12,10,-(0.25*m0r[1]*u_selfr[11])-0.25*m0r[1]*u_otherr[11]-0.25*m0r[0]*u_selfr[10]-0.25*m0r[0]*u_otherr[10]+m1r[10]-0.25*m0r[3]*u_selfr[9]-0.25*m0r[3]*u_otherr[9]-0.25*m0r[2]*u_selfr[8]-0.25*m0r[2]*u_otherr[8]); 
  gkyl_mat_set(A,12,11,-(0.25*m0r[0]*u_selfr[11])-0.25*m0r[0]*u_otherr[11]+m1r[11]-0.25*m0r[1]*u_selfr[10]-0.25*m0r[1]*u_otherr[10]-0.25*m0r[2]*u_selfr[9]-0.25*m0r[2]*u_otherr[9]-0.25*m0r[3]*u_selfr[8]-0.25*m0r[3]*u_otherr[8]); 
  gkyl_mat_set(A,13,8,-(0.25*m0r[2]*u_selfr[11])-0.25*m0r[2]*u_otherr[11]-0.25*m0r[3]*u_selfr[10]-0.25*m0r[3]*u_otherr[10]-0.25*m0r[0]*u_selfr[9]-0.25*m0r[0]*u_otherr[9]+m1r[9]-0.25*m0r[1]*u_selfr[8]-0.25*m0r[1]*u_otherr[8]); 
  gkyl_mat_set(A,13,9,-(0.45*m0r[3]*u_selfr[11])-0.45*m0r[3]*u_otherr[11]-0.25*m0r[2]*u_selfr[10]-0.25*m0r[2]*u_otherr[10]-0.45*m0r[1]*u_selfr[9]-0.45*m0r[1]*u_otherr[9]-0.25*m0r[0]*u_selfr[8]-0.25*m0r[0]*u_otherr[8]+m1r[8]); 
  gkyl_mat_set(A,13,10,-(0.25*m0r[0]*u_selfr[11])-0.25*m0r[0]*u_otherr[11]+m1r[11]-0.25*m0r[1]*u_selfr[10]-0.25*m0r[1]*u_otherr[10]-0.25*m0r[2]*u_selfr[9]-0.25*m0r[2]*u_otherr[9]-0.25*m0r[3]*u_selfr[8]-0.25*m0r[3]*u_otherr[8]); 
  gkyl_mat_set(A,13,11,-(0.45*m0r[1]*u_selfr[11])-0.45*m0r[1]*u_otherr[11]-0.25*m0r[0]*u_selfr[10]-0.25*m0r[0]*u_otherr[10]+m1r[10]-0.45*m0r[3]*u_selfr[9]-0.45*m0r[3]*u_otherr[9]-0.25*m0r[2]*u_selfr[8]-0.25*m0r[2]*u_otherr[8]); 
  gkyl_mat_set(A,14,8,-(0.25*m0r[1]*u_selfr[11])-0.25*m0r[1]*u_otherr[11]-0.25*m0r[0]*u_selfr[10]-0.25*m0r[0]*u_otherr[10]+m1r[10]-0.25*m0r[3]*u_selfr[9]-0.25*m0r[3]*u_otherr[9]-0.25*m0r[2]*u_selfr[8]-0.25*m0r[2]*u_otherr[8]); 
  gkyl_mat_set(A,14,9,-(0.25*m0r[0]*u_selfr[11])-0.25*m0r[0]*u_otherr[11]+m1r[11]-0.25*m0r[1]*u_selfr[10]-0.25*m0r[1]*u_otherr[10]-0.25*m0r[2]*u_selfr[9]-0.25*m0r[2]*u_otherr[9]-0.25*m0r[3]*u_selfr[8]-0.25*m0r[3]*u_otherr[8]); 
  gkyl_mat_set(A,14,10,-(0.45*m0r[3]*u_selfr[11])-0.45*m0r[3]*u_otherr[11]-0.45*m0r[2]*u_selfr[10]-0.45*m0r[2]*u_otherr[10]-0.25*m0r[1]*u_selfr[9]-0.25*m0r[1]*u_otherr[9]-0.25*m0r[0]*u_selfr[8]-0.25*m0r[0]*u_otherr[8]+m1r[8]); 
  gkyl_mat_set(A,14,11,-(0.45*m0r[2]*u_selfr[11])-0.45*m0r[2]*u_otherr[11]-0.45*m0r[3]*u_selfr[10]-0.45*m0r[3]*u_otherr[10]-0.25*m0r[0]*u_selfr[9]-0.25*m0r[0]*u_otherr[9]+m1r[9]-0.25*m0r[1]*u_selfr[8]-0.25*m0r[1]*u_otherr[8]); 
  gkyl_mat_set(A,15,8,-(0.25*m0r[0]*u_selfr[11])-0.25*m0r[0]*u_otherr[11]+m1r[11]-0.25*m0r[1]*u_selfr[10]-0.25*m0r[1]*u_otherr[10]-0.25*m0r[2]*u_selfr[9]-0.25*m0r[2]*u_otherr[9]-0.25*m0r[3]*u_selfr[8]-0.25*m0r[3]*u_otherr[8]); 
  gkyl_mat_set(A,15,9,-(0.45*m0r[1]*u_selfr[11])-0.45*m0r[1]*u_otherr[11]-0.25*m0r[0]*u_selfr[10]-0.25*m0r[0]*u_otherr[10]+m1r[10]-0.45*m0r[3]*u_selfr[9]-0.45*m0r[3]*u_otherr[9]-0.25*m0r[2]*u_selfr[8]-0.25*m0r[2]*u_otherr[8]); 
  gkyl_mat_set(A,15,10,-(0.45*m0r[2]*u_selfr[11])-0.45*m0r[2]*u_otherr[11]-0.45*m0r[3]*u_selfr[10]-0.45*m0r[3]*u_otherr[10]-0.25*m0r[0]*u_selfr[9]-0.25*m0r[0]*u_otherr[9]+m1r[9]-0.25*m0r[1]*u_selfr[8]-0.25*m0r[1]*u_otherr[8]); 
  gkyl_mat_set(A,15,11,-(0.81*m0r[3]*u_selfr[11])-0.81*m0r[3]*u_otherr[11]-0.45*m0r[2]*u_selfr[10]-0.45*m0r[2]*u_otherr[10]-0.45*m0r[1]*u_selfr[9]-0.45*m0r[1]*u_otherr[9]-0.25*m0r[0]*u_selfr[8]-0.25*m0r[0]*u_otherr[8]+m1r[8]); 
 
  momRHS[8] += -(0.5*((alphaE[3]*u_selfr[11]-1.0*alphaE[3]*u_otherr[11]+alphaE[2]*u_selfr[10]-1.0*alphaE[2]*u_otherr[10]+alphaE[1]*u_selfr[9]-1.0*alphaE[1]*u_otherr[9]+alphaE[0]*u_selfr[8]-1.0*alphaE[0]*u_otherr[8])*m_sumDms-4.0*m1r[8])); 
  momRHS[9] += -(0.5*((alphaE[2]*u_selfr[11]-1.0*alphaE[2]*u_otherr[11]+alphaE[3]*u_selfr[10]-1.0*alphaE[3]*u_otherr[10]+alphaE[0]*u_selfr[9]-1.0*alphaE[0]*u_otherr[9]+alphaE[1]*u_selfr[8]-1.0*alphaE[1]*u_otherr[8])*m_sumDms-4.0*m1r[9])); 
  momRHS[10] += -(0.5*((alphaE[1]*u_selfr[11]-1.0*alphaE[1]*u_otherr[11]+alphaE[0]*u_selfr[10]-1.0*alphaE[0]*u_otherr[10]+alphaE[3]*u_selfr[9]-1.0*alphaE[3]*u_otherr[9]+alphaE[2]*u_selfr[8]-1.0*alphaE[2]*u_otherr[8])*m_sumDms-4.0*m1r[10])); 
  momRHS[11] += -(0.5*((alphaE[0]*u_selfr[11]-1.0*alphaE[0]*u_otherr[11]+alphaE[1]*u_selfr[10]-1.0*alphaE[1]*u_otherr[10]+alphaE[2]*u_selfr[9]-1.0*alphaE[2]*u_otherr[9]+alphaE[3]*u_selfr[8]-1.0*alphaE[3]*u_otherr[8])*m_sumDms-4.0*m1r[11])); 
 
  double ucMSelf[4] = {0.0}; 
  double ucMOther[4] = {0.0}; 
  for (int vd=0; vd<3; vd++) 
  { 
    int a0 = 4*vd; 
    // Contribution to dot-product from weak multiplication of vd component. 
    ucMSelf[0] += 0.5*cMr[a0+3]*u_selfr[a0+3]+0.5*cMr[a0+2]*u_selfr[a0+2]+0.5*cMr[a0+1]*u_selfr[a0+1]+0.5*cMr[a0]*u_selfr[a0]; 
    ucMSelf[1] += 0.5*cMr[a0+2]*u_selfr[a0+3]+0.5*u_selfr[a0+2]*cMr[a0+3]+0.5*cMr[a0]*u_selfr[a0+1]+0.5*u_selfr[a0]*cMr[a0+1]; 
    ucMSelf[2] += 0.5*cMr[a0+1]*u_selfr[a0+3]+0.5*u_selfr[a0+1]*cMr[a0+3]+0.5*cMr[a0]*u_selfr[a0+2]+0.5*u_selfr[a0]*cMr[a0+2]; 
    ucMSelf[3] += 0.5*cMr[a0]*u_selfr[a0+3]+0.5*u_selfr[a0]*cMr[a0+3]+0.5*cMr[a0+1]*u_selfr[a0+2]+0.5*u_selfr[a0+1]*cMr[a0+2]; 
    ucMOther[0] += 0.5*cMr[a0+3]*u_otherr[a0+3]+0.5*cMr[a0+2]*u_otherr[a0+2]+0.5*cMr[a0+1]*u_otherr[a0+1]+0.5*cMr[a0]*u_otherr[a0]; 
    ucMOther[1] += 0.5*cMr[a0+2]*u_otherr[a0+3]+0.5*u_otherr[a0+2]*cMr[a0+3]+0.5*cMr[a0]*u_otherr[a0+1]+0.5*u_otherr[a0]*cMr[a0+1]; 
    ucMOther[2] += 0.5*cMr[a0+1]*u_otherr[a0+3]+0.5*u_otherr[a0+1]*cMr[a0+3]+0.5*cMr[a0]*u_otherr[a0+2]+0.5*u_otherr[a0]*cMr[a0+2]; 
    ucMOther[3] += 0.5*cMr[a0]*u_otherr[a0+3]+0.5*u_otherr[a0]*cMr[a0+3]+0.5*cMr[a0+1]*u_otherr[a0+2]+0.5*u_otherr[a0+1]*cMr[a0+2]; 
  } 
 
  // ... Block from correction to (self) 2nd moment of collision operator ... // 
  gkyl_mat_set(A,12,12,0.5*ucMSelf[0]+0.5*ucMOther[0]+3.0*m0r[0]-1.0*cEr[0]); 
  gkyl_mat_set(A,12,13,0.5*ucMSelf[1]+0.5*ucMOther[1]+3.0*m0r[1]-1.0*cEr[1]); 
  gkyl_mat_set(A,12,14,0.5*ucMSelf[2]+0.5*ucMOther[2]+3.0*m0r[2]-1.0*cEr[2]); 
  gkyl_mat_set(A,12,15,0.5*ucMSelf[3]+0.5*ucMOther[3]+3.0*m0r[3]-1.0*cEr[3]); 
  gkyl_mat_set(A,13,12,0.5*ucMSelf[1]+0.5*ucMOther[1]+3.0*m0r[1]-1.0*cEr[1]); 
  gkyl_mat_set(A,13,13,0.5*ucMSelf[0]+0.5*ucMOther[0]+3.0*m0r[0]-1.0*cEr[0]); 
  gkyl_mat_set(A,13,14,0.5*ucMSelf[3]+0.5*ucMOther[3]+3.0*m0r[3]-1.0*cEr[3]); 
  gkyl_mat_set(A,13,15,0.5*ucMSelf[2]+0.5*ucMOther[2]+3.0*m0r[2]-1.0*cEr[2]); 
  gkyl_mat_set(A,14,12,0.5*ucMSelf[2]+0.5*ucMOther[2]+3.0*m0r[2]-1.0*cEr[2]); 
  gkyl_mat_set(A,14,13,0.5*ucMSelf[3]+0.5*ucMOther[3]+3.0*m0r[3]-1.0*cEr[3]); 
  gkyl_mat_set(A,14,14,0.5*ucMSelf[0]+0.5*ucMOther[0]+3.0*m0r[0]-1.0*cEr[0]); 
  gkyl_mat_set(A,14,15,0.5*ucMSelf[1]+0.5*ucMOther[1]+3.0*m0r[1]-1.0*cEr[1]); 
  gkyl_mat_set(A,15,12,0.5*ucMSelf[3]+0.5*ucMOther[3]+3.0*m0r[3]-1.0*cEr[3]); 
  gkyl_mat_set(A,15,13,0.5*ucMSelf[2]+0.5*ucMOther[2]+3.0*m0r[2]-1.0*cEr[2]); 
  gkyl_mat_set(A,15,14,0.5*ucMSelf[1]+0.5*ucMOther[1]+3.0*m0r[1]-1.0*cEr[1]); 
  gkyl_mat_set(A,15,15,0.5*ucMSelf[0]+0.5*ucMOther[0]+3.0*m0r[0]-1.0*cEr[0]); 
 
  double uM1Self[4] = {0.0}; 
  double uM1Other[4] = {0.0}; 
  double uSumSq[4] = {0.0}; 
  for (int vd=0; vd<3; vd++) 
  { 
    int a0 = 4*vd; 
    // Dot product terms in energy equation RHS. 
    uM1Self[0] += 0.5*m1r[a0+3]*u_selfr[a0+3]+0.5*m1r[a0+2]*u_selfr[a0+2]+0.5*m1r[a0+1]*u_selfr[a0+1]+0.5*m1r[a0]*u_selfr[a0]; 
    uM1Self[1] += 0.5*m1r[a0+2]*u_selfr[a0+3]+0.5*u_selfr[a0+2]*m1r[a0+3]+0.5*m1r[a0]*u_selfr[a0+1]+0.5*u_selfr[a0]*m1r[a0+1]; 
    uM1Self[2] += 0.5*m1r[a0+1]*u_selfr[a0+3]+0.5*u_selfr[a0+1]*m1r[a0+3]+0.5*m1r[a0]*u_selfr[a0+2]+0.5*u_selfr[a0]*m1r[a0+2]; 
    uM1Self[3] += 0.5*m1r[a0]*u_selfr[a0+3]+0.5*u_selfr[a0]*m1r[a0+3]+0.5*m1r[a0+1]*u_selfr[a0+2]+0.5*u_selfr[a0+1]*m1r[a0+2]; 
    uM1Other[0] += 0.5*m1r[a0+3]*u_otherr[a0+3]+0.5*m1r[a0+2]*u_otherr[a0+2]+0.5*m1r[a0+1]*u_otherr[a0+1]+0.5*m1r[a0]*u_otherr[a0]; 
    uM1Other[1] += 0.5*m1r[a0+2]*u_otherr[a0+3]+0.5*u_otherr[a0+2]*m1r[a0+3]+0.5*m1r[a0]*u_otherr[a0+1]+0.5*u_otherr[a0]*m1r[a0+1]; 
    uM1Other[2] += 0.5*m1r[a0+1]*u_otherr[a0+3]+0.5*u_otherr[a0+1]*m1r[a0+3]+0.5*m1r[a0]*u_otherr[a0+2]+0.5*u_otherr[a0]*m1r[a0+2]; 
    uM1Other[3] += 0.5*m1r[a0]*u_otherr[a0+3]+0.5*u_otherr[a0]*m1r[a0+3]+0.5*m1r[a0+1]*u_otherr[a0+2]+0.5*u_otherr[a0+1]*m1r[a0+2]; 
  const double u_selfr0R2 = pow(u_selfr[a0],2);
  const double u_selfr1R2 = pow(u_selfr[a0+1],2);
  const double u_selfr2R2 = pow(u_selfr[a0+2],2);
  const double u_selfr3R2 = pow(u_selfr[a0+3],2);
  const double u_otherr0R2 = pow(u_otherr[a0],2);
  const double u_otherr1R2 = pow(u_otherr[a0+1],2);
  const double u_otherr2R2 = pow(u_otherr[a0+2],2);
  const double u_otherr3R2 = pow(u_otherr[a0+3],2);

  uSumSq[0] += 0.5*u_selfr3R2-1.0*u_otherr[a0+3]*u_selfr[a0+3]+0.5*u_otherr3R2+0.5*u_selfr2R2-1.0*u_otherr[a0+2]*u_selfr[a0+2]+0.5*u_otherr2R2+0.5*u_selfr1R2-1.0*u_otherr[a0+1]*u_selfr[a0+1]+0.5*u_otherr1R2+0.5*u_selfr0R2-1.0*u_otherr[a0]*u_selfr[a0]+0.5*u_otherr0R2; 
  uSumSq[1] += u_selfr[a0+2]*u_selfr[a0+3]-1.0*u_otherr[a0+2]*u_selfr[a0+3]-1.0*u_selfr[a0+2]*u_otherr[a0+3]+u_otherr[a0+2]*u_otherr[a0+3]+u_selfr[a0]*u_selfr[a0+1]-1.0*u_otherr[a0]*u_selfr[a0+1]-1.0*u_selfr[a0]*u_otherr[a0+1]+u_otherr[a0]*u_otherr[a0+1]; 
  uSumSq[2] += u_selfr[a0+1]*u_selfr[a0+3]-1.0*u_otherr[a0+1]*u_selfr[a0+3]-1.0*u_selfr[a0+1]*u_otherr[a0+3]+u_otherr[a0+1]*u_otherr[a0+3]+u_selfr[a0]*u_selfr[a0+2]-1.0*u_otherr[a0]*u_selfr[a0+2]-1.0*u_selfr[a0]*u_otherr[a0+2]+u_otherr[a0]*u_otherr[a0+2]; 
  uSumSq[3] += u_selfr[a0]*u_selfr[a0+3]-1.0*u_otherr[a0]*u_selfr[a0+3]-1.0*u_selfr[a0]*u_otherr[a0+3]+u_otherr[a0]*u_otherr[a0+3]+u_selfr[a0+1]*u_selfr[a0+2]-1.0*u_otherr[a0+1]*u_selfr[a0+2]-1.0*u_selfr[a0+1]*u_otherr[a0+2]+u_otherr[a0+1]*u_otherr[a0+2]; 
  } 
 
  double m_diff = m_other-m_self;
  double enRHS[4] = {0.0}; 
  enRHS[0] = (1.5*alphaE[3]*vtsq_other[3]*m_other)/m_self+(1.5*alphaE[2]*vtsq_other[2]*m_other)/m_self+(1.5*alphaE[1]*vtsq_other[1]*m_other)/m_self+(1.5*alphaE[0]*vtsq_other[0]*m_other)/m_self+(0.25*alphaE[3]*uSumSq[3]*m_diff)/m_self+(0.25*alphaE[2]*uSumSq[2]*m_diff)/m_self+(0.25*alphaE[1]*uSumSq[1]*m_diff)/m_self+(0.25*alphaE[0]*uSumSq[0]*m_diff)/m_self-1.5*alphaE[3]*vtsq_self[3]-1.5*alphaE[2]*vtsq_self[2]-1.5*alphaE[1]*vtsq_self[1]-1.5*alphaE[0]*vtsq_self[0]-1.0*uM1Self[0]-1.0*uM1Other[0]+4.0*m2r[0]; 
  enRHS[1] = (1.5*alphaE[2]*vtsq_other[3]*m_other)/m_self+(1.5*vtsq_other[2]*alphaE[3]*m_other)/m_self+(1.5*alphaE[0]*vtsq_other[1]*m_other)/m_self+(1.5*vtsq_other[0]*alphaE[1]*m_other)/m_self+(0.25*alphaE[2]*uSumSq[3]*m_diff)/m_self+(0.25*uSumSq[2]*alphaE[3]*m_diff)/m_self+(0.25*alphaE[0]*uSumSq[1]*m_diff)/m_self+(0.25*uSumSq[0]*alphaE[1]*m_diff)/m_self-1.5*alphaE[2]*vtsq_self[3]-1.5*vtsq_self[2]*alphaE[3]-1.5*alphaE[0]*vtsq_self[1]-1.0*uM1Self[1]-1.0*uM1Other[1]+4.0*m2r[1]-1.5*vtsq_self[0]*alphaE[1]; 
  enRHS[2] = (1.5*alphaE[1]*vtsq_other[3]*m_other)/m_self+(1.5*vtsq_other[1]*alphaE[3]*m_other)/m_self+(1.5*alphaE[0]*vtsq_other[2]*m_other)/m_self+(1.5*vtsq_other[0]*alphaE[2]*m_other)/m_self+(0.25*alphaE[1]*uSumSq[3]*m_diff)/m_self+(0.25*uSumSq[1]*alphaE[3]*m_diff)/m_self+(0.25*alphaE[0]*uSumSq[2]*m_diff)/m_self+(0.25*uSumSq[0]*alphaE[2]*m_diff)/m_self-1.5*alphaE[1]*vtsq_self[3]-1.5*vtsq_self[1]*alphaE[3]-1.5*alphaE[0]*vtsq_self[2]-1.0*uM1Self[2]-1.0*uM1Other[2]+4.0*m2r[2]-1.5*vtsq_self[0]*alphaE[2]; 
  enRHS[3] = (1.5*alphaE[0]*vtsq_other[3]*m_other)/m_self+(1.5*vtsq_other[0]*alphaE[3]*m_other)/m_self+(1.5*alphaE[1]*vtsq_other[2]*m_other)/m_self+(1.5*vtsq_other[1]*alphaE[2]*m_other)/m_self+(0.25*alphaE[0]*uSumSq[3]*m_diff)/m_self+(0.25*uSumSq[0]*alphaE[3]*m_diff)/m_self+(0.25*alphaE[1]*uSumSq[2]*m_diff)/m_self+(0.25*uSumSq[1]*alphaE[2]*m_diff)/m_self-1.5*alphaE[0]*vtsq_self[3]-1.0*uM1Self[3]-1.0*uM1Other[3]+4.0*m2r[3]-1.5*vtsq_self[0]*alphaE[3]-1.5*alphaE[1]*vtsq_self[2]-1.5*vtsq_self[1]*alphaE[2]; 
 
  gkyl_mat_set(rhs,0,0,momRHS[0]); 
  gkyl_mat_set(rhs,1,0,momRHS[1]); 
  gkyl_mat_set(rhs,2,0,momRHS[2]); 
  gkyl_mat_set(rhs,3,0,momRHS[3]); 
  gkyl_mat_set(rhs,4,0,momRHS[4]); 
  gkyl_mat_set(rhs,5,0,momRHS[5]); 
  gkyl_mat_set(rhs,6,0,momRHS[6]); 
  gkyl_mat_set(rhs,7,0,momRHS[7]); 
  gkyl_mat_set(rhs,8,0,momRHS[8]); 
  gkyl_mat_set(rhs,9,0,momRHS[9]); 
  gkyl_mat_set(rhs,10,0,momRHS[10]); 
  gkyl_mat_set(rhs,11,0,momRHS[11]); 
  gkyl_mat_set(rhs,12,0,enRHS[0]); 
  gkyl_mat_set(rhs,13,0,enRHS[1]); 
  gkyl_mat_set(rhs,14,0,enRHS[2]); 
  gkyl_mat_set(rhs,15,0,enRHS[3]); 
} 
 
