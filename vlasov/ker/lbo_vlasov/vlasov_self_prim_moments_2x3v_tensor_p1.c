#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_self_prim_moments_2x3v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections, const double *nu) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // moms:                 moments of the distribution function (Zeroth, First, and Second in single array). 
  // boundary_corrections: boundary corrections to u and vtSq. 
  // nu:                   collision frequency. 
 
  double m0r[4] = {0.0}; 
  double m1r[12] = {0.0}; 
  double cMr[12] = {0.0}; 
  double cEr[4] = {0.0}; 
  
  if (nu[0] > 0.0) { 
  
  // If m0 or m2 is below zero at a corner, use cell averages.
  bool notCellAvg = true;
  if (notCellAvg && (0.5*(3.0*moms[3]-1.7320508075688772*(moms[2]+moms[1])+moms[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms[3]+1.7320508075688772*moms[2]-1.7320508075688772*moms[1]-1.0*moms[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms[3]-1.7320508075688772*moms[2]+1.7320508075688772*moms[1]-1.0*moms[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*moms[3]+1.7320508075688772*(moms[2]+moms[1])+moms[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*moms[19]-1.7320508075688772*(moms[18]+moms[17])+moms[16]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms[19]+1.7320508075688772*moms[18]-1.7320508075688772*moms[17]-1.0*moms[16])) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.0*moms[19]-1.7320508075688772*moms[18]+1.7320508075688772*moms[17]-1.0*moms[16])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.0*moms[19]+1.7320508075688772*(moms[18]+moms[17])+moms[16]) < 0)) notCellAvg = false; 
 
  if (notCellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m0r[2] = moms[2]; 
    m0r[3] = moms[3]; 
    m1r[0] = moms[4]; 
    m1r[1] = moms[5]; 
    m1r[2] = moms[6]; 
    m1r[3] = moms[7]; 
    m1r[4] = moms[8]; 
    m1r[5] = moms[9]; 
    m1r[6] = moms[10]; 
    m1r[7] = moms[11]; 
    m1r[8] = moms[12]; 
    m1r[9] = moms[13]; 
    m1r[10] = moms[14]; 
    m1r[11] = moms[15]; 
    gkyl_mat_set(rhs,0,0,moms[4]); 
    gkyl_mat_set(rhs,1,0,moms[5]); 
    gkyl_mat_set(rhs,2,0,moms[6]); 
    gkyl_mat_set(rhs,3,0,moms[7]); 
    gkyl_mat_set(rhs,4,0,moms[8]); 
    gkyl_mat_set(rhs,5,0,moms[9]); 
    gkyl_mat_set(rhs,6,0,moms[10]); 
    gkyl_mat_set(rhs,7,0,moms[11]); 
    gkyl_mat_set(rhs,8,0,moms[12]); 
    gkyl_mat_set(rhs,9,0,moms[13]); 
    gkyl_mat_set(rhs,10,0,moms[14]); 
    gkyl_mat_set(rhs,11,0,moms[15]); 
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
    gkyl_mat_set(rhs,12,0,2.0*moms[16]); 
    gkyl_mat_set(rhs,13,0,2.0*moms[17]); 
    gkyl_mat_set(rhs,14,0,2.0*moms[18]); 
    gkyl_mat_set(rhs,15,0,2.0*moms[19]); 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m0r[3] = 0.0; 
    m1r[0] = moms[4]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    m1r[3] = 0.0; 
    gkyl_mat_set(rhs,0,0,moms[4]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    cMr[3] = 0.0; 
    m1r[4] = moms[8]; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    gkyl_mat_set(rhs,0,0,moms[8]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = 0.0; 
    cMr[6] = 0.0; 
    cMr[7] = 0.0; 
    m1r[8] = moms[12]; 
    m1r[9] = 0.0; 
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    gkyl_mat_set(rhs,0,0,moms[12]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    cMr[8] = boundary_corrections[8]; 
    cMr[9] = 0.0; 
    cMr[10] = 0.0; 
    cMr[11] = 0.0; 
    cEr[0] = boundary_corrections[12]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    gkyl_mat_set(rhs,12,0,2.0*moms[16]); 
    gkyl_mat_set(rhs,13,0,0.0); 
    gkyl_mat_set(rhs,14,0,0.0); 
    gkyl_mat_set(rhs,15,0,0.0); 
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
    gkyl_mat_set(rhs,0,0,1.0); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    cMr[0] = 0.0; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    cMr[3] = 0.0; 
    m1r[4] = 1.0; 
    m1r[5] = 0.0; 
    m1r[6] = 0.0; 
    m1r[7] = 0.0; 
    gkyl_mat_set(rhs,0,0,1.0); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    cMr[6] = 0.0; 
    cMr[7] = 0.0; 
    m1r[8] = 1.0; 
    m1r[9] = 0.0; 
    m1r[10] = 0.0; 
    m1r[11] = 0.0; 
    gkyl_mat_set(rhs,0,0,1.0); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
    cMr[8] = 0.0; 
    cMr[9] = 0.0; 
    cMr[10] = 0.0; 
    cMr[11] = 0.0; 
    cEr[0] = 0.0; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    gkyl_mat_set(rhs,12,0,1.0); 
    gkyl_mat_set(rhs,13,0,0.0); 
    gkyl_mat_set(rhs,14,0,0.0); 
    gkyl_mat_set(rhs,15,0,0.0); 
  
  }
  
  // ....... Block from weak multiply of ux and m0  .......... // 
  gkyl_mat_set(A,0,0,0.5*m0r[0]); 
  gkyl_mat_set(A,0,1,0.5*m0r[1]); 
  gkyl_mat_set(A,0,2,0.5*m0r[2]); 
  gkyl_mat_set(A,0,3,0.5*m0r[3]); 
  gkyl_mat_set(A,1,0,0.5*m0r[1]); 
  gkyl_mat_set(A,1,1,0.5*m0r[0]); 
  gkyl_mat_set(A,1,2,0.5*m0r[3]); 
  gkyl_mat_set(A,1,3,0.5*m0r[2]); 
  gkyl_mat_set(A,2,0,0.5*m0r[2]); 
  gkyl_mat_set(A,2,1,0.5*m0r[3]); 
  gkyl_mat_set(A,2,2,0.5*m0r[0]); 
  gkyl_mat_set(A,2,3,0.5*m0r[1]); 
  gkyl_mat_set(A,3,0,0.5*m0r[3]); 
  gkyl_mat_set(A,3,1,0.5*m0r[2]); 
  gkyl_mat_set(A,3,2,0.5*m0r[1]); 
  gkyl_mat_set(A,3,3,0.5*m0r[0]); 
 
  // ....... Block from correction to ux .......... // 
  gkyl_mat_set(A,0,12,-(0.5*cMr[0])); 
  gkyl_mat_set(A,0,13,-(0.5*cMr[1])); 
  gkyl_mat_set(A,0,14,-(0.5*cMr[2])); 
  gkyl_mat_set(A,0,15,-(0.5*cMr[3])); 
  gkyl_mat_set(A,1,12,-(0.5*cMr[1])); 
  gkyl_mat_set(A,1,13,-(0.5*cMr[0])); 
  gkyl_mat_set(A,1,14,-(0.5*cMr[3])); 
  gkyl_mat_set(A,1,15,-(0.5*cMr[2])); 
  gkyl_mat_set(A,2,12,-(0.5*cMr[2])); 
  gkyl_mat_set(A,2,13,-(0.5*cMr[3])); 
  gkyl_mat_set(A,2,14,-(0.5*cMr[0])); 
  gkyl_mat_set(A,2,15,-(0.5*cMr[1])); 
  gkyl_mat_set(A,3,12,-(0.5*cMr[3])); 
  gkyl_mat_set(A,3,13,-(0.5*cMr[2])); 
  gkyl_mat_set(A,3,14,-(0.5*cMr[1])); 
  gkyl_mat_set(A,3,15,-(0.5*cMr[0])); 
 
  // ....... Block from weak multiply of ux and m1x  .......... // 
  gkyl_mat_set(A,12,0,0.5*m1r[0]); 
  gkyl_mat_set(A,12,1,0.5*m1r[1]); 
  gkyl_mat_set(A,12,2,0.5*m1r[2]); 
  gkyl_mat_set(A,12,3,0.5*m1r[3]); 
  gkyl_mat_set(A,13,0,0.5*m1r[1]); 
  gkyl_mat_set(A,13,1,0.5*m1r[0]); 
  gkyl_mat_set(A,13,2,0.5*m1r[3]); 
  gkyl_mat_set(A,13,3,0.5*m1r[2]); 
  gkyl_mat_set(A,14,0,0.5*m1r[2]); 
  gkyl_mat_set(A,14,1,0.5*m1r[3]); 
  gkyl_mat_set(A,14,2,0.5*m1r[0]); 
  gkyl_mat_set(A,14,3,0.5*m1r[1]); 
  gkyl_mat_set(A,15,0,0.5*m1r[3]); 
  gkyl_mat_set(A,15,1,0.5*m1r[2]); 
  gkyl_mat_set(A,15,2,0.5*m1r[1]); 
  gkyl_mat_set(A,15,3,0.5*m1r[0]); 
 
  // ....... Block from weak multiply of uy and m0  .......... // 
  gkyl_mat_set(A,4,4,0.5*m0r[0]); 
  gkyl_mat_set(A,4,5,0.5*m0r[1]); 
  gkyl_mat_set(A,4,6,0.5*m0r[2]); 
  gkyl_mat_set(A,4,7,0.5*m0r[3]); 
  gkyl_mat_set(A,5,4,0.5*m0r[1]); 
  gkyl_mat_set(A,5,5,0.5*m0r[0]); 
  gkyl_mat_set(A,5,6,0.5*m0r[3]); 
  gkyl_mat_set(A,5,7,0.5*m0r[2]); 
  gkyl_mat_set(A,6,4,0.5*m0r[2]); 
  gkyl_mat_set(A,6,5,0.5*m0r[3]); 
  gkyl_mat_set(A,6,6,0.5*m0r[0]); 
  gkyl_mat_set(A,6,7,0.5*m0r[1]); 
  gkyl_mat_set(A,7,4,0.5*m0r[3]); 
  gkyl_mat_set(A,7,5,0.5*m0r[2]); 
  gkyl_mat_set(A,7,6,0.5*m0r[1]); 
  gkyl_mat_set(A,7,7,0.5*m0r[0]); 
 
  // ....... Block from correction to uy .......... // 
  gkyl_mat_set(A,4,12,-(0.5*cMr[4])); 
  gkyl_mat_set(A,4,13,-(0.5*cMr[5])); 
  gkyl_mat_set(A,4,14,-(0.5*cMr[6])); 
  gkyl_mat_set(A,4,15,-(0.5*cMr[7])); 
  gkyl_mat_set(A,5,12,-(0.5*cMr[5])); 
  gkyl_mat_set(A,5,13,-(0.5*cMr[4])); 
  gkyl_mat_set(A,5,14,-(0.5*cMr[7])); 
  gkyl_mat_set(A,5,15,-(0.5*cMr[6])); 
  gkyl_mat_set(A,6,12,-(0.5*cMr[6])); 
  gkyl_mat_set(A,6,13,-(0.5*cMr[7])); 
  gkyl_mat_set(A,6,14,-(0.5*cMr[4])); 
  gkyl_mat_set(A,6,15,-(0.5*cMr[5])); 
  gkyl_mat_set(A,7,12,-(0.5*cMr[7])); 
  gkyl_mat_set(A,7,13,-(0.5*cMr[6])); 
  gkyl_mat_set(A,7,14,-(0.5*cMr[5])); 
  gkyl_mat_set(A,7,15,-(0.5*cMr[4])); 
 
  // ....... Block from weak multiply of uy and m1y  .......... // 
  gkyl_mat_set(A,12,4,0.5*m1r[4]); 
  gkyl_mat_set(A,12,5,0.5*m1r[5]); 
  gkyl_mat_set(A,12,6,0.5*m1r[6]); 
  gkyl_mat_set(A,12,7,0.5*m1r[7]); 
  gkyl_mat_set(A,13,4,0.5*m1r[5]); 
  gkyl_mat_set(A,13,5,0.5*m1r[4]); 
  gkyl_mat_set(A,13,6,0.5*m1r[7]); 
  gkyl_mat_set(A,13,7,0.5*m1r[6]); 
  gkyl_mat_set(A,14,4,0.5*m1r[6]); 
  gkyl_mat_set(A,14,5,0.5*m1r[7]); 
  gkyl_mat_set(A,14,6,0.5*m1r[4]); 
  gkyl_mat_set(A,14,7,0.5*m1r[5]); 
  gkyl_mat_set(A,15,4,0.5*m1r[7]); 
  gkyl_mat_set(A,15,5,0.5*m1r[6]); 
  gkyl_mat_set(A,15,6,0.5*m1r[5]); 
  gkyl_mat_set(A,15,7,0.5*m1r[4]); 
 
  // ....... Block from weak multiply of uz and m0  .......... // 
  gkyl_mat_set(A,8,8,0.5*m0r[0]); 
  gkyl_mat_set(A,8,9,0.5*m0r[1]); 
  gkyl_mat_set(A,8,10,0.5*m0r[2]); 
  gkyl_mat_set(A,8,11,0.5*m0r[3]); 
  gkyl_mat_set(A,9,8,0.5*m0r[1]); 
  gkyl_mat_set(A,9,9,0.5*m0r[0]); 
  gkyl_mat_set(A,9,10,0.5*m0r[3]); 
  gkyl_mat_set(A,9,11,0.5*m0r[2]); 
  gkyl_mat_set(A,10,8,0.5*m0r[2]); 
  gkyl_mat_set(A,10,9,0.5*m0r[3]); 
  gkyl_mat_set(A,10,10,0.5*m0r[0]); 
  gkyl_mat_set(A,10,11,0.5*m0r[1]); 
  gkyl_mat_set(A,11,8,0.5*m0r[3]); 
  gkyl_mat_set(A,11,9,0.5*m0r[2]); 
  gkyl_mat_set(A,11,10,0.5*m0r[1]); 
  gkyl_mat_set(A,11,11,0.5*m0r[0]); 
 
  // ....... Block from correction to uz .......... // 
  gkyl_mat_set(A,8,12,-(0.5*cMr[8])); 
  gkyl_mat_set(A,8,13,-(0.5*cMr[9])); 
  gkyl_mat_set(A,8,14,-(0.5*cMr[10])); 
  gkyl_mat_set(A,8,15,-(0.5*cMr[11])); 
  gkyl_mat_set(A,9,12,-(0.5*cMr[9])); 
  gkyl_mat_set(A,9,13,-(0.5*cMr[8])); 
  gkyl_mat_set(A,9,14,-(0.5*cMr[11])); 
  gkyl_mat_set(A,9,15,-(0.5*cMr[10])); 
  gkyl_mat_set(A,10,12,-(0.5*cMr[10])); 
  gkyl_mat_set(A,10,13,-(0.5*cMr[11])); 
  gkyl_mat_set(A,10,14,-(0.5*cMr[8])); 
  gkyl_mat_set(A,10,15,-(0.5*cMr[9])); 
  gkyl_mat_set(A,11,12,-(0.5*cMr[11])); 
  gkyl_mat_set(A,11,13,-(0.5*cMr[10])); 
  gkyl_mat_set(A,11,14,-(0.5*cMr[9])); 
  gkyl_mat_set(A,11,15,-(0.5*cMr[8])); 
 
  // ....... Block from weak multiply of uz and m1z  .......... // 
  gkyl_mat_set(A,12,8,0.5*m1r[8]); 
  gkyl_mat_set(A,12,9,0.5*m1r[9]); 
  gkyl_mat_set(A,12,10,0.5*m1r[10]); 
  gkyl_mat_set(A,12,11,0.5*m1r[11]); 
  gkyl_mat_set(A,13,8,0.5*m1r[9]); 
  gkyl_mat_set(A,13,9,0.5*m1r[8]); 
  gkyl_mat_set(A,13,10,0.5*m1r[11]); 
  gkyl_mat_set(A,13,11,0.5*m1r[10]); 
  gkyl_mat_set(A,14,8,0.5*m1r[10]); 
  gkyl_mat_set(A,14,9,0.5*m1r[11]); 
  gkyl_mat_set(A,14,10,0.5*m1r[8]); 
  gkyl_mat_set(A,14,11,0.5*m1r[9]); 
  gkyl_mat_set(A,15,8,0.5*m1r[11]); 
  gkyl_mat_set(A,15,9,0.5*m1r[10]); 
  gkyl_mat_set(A,15,10,0.5*m1r[9]); 
  gkyl_mat_set(A,15,11,0.5*m1r[8]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,12,12,1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,12,13,1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,12,14,1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,12,15,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,13,12,1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,13,13,1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,13,14,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,13,15,1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,14,12,1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,14,13,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,14,14,1.5*m0r[0]-0.5*cEr[0]); 
  gkyl_mat_set(A,14,15,1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,15,12,1.5*m0r[3]-0.5*cEr[3]); 
  gkyl_mat_set(A,15,13,1.5*m0r[2]-0.5*cEr[2]); 
  gkyl_mat_set(A,15,14,1.5*m0r[1]-0.5*cEr[1]); 
  gkyl_mat_set(A,15,15,1.5*m0r[0]-0.5*cEr[0]); 
 
} 
 
