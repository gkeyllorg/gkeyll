#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_self_prim_moments_1x1v_tensor_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections, const double *nu) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // moms:                 moments of the distribution function (Zeroth, First, and Second in single array). 
  // boundary_corrections: boundary corrections to u and vtSq. 
  // nu:                   collision frequency. 
 
  double m0r[2] = {0.0}; 
  double m1r[2] = {0.0}; 
  double cMr[2] = {0.0}; 
  double cEr[2] = {0.0}; 
  
  if (nu[0] > 0.0) { 
  
  // If m0 or m2 is below zero at a corner, use cell averages.
  bool notCellAvg = true;
  if (notCellAvg && (-(0.5*(2.4494897427831783*moms[1]-1.4142135623730951*moms[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*moms[1]+1.4142135623730951*moms[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(2.4494897427831783*moms[5]-1.4142135623730951*moms[4])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(2.4494897427831783*moms[5]+1.4142135623730951*moms[4]) < 0)) notCellAvg = false; 
 
  if (notCellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m1r[0] = moms[2]; 
    m1r[1] = moms[3]; 
    gkyl_mat_set(rhs,0,0,moms[2]); 
    gkyl_mat_set(rhs,1,0,moms[3]); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cEr[0] = boundary_corrections[2]; 
    cEr[1] = boundary_corrections[3]; 
    gkyl_mat_set(rhs,2,0,2.0*moms[4]); 
    gkyl_mat_set(rhs,3,0,2.0*moms[5]); 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = 0.0; 
    m1r[0] = moms[2]; 
    m1r[1] = 0.0; 
    gkyl_mat_set(rhs,0,0,moms[2]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    cEr[0] = boundary_corrections[2]; 
    cEr[1] = 0.0; 
    gkyl_mat_set(rhs,2,0,2.0*moms[4]); 
    gkyl_mat_set(rhs,3,0,0.0); 
  }
 
  } else { 
  
    m0r[0] = 1.0; 
    m0r[1] = 0.0; 
    m1r[0] = 1.0; 
    m1r[1] = 0.0; 
    gkyl_mat_set(rhs,0,0,1.0); 
    gkyl_mat_set(rhs,1,0,0.0); 
    cMr[0] = 0.0; 
    cMr[1] = 0.0; 
    cEr[0] = 0.0; 
    cEr[1] = 0.0; 
    gkyl_mat_set(rhs,2,0,1.0); 
    gkyl_mat_set(rhs,3,0,0.0); 
  
  }
  
  // ....... Block from weak multiply of ux and m0  .......... // 
  gkyl_mat_set(A,0,0,0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,0,1,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,1,0,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,1,1,0.7071067811865475*m0r[0]); 
 
  // ....... Block from correction to ux .......... // 
  gkyl_mat_set(A,0,2,-(0.7071067811865475*cMr[0])); 
  gkyl_mat_set(A,0,3,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,1,2,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,1,3,-(0.7071067811865475*cMr[0])); 
 
  // ....... Block from weak multiply of ux and m1x  .......... // 
  gkyl_mat_set(A,2,0,0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,2,1,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,3,0,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,3,1,0.7071067811865475*m1r[0]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,2,2,0.7071067811865475*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,2,3,0.7071067811865475*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,3,2,0.7071067811865475*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,3,3,0.7071067811865475*m0r[0]-0.7071067811865475*cEr[0]); 
 
} 
 
