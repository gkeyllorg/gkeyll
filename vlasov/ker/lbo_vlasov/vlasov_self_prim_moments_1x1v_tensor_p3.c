#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_self_prim_moments_1x1v_tensor_p3(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections, const double *nu) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // moms:                 moments of the distribution function (Zeroth, First, and Second in single array). 
  // boundary_corrections: boundary corrections to u and vtSq. 
  // nu:                   collision frequency. 
 
  double m0r[4] = {0.0}; 
  double m1r[4] = {0.0}; 
  double cMr[4] = {0.0}; 
  double cEr[4] = {0.0}; 
  
  if (nu[0] > 0.0) { 
  
  // If m0 or m2 is below zero at a corner, use cell averages.
  bool notCellAvg = true;
  if (notCellAvg && (-(0.5*(3.7416573867739418*moms[3]-3.1622776601683795*moms[2]+2.4494897427831783*moms[1]-1.4142135623730951*moms[0])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.7416573867739418*moms[3]+3.1622776601683795*moms[2]+2.4494897427831783*moms[1]+1.4142135623730951*moms[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (-(0.5*(3.7416573867739418*moms[11]-3.1622776601683795*moms[10]+2.4494897427831783*moms[9]-1.4142135623730951*moms[8])) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.5*(3.7416573867739418*moms[11]+3.1622776601683795*moms[10]+2.4494897427831783*moms[9]+1.4142135623730951*moms[8]) < 0)) notCellAvg = false; 
 
  if (notCellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m0r[2] = moms[2]; 
    m0r[3] = moms[3]; 
    m1r[0] = moms[4]; 
    m1r[1] = moms[5]; 
    m1r[2] = moms[6]; 
    m1r[3] = moms[7]; 
    gkyl_mat_set(rhs,0,0,moms[4]); 
    gkyl_mat_set(rhs,1,0,moms[5]); 
    gkyl_mat_set(rhs,2,0,moms[6]); 
    gkyl_mat_set(rhs,3,0,moms[7]); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = boundary_corrections[3]; 
    cEr[0] = boundary_corrections[4]; 
    cEr[1] = boundary_corrections[5]; 
    cEr[2] = boundary_corrections[6]; 
    cEr[3] = boundary_corrections[7]; 
    gkyl_mat_set(rhs,4,0,2.0*moms[8]); 
    gkyl_mat_set(rhs,5,0,2.0*moms[9]); 
    gkyl_mat_set(rhs,6,0,2.0*moms[10]); 
    gkyl_mat_set(rhs,7,0,2.0*moms[11]); 
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
    cEr[0] = boundary_corrections[4]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    gkyl_mat_set(rhs,4,0,2.0*moms[8]); 
    gkyl_mat_set(rhs,5,0,0.0); 
    gkyl_mat_set(rhs,6,0,0.0); 
    gkyl_mat_set(rhs,7,0,0.0); 
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
    cEr[0] = 0.0; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    cEr[3] = 0.0; 
    gkyl_mat_set(rhs,4,0,1.0); 
    gkyl_mat_set(rhs,5,0,0.0); 
    gkyl_mat_set(rhs,6,0,0.0); 
    gkyl_mat_set(rhs,7,0,0.0); 
  
  }
  
  // ....... Block from weak multiply of ux and m0  .......... // 
  gkyl_mat_set(A,0,0,0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,0,1,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,0,2,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,0,3,0.7071067811865475*m0r[3]); 
  gkyl_mat_set(A,1,0,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,1,1,0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,1,2,0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,1,3,0.6210590034081186*m0r[2]); 
  gkyl_mat_set(A,2,0,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,2,1,0.6210590034081186*m0r[3]+0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,2,2,0.45175395145262565*m0r[2]+0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,2,3,0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]); 
  gkyl_mat_set(A,3,0,0.7071067811865475*m0r[3]); 
  gkyl_mat_set(A,3,1,0.6210590034081186*m0r[2]); 
  gkyl_mat_set(A,3,2,0.421637021355784*m0r[3]+0.6210590034081186*m0r[1]); 
  gkyl_mat_set(A,3,3,0.421637021355784*m0r[2]+0.7071067811865475*m0r[0]); 
 
  // ....... Block from correction to ux .......... // 
  gkyl_mat_set(A,0,4,-(0.7071067811865475*cMr[0])); 
  gkyl_mat_set(A,0,5,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,0,6,-(0.7071067811865475*cMr[2])); 
  gkyl_mat_set(A,0,7,-(0.7071067811865475*cMr[3])); 
  gkyl_mat_set(A,1,4,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,1,5,-(0.6324555320336759*cMr[2])-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,1,6,-(0.6210590034081186*cMr[3])-0.6324555320336759*cMr[1]); 
  gkyl_mat_set(A,1,7,-(0.6210590034081186*cMr[2])); 
  gkyl_mat_set(A,2,4,-(0.7071067811865475*cMr[2])); 
  gkyl_mat_set(A,2,5,-(0.6210590034081186*cMr[3])-0.6324555320336759*cMr[1]); 
  gkyl_mat_set(A,2,6,-(0.45175395145262565*cMr[2])-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,2,7,-(0.421637021355784*cMr[3])-0.6210590034081186*cMr[1]); 
  gkyl_mat_set(A,3,4,-(0.7071067811865475*cMr[3])); 
  gkyl_mat_set(A,3,5,-(0.6210590034081186*cMr[2])); 
  gkyl_mat_set(A,3,6,-(0.421637021355784*cMr[3])-0.6210590034081186*cMr[1]); 
  gkyl_mat_set(A,3,7,-(0.421637021355784*cMr[2])-0.7071067811865475*cMr[0]); 
 
  // ....... Block from weak multiply of ux and m1x  .......... // 
  gkyl_mat_set(A,4,0,0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,4,1,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,4,2,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,4,3,0.7071067811865475*m1r[3]); 
  gkyl_mat_set(A,5,0,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,5,1,0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,5,2,0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,5,3,0.6210590034081186*m1r[2]); 
  gkyl_mat_set(A,6,0,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,6,1,0.6210590034081186*m1r[3]+0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,6,2,0.45175395145262565*m1r[2]+0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,6,3,0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]); 
  gkyl_mat_set(A,7,0,0.7071067811865475*m1r[3]); 
  gkyl_mat_set(A,7,1,0.6210590034081186*m1r[2]); 
  gkyl_mat_set(A,7,2,0.421637021355784*m1r[3]+0.6210590034081186*m1r[1]); 
  gkyl_mat_set(A,7,3,0.421637021355784*m1r[2]+0.7071067811865475*m1r[0]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,4,4,0.7071067811865475*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,4,5,0.7071067811865475*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,4,6,0.7071067811865475*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,4,7,0.7071067811865475*m0r[3]-0.7071067811865475*cEr[3]); 
  gkyl_mat_set(A,5,4,0.7071067811865475*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,5,5,0.6324555320336759*m0r[2]-0.6324555320336759*cEr[2]+0.7071067811865475*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,5,6,0.6210590034081186*m0r[3]-0.6210590034081186*cEr[3]+0.6324555320336759*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,5,7,0.6210590034081186*m0r[2]-0.6210590034081186*cEr[2]); 
  gkyl_mat_set(A,6,4,0.7071067811865475*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,6,5,0.6210590034081186*m0r[3]-0.6210590034081186*cEr[3]+0.6324555320336759*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,6,6,0.45175395145262565*m0r[2]-0.45175395145262565*cEr[2]+0.7071067811865475*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,6,7,0.421637021355784*m0r[3]-0.421637021355784*cEr[3]+0.6210590034081186*m0r[1]-0.6210590034081186*cEr[1]); 
  gkyl_mat_set(A,7,4,0.7071067811865475*m0r[3]-0.7071067811865475*cEr[3]); 
  gkyl_mat_set(A,7,5,0.6210590034081186*m0r[2]-0.6210590034081186*cEr[2]); 
  gkyl_mat_set(A,7,6,0.421637021355784*m0r[3]-0.421637021355784*cEr[3]+0.6210590034081186*m0r[1]-0.6210590034081186*cEr[1]); 
  gkyl_mat_set(A,7,7,0.421637021355784*m0r[2]-0.421637021355784*cEr[2]+0.7071067811865475*m0r[0]-0.7071067811865475*cEr[0]); 
 
} 
 
