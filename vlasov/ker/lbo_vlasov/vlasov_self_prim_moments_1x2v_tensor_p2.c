#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_self_prim_moments_1x2v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections, const double *nu) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // moms:                 moments of the distribution function (Zeroth, First, and Second in single array). 
  // boundary_corrections: boundary corrections to u and vtSq. 
  // nu:                   collision frequency. 
 
  double m0r[3] = {0.0}; 
  double m1r[6] = {0.0}; 
  double cMr[6] = {0.0}; 
  double cEr[3] = {0.0}; 
  
  if (nu[0] > 0.0) { 
  
  // If m0 or m2 is below zero at a corner, use cell averages.
  bool notCellAvg = true;
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[2]-1.7320508075688772*moms[1]+moms[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[2]+1.7320508075688772*moms[1]+moms[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[11]-1.7320508075688772*moms[10]+moms[9]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[11]+1.7320508075688772*moms[10]+moms[9]) < 0)) notCellAvg = false; 
 
  if (notCellAvg) { 
    m0r[0] = moms[0]; 
    m0r[1] = moms[1]; 
    m0r[2] = moms[2]; 
    m1r[0] = moms[3]; 
    m1r[1] = moms[4]; 
    m1r[2] = moms[5]; 
    m1r[3] = moms[6]; 
    m1r[4] = moms[7]; 
    m1r[5] = moms[8]; 
    gkyl_mat_set(rhs,0,0,moms[3]); 
    gkyl_mat_set(rhs,1,0,moms[4]); 
    gkyl_mat_set(rhs,2,0,moms[5]); 
    gkyl_mat_set(rhs,3,0,moms[6]); 
    gkyl_mat_set(rhs,4,0,moms[7]); 
    gkyl_mat_set(rhs,5,0,moms[8]); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = boundary_corrections[5]; 
    cEr[0] = boundary_corrections[6]; 
    cEr[1] = boundary_corrections[7]; 
    cEr[2] = boundary_corrections[8]; 
    gkyl_mat_set(rhs,6,0,2.0*moms[9]); 
    gkyl_mat_set(rhs,7,0,2.0*moms[10]); 
    gkyl_mat_set(rhs,8,0,2.0*moms[11]); 
  } else { 
    m0r[0] = moms[0]; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = moms[3]; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    gkyl_mat_set(rhs,0,0,moms[3]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    m1r[3] = moms[6]; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    gkyl_mat_set(rhs,0,0,moms[6]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    cEr[0] = boundary_corrections[6]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    gkyl_mat_set(rhs,6,0,2.0*moms[9]); 
    gkyl_mat_set(rhs,7,0,0.0); 
    gkyl_mat_set(rhs,8,0,0.0); 
  }
 
  } else { 
  
    m0r[0] = 1.0; 
    m0r[1] = 0.0; 
    m0r[2] = 0.0; 
    m1r[0] = 1.0; 
    m1r[1] = 0.0; 
    m1r[2] = 0.0; 
    gkyl_mat_set(rhs,0,0,1.0); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    cMr[0] = 0.0; 
    cMr[1] = 0.0; 
    cMr[2] = 0.0; 
    m1r[3] = 1.0; 
    m1r[4] = 0.0; 
    m1r[5] = 0.0; 
    gkyl_mat_set(rhs,0,0,1.0); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    cMr[3] = 0.0; 
    cMr[4] = 0.0; 
    cMr[5] = 0.0; 
    cEr[0] = 0.0; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    gkyl_mat_set(rhs,6,0,1.0); 
    gkyl_mat_set(rhs,7,0,0.0); 
    gkyl_mat_set(rhs,8,0,0.0); 
  
  }
  
  // ....... Block from weak multiply of ux and m0  .......... // 
  gkyl_mat_set(A,0,0,0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,0,1,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,0,2,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,1,0,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,1,1,0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,1,2,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,2,0,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,2,1,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,2,2,0.45175395145262565*m0r[2]+0.7071067811865475*m0r[0]); 
 
  // ....... Block from correction to ux .......... // 
  gkyl_mat_set(A,0,6,-(0.7071067811865475*cMr[0])); 
  gkyl_mat_set(A,0,7,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,0,8,-(0.7071067811865475*cMr[2])); 
  gkyl_mat_set(A,1,6,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,1,7,-(0.6324555320336759*cMr[2])-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,1,8,-(0.6324555320336759*cMr[1])); 
  gkyl_mat_set(A,2,6,-(0.7071067811865475*cMr[2])); 
  gkyl_mat_set(A,2,7,-(0.6324555320336759*cMr[1])); 
  gkyl_mat_set(A,2,8,-(0.45175395145262565*cMr[2])-0.7071067811865475*cMr[0]); 
 
  // ....... Block from weak multiply of ux and m1x  .......... // 
  gkyl_mat_set(A,6,0,0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,6,1,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,6,2,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,7,0,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,7,1,0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,7,2,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,8,0,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,8,1,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,8,2,0.45175395145262565*m1r[2]+0.7071067811865475*m1r[0]); 
 
  // ....... Block from weak multiply of uy and m0  .......... // 
  gkyl_mat_set(A,3,3,0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,3,4,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,3,5,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,4,3,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,4,4,0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,4,5,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,5,3,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,5,4,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,5,5,0.45175395145262565*m0r[2]+0.7071067811865475*m0r[0]); 
 
  // ....... Block from correction to uy .......... // 
  gkyl_mat_set(A,3,6,-(0.7071067811865475*cMr[3])); 
  gkyl_mat_set(A,3,7,-(0.7071067811865475*cMr[4])); 
  gkyl_mat_set(A,3,8,-(0.7071067811865475*cMr[5])); 
  gkyl_mat_set(A,4,6,-(0.7071067811865475*cMr[4])); 
  gkyl_mat_set(A,4,7,-(0.6324555320336759*cMr[5])-0.7071067811865475*cMr[3]); 
  gkyl_mat_set(A,4,8,-(0.6324555320336759*cMr[4])); 
  gkyl_mat_set(A,5,6,-(0.7071067811865475*cMr[5])); 
  gkyl_mat_set(A,5,7,-(0.6324555320336759*cMr[4])); 
  gkyl_mat_set(A,5,8,-(0.45175395145262565*cMr[5])-0.7071067811865475*cMr[3]); 
 
  // ....... Block from weak multiply of uy and m1y  .......... // 
  gkyl_mat_set(A,6,3,0.7071067811865475*m1r[3]); 
  gkyl_mat_set(A,6,4,0.7071067811865475*m1r[4]); 
  gkyl_mat_set(A,6,5,0.7071067811865475*m1r[5]); 
  gkyl_mat_set(A,7,3,0.7071067811865475*m1r[4]); 
  gkyl_mat_set(A,7,4,0.6324555320336759*m1r[5]+0.7071067811865475*m1r[3]); 
  gkyl_mat_set(A,7,5,0.6324555320336759*m1r[4]); 
  gkyl_mat_set(A,8,3,0.7071067811865475*m1r[5]); 
  gkyl_mat_set(A,8,4,0.6324555320336759*m1r[4]); 
  gkyl_mat_set(A,8,5,0.45175395145262565*m1r[5]+0.7071067811865475*m1r[3]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,6,6,1.4142135623730951*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,6,7,1.4142135623730951*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,6,8,1.4142135623730951*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,7,6,1.4142135623730951*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,7,7,1.264911064067352*m0r[2]-0.6324555320336759*cEr[2]+1.4142135623730951*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,7,8,1.264911064067352*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,8,6,1.4142135623730951*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,8,7,1.264911064067352*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,8,8,0.9035079029052515*m0r[2]-0.45175395145262565*cEr[2]+1.4142135623730951*m0r[0]-0.7071067811865475*cEr[0]); 
 
} 
 
