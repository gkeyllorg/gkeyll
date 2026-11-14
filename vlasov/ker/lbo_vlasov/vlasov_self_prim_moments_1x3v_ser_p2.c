#include <gkyl_prim_lbo_vlasov_kernels.h> 
 
GKYL_CU_DH void vlasov_self_prim_moments_1x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections, const double *nu) 
{ 
  // A:                    Matrix to be inverted to solve Ax = rhs (set by this function). 
  // rhs:                  right-hand side of Ax = rhs (set by this function). 
  // moms:                 moments of the distribution function (Zeroth, First, and Second in single array). 
  // boundary_corrections: boundary corrections to u and vtSq. 
  // nu:                   collision frequency. 
 
  double m0r[3] = {0.0}; 
  double m1r[9] = {0.0}; 
  double cMr[9] = {0.0}; 
  double cEr[3] = {0.0}; 
  
  if (nu[0] > 0.0) { 
  
  // If m0 or m2 is below zero at a corner, use cell averages.
  bool notCellAvg = true;
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[2]-1.7320508075688772*moms[1]+moms[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[2]+1.7320508075688772*moms[1]+moms[0]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[14]-1.7320508075688772*moms[13]+moms[12]) < 0)) notCellAvg = false; 
  if (notCellAvg && (0.7071067811865475*(2.23606797749979*moms[14]+1.7320508075688772*moms[13]+moms[12]) < 0)) notCellAvg = false; 
 
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
    m1r[6] = moms[9]; 
    m1r[7] = moms[10]; 
    m1r[8] = moms[11]; 
    gkyl_mat_set(rhs,0,0,moms[3]); 
    gkyl_mat_set(rhs,1,0,moms[4]); 
    gkyl_mat_set(rhs,2,0,moms[5]); 
    gkyl_mat_set(rhs,3,0,moms[6]); 
    gkyl_mat_set(rhs,4,0,moms[7]); 
    gkyl_mat_set(rhs,5,0,moms[8]); 
    gkyl_mat_set(rhs,6,0,moms[9]); 
    gkyl_mat_set(rhs,7,0,moms[10]); 
    gkyl_mat_set(rhs,8,0,moms[11]); 
    cMr[0] = boundary_corrections[0]; 
    cMr[1] = boundary_corrections[1]; 
    cMr[2] = boundary_corrections[2]; 
    cMr[3] = boundary_corrections[3]; 
    cMr[4] = boundary_corrections[4]; 
    cMr[5] = boundary_corrections[5]; 
    cMr[6] = boundary_corrections[6]; 
    cMr[7] = boundary_corrections[7]; 
    cMr[8] = boundary_corrections[8]; 
    cEr[0] = boundary_corrections[9]; 
    cEr[1] = boundary_corrections[10]; 
    cEr[2] = boundary_corrections[11]; 
    gkyl_mat_set(rhs,9,0,2.0*moms[12]); 
    gkyl_mat_set(rhs,10,0,2.0*moms[13]); 
    gkyl_mat_set(rhs,11,0,2.0*moms[14]); 
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
    m1r[6] = moms[9]; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    gkyl_mat_set(rhs,0,0,moms[9]); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    cMr[6] = boundary_corrections[6]; 
    cMr[7] = 0.0; 
    cMr[8] = 0.0; 
    cEr[0] = boundary_corrections[9]; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    gkyl_mat_set(rhs,9,0,2.0*moms[12]); 
    gkyl_mat_set(rhs,10,0,0.0); 
    gkyl_mat_set(rhs,11,0,0.0); 
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
    m1r[6] = 1.0; 
    m1r[7] = 0.0; 
    m1r[8] = 0.0; 
    gkyl_mat_set(rhs,0,0,1.0); 
    gkyl_mat_set(rhs,1,0,0.0); 
    gkyl_mat_set(rhs,2,0,0.0); 
    cMr[6] = 0.0; 
    cMr[7] = 0.0; 
    cMr[8] = 0.0; 
    cEr[0] = 0.0; 
    cEr[1] = 0.0; 
    cEr[2] = 0.0; 
    gkyl_mat_set(rhs,9,0,1.0); 
    gkyl_mat_set(rhs,10,0,0.0); 
    gkyl_mat_set(rhs,11,0,0.0); 
  
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
  gkyl_mat_set(A,0,9,-(0.7071067811865475*cMr[0])); 
  gkyl_mat_set(A,0,10,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,0,11,-(0.7071067811865475*cMr[2])); 
  gkyl_mat_set(A,1,9,-(0.7071067811865475*cMr[1])); 
  gkyl_mat_set(A,1,10,-(0.6324555320336759*cMr[2])-0.7071067811865475*cMr[0]); 
  gkyl_mat_set(A,1,11,-(0.6324555320336759*cMr[1])); 
  gkyl_mat_set(A,2,9,-(0.7071067811865475*cMr[2])); 
  gkyl_mat_set(A,2,10,-(0.6324555320336759*cMr[1])); 
  gkyl_mat_set(A,2,11,-(0.45175395145262565*cMr[2])-0.7071067811865475*cMr[0]); 
 
  // ....... Block from weak multiply of ux and m1x  .......... // 
  gkyl_mat_set(A,9,0,0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,9,1,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,9,2,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,10,0,0.7071067811865475*m1r[1]); 
  gkyl_mat_set(A,10,1,0.6324555320336759*m1r[2]+0.7071067811865475*m1r[0]); 
  gkyl_mat_set(A,10,2,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,11,0,0.7071067811865475*m1r[2]); 
  gkyl_mat_set(A,11,1,0.6324555320336759*m1r[1]); 
  gkyl_mat_set(A,11,2,0.45175395145262565*m1r[2]+0.7071067811865475*m1r[0]); 
 
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
  gkyl_mat_set(A,3,9,-(0.7071067811865475*cMr[3])); 
  gkyl_mat_set(A,3,10,-(0.7071067811865475*cMr[4])); 
  gkyl_mat_set(A,3,11,-(0.7071067811865475*cMr[5])); 
  gkyl_mat_set(A,4,9,-(0.7071067811865475*cMr[4])); 
  gkyl_mat_set(A,4,10,-(0.6324555320336759*cMr[5])-0.7071067811865475*cMr[3]); 
  gkyl_mat_set(A,4,11,-(0.6324555320336759*cMr[4])); 
  gkyl_mat_set(A,5,9,-(0.7071067811865475*cMr[5])); 
  gkyl_mat_set(A,5,10,-(0.6324555320336759*cMr[4])); 
  gkyl_mat_set(A,5,11,-(0.45175395145262565*cMr[5])-0.7071067811865475*cMr[3]); 
 
  // ....... Block from weak multiply of uy and m1y  .......... // 
  gkyl_mat_set(A,9,3,0.7071067811865475*m1r[3]); 
  gkyl_mat_set(A,9,4,0.7071067811865475*m1r[4]); 
  gkyl_mat_set(A,9,5,0.7071067811865475*m1r[5]); 
  gkyl_mat_set(A,10,3,0.7071067811865475*m1r[4]); 
  gkyl_mat_set(A,10,4,0.6324555320336759*m1r[5]+0.7071067811865475*m1r[3]); 
  gkyl_mat_set(A,10,5,0.6324555320336759*m1r[4]); 
  gkyl_mat_set(A,11,3,0.7071067811865475*m1r[5]); 
  gkyl_mat_set(A,11,4,0.6324555320336759*m1r[4]); 
  gkyl_mat_set(A,11,5,0.45175395145262565*m1r[5]+0.7071067811865475*m1r[3]); 
 
  // ....... Block from weak multiply of uz and m0  .......... // 
  gkyl_mat_set(A,6,6,0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,6,7,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,6,8,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,7,6,0.7071067811865475*m0r[1]); 
  gkyl_mat_set(A,7,7,0.6324555320336759*m0r[2]+0.7071067811865475*m0r[0]); 
  gkyl_mat_set(A,7,8,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,8,6,0.7071067811865475*m0r[2]); 
  gkyl_mat_set(A,8,7,0.6324555320336759*m0r[1]); 
  gkyl_mat_set(A,8,8,0.45175395145262565*m0r[2]+0.7071067811865475*m0r[0]); 
 
  // ....... Block from correction to uz .......... // 
  gkyl_mat_set(A,6,9,-(0.7071067811865475*cMr[6])); 
  gkyl_mat_set(A,6,10,-(0.7071067811865475*cMr[7])); 
  gkyl_mat_set(A,6,11,-(0.7071067811865475*cMr[8])); 
  gkyl_mat_set(A,7,9,-(0.7071067811865475*cMr[7])); 
  gkyl_mat_set(A,7,10,-(0.6324555320336759*cMr[8])-0.7071067811865475*cMr[6]); 
  gkyl_mat_set(A,7,11,-(0.6324555320336759*cMr[7])); 
  gkyl_mat_set(A,8,9,-(0.7071067811865475*cMr[8])); 
  gkyl_mat_set(A,8,10,-(0.6324555320336759*cMr[7])); 
  gkyl_mat_set(A,8,11,-(0.45175395145262565*cMr[8])-0.7071067811865475*cMr[6]); 
 
  // ....... Block from weak multiply of uz and m1z  .......... // 
  gkyl_mat_set(A,9,6,0.7071067811865475*m1r[6]); 
  gkyl_mat_set(A,9,7,0.7071067811865475*m1r[7]); 
  gkyl_mat_set(A,9,8,0.7071067811865475*m1r[8]); 
  gkyl_mat_set(A,10,6,0.7071067811865475*m1r[7]); 
  gkyl_mat_set(A,10,7,0.6324555320336759*m1r[8]+0.7071067811865475*m1r[6]); 
  gkyl_mat_set(A,10,8,0.6324555320336759*m1r[7]); 
  gkyl_mat_set(A,11,6,0.7071067811865475*m1r[8]); 
  gkyl_mat_set(A,11,7,0.6324555320336759*m1r[7]); 
  gkyl_mat_set(A,11,8,0.45175395145262565*m1r[8]+0.7071067811865475*m1r[6]); 
 
  // ....... Block from correction to vtSq .......... // 
  gkyl_mat_set(A,9,9,2.1213203435596424*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,9,10,2.1213203435596424*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,9,11,2.1213203435596424*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,10,9,2.1213203435596424*m0r[1]-0.7071067811865475*cEr[1]); 
  gkyl_mat_set(A,10,10,1.8973665961010278*m0r[2]-0.6324555320336759*cEr[2]+2.1213203435596424*m0r[0]-0.7071067811865475*cEr[0]); 
  gkyl_mat_set(A,10,11,1.8973665961010278*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,11,9,2.1213203435596424*m0r[2]-0.7071067811865475*cEr[2]); 
  gkyl_mat_set(A,11,10,1.8973665961010278*m0r[1]-0.6324555320336759*cEr[1]); 
  gkyl_mat_set(A,11,11,1.355261854357877*m0r[2]-0.45175395145262565*cEr[2]+2.1213203435596424*m0r[0]-0.7071067811865475*cEr[0]); 
 
} 
 
