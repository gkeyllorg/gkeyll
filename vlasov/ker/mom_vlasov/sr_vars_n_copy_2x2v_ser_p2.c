#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_basis_ser_2x_p2_exp_sq.h> 
#include <gkyl_basis_ser_2x_p2_sqrt.h> 
#include <gkyl_basis_ser_2x_p1_sqrt.h> 
GKYL_CU_DH void sr_vars_n_copy_2x2v_ser_p2(int count, struct gkyl_nmat *x, const double *M0, double* GKYL_RESTRICT n) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // x:     Input solution vector. 
  // M0:    Lab frame density = Gamma*n.
  // n:     Rest-frame density computed as Gamma_inv*M0 where Gamma_inv = sqrt(1 - |V_drift|^2). 
 
  struct gkyl_mat x0 = gkyl_nmat_get(x, count+0); 
  double V_0[8] = {0.0}; 
  V_0[0] = gkyl_mat_get(&x0,0,0); 
  V_0[1] = gkyl_mat_get(&x0,1,0); 
  V_0[2] = gkyl_mat_get(&x0,2,0); 
  V_0[3] = gkyl_mat_get(&x0,3,0); 
  V_0[4] = gkyl_mat_get(&x0,4,0); 
  V_0[5] = gkyl_mat_get(&x0,5,0); 
  V_0[6] = gkyl_mat_get(&x0,6,0); 
  V_0[7] = gkyl_mat_get(&x0,7,0); 
  double V_0_sq[8] = {0.0}; 
  ser_2x_p2_exp_sq(V_0, V_0_sq); 
 
  struct gkyl_mat x1 = gkyl_nmat_get(x, count+1); 
  double V_1[8] = {0.0}; 
  V_1[0] = gkyl_mat_get(&x1,0,0); 
  V_1[1] = gkyl_mat_get(&x1,1,0); 
  V_1[2] = gkyl_mat_get(&x1,2,0); 
  V_1[3] = gkyl_mat_get(&x1,3,0); 
  V_1[4] = gkyl_mat_get(&x1,4,0); 
  V_1[5] = gkyl_mat_get(&x1,5,0); 
  V_1[6] = gkyl_mat_get(&x1,6,0); 
  V_1[7] = gkyl_mat_get(&x1,7,0); 
  double V_1_sq[8] = {0.0}; 
  ser_2x_p2_exp_sq(V_1, V_1_sq); 
 
  double Gamma2_inv[8] = {0.0}; 
  Gamma2_inv[0] = -(1.0*V_1_sq[0])-1.0*V_0_sq[0]+2.0; 
  Gamma2_inv[1] = -(1.0*V_1_sq[1])-1.0*V_0_sq[1]; 
  Gamma2_inv[2] = -(1.0*V_1_sq[2])-1.0*V_0_sq[2]; 
  Gamma2_inv[3] = -(1.0*V_1_sq[3])-1.0*V_0_sq[3]; 
  Gamma2_inv[4] = -(1.0*V_1_sq[4])-1.0*V_0_sq[4]; 
  Gamma2_inv[5] = -(1.0*V_1_sq[5])-1.0*V_0_sq[5]; 
  Gamma2_inv[6] = -(1.0*V_1_sq[6])-1.0*V_0_sq[6]; 
  Gamma2_inv[7] = -(1.0*V_1_sq[7])-1.0*V_0_sq[7]; 

  int cell_avg = 0;
  if (-(0.5999999999999995*Gamma2_inv[7])-0.5999999999999999*Gamma2_inv[6]+0.4472135954999579*Gamma2_inv[5]+0.4472135954999579*Gamma2_inv[4]+0.9*Gamma2_inv[3]-0.6708203932499369*Gamma2_inv[2]-0.6708203932499369*Gamma2_inv[1]+0.5*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (0.75*Gamma2_inv[7]-0.5590169943749475*Gamma2_inv[5]+0.4472135954999579*Gamma2_inv[4]-0.6708203932499369*Gamma2_inv[1]+0.5*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (-(0.5999999999999995*Gamma2_inv[7])+0.5999999999999999*Gamma2_inv[6]+0.4472135954999579*Gamma2_inv[5]+0.4472135954999579*Gamma2_inv[4]-0.9*Gamma2_inv[3]+0.6708203932499369*Gamma2_inv[2]-0.6708203932499369*Gamma2_inv[1]+0.5*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (0.75*Gamma2_inv[6]+0.4472135954999579*Gamma2_inv[5]-0.5590169943749475*Gamma2_inv[4]-0.6708203932499369*Gamma2_inv[2]+0.5*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (-(0.5590169943749475*Gamma2_inv[5])-0.5590169943749475*Gamma2_inv[4]+0.5*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (-(0.75*Gamma2_inv[6])+0.4472135954999579*Gamma2_inv[5]-0.5590169943749475*Gamma2_inv[4]+0.6708203932499369*Gamma2_inv[2]+0.5*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (0.5999999999999995*Gamma2_inv[7]-0.5999999999999999*Gamma2_inv[6]+0.4472135954999579*Gamma2_inv[5]+0.4472135954999579*Gamma2_inv[4]-0.9*Gamma2_inv[3]-0.6708203932499369*Gamma2_inv[2]+0.6708203932499369*Gamma2_inv[1]+0.5*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (-(0.75*Gamma2_inv[7])-0.5590169943749475*Gamma2_inv[5]+0.4472135954999579*Gamma2_inv[4]+0.6708203932499369*Gamma2_inv[1]+0.5*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (0.5999999999999995*Gamma2_inv[7]+0.5999999999999999*Gamma2_inv[6]+0.4472135954999579*Gamma2_inv[5]+0.4472135954999579*Gamma2_inv[4]+0.9*Gamma2_inv[3]+0.6708203932499369*Gamma2_inv[2]+0.6708203932499369*Gamma2_inv[1]+0.5*Gamma2_inv[0] < 0.0) cell_avg = 1; 
  if (cell_avg) { 
    double Gamma2_inv_lobatto[4] = {0.0}; 
    double Gamma2_inv_p1[4] = {0.0}; 
    double Gamma_inv_p1[4] = {0.0}; 
    double V0_quad = 0.0; 
    double V1_quad = 0.0; 

    V0_quad = -(1.9364916731037085*V_0[7])-1.9364916731037085*V_0[6]+1.118033988749895*V_0[5]+1.118033988749895*V_0[4]+1.5*V_0[3]-0.8660254037844386*V_0[2]-0.8660254037844386*V_0[1]+0.5*V_0[0]; 
    V1_quad = -(1.9364916731037085*V_1[7])-1.9364916731037085*V_1[6]+1.118033988749895*V_1[5]+1.118033988749895*V_1[4]+1.5*V_1[3]-0.8660254037844386*V_1[2]-0.8660254037844386*V_1[1]+0.5*V_1[0]; 
    Gamma2_inv_lobatto[0] = 1.0 - V0_quad*V0_quad - V1_quad*V1_quad; 
    if (Gamma2_inv_lobatto[0] < 0.0) Gamma2_inv_lobatto[0] = 1.0e-16; 

    V0_quad = 1.9364916731037085*V_0[7]-1.9364916731037085*V_0[6]+1.118033988749895*V_0[5]+1.118033988749895*V_0[4]-1.5*V_0[3]-0.8660254037844386*V_0[2]+0.8660254037844386*V_0[1]+0.5*V_0[0]; 
    V1_quad = 1.9364916731037085*V_1[7]-1.9364916731037085*V_1[6]+1.118033988749895*V_1[5]+1.118033988749895*V_1[4]-1.5*V_1[3]-0.8660254037844386*V_1[2]+0.8660254037844386*V_1[1]+0.5*V_1[0]; 
    Gamma2_inv_lobatto[1] = 1.0 - V0_quad*V0_quad - V1_quad*V1_quad; 
    if (Gamma2_inv_lobatto[1] < 0.0) Gamma2_inv_lobatto[1] = 1.0e-16; 

    V0_quad = -(1.9364916731037085*V_0[7])+1.9364916731037085*V_0[6]+1.118033988749895*V_0[5]+1.118033988749895*V_0[4]-1.5*V_0[3]+0.8660254037844386*V_0[2]-0.8660254037844386*V_0[1]+0.5*V_0[0]; 
    V1_quad = -(1.9364916731037085*V_1[7])+1.9364916731037085*V_1[6]+1.118033988749895*V_1[5]+1.118033988749895*V_1[4]-1.5*V_1[3]+0.8660254037844386*V_1[2]-0.8660254037844386*V_1[1]+0.5*V_1[0]; 
    Gamma2_inv_lobatto[2] = 1.0 - V0_quad*V0_quad - V1_quad*V1_quad; 
    if (Gamma2_inv_lobatto[2] < 0.0) Gamma2_inv_lobatto[2] = 1.0e-16; 

    V0_quad = 1.9364916731037085*V_0[7]+1.9364916731037085*V_0[6]+1.118033988749895*V_0[5]+1.118033988749895*V_0[4]+1.5*V_0[3]+0.8660254037844386*V_0[2]+0.8660254037844386*V_0[1]+0.5*V_0[0]; 
    V1_quad = 1.9364916731037085*V_1[7]+1.9364916731037085*V_1[6]+1.118033988749895*V_1[5]+1.118033988749895*V_1[4]+1.5*V_1[3]+0.8660254037844386*V_1[2]+0.8660254037844386*V_1[1]+0.5*V_1[0]; 
    Gamma2_inv_lobatto[3] = 1.0 - V0_quad*V0_quad - V1_quad*V1_quad; 
    if (Gamma2_inv_lobatto[3] < 0.0) Gamma2_inv_lobatto[3] = 1.0e-16; 

    Gamma2_inv_p1[0] = 0.5*Gamma2_inv_lobatto[3]+0.5*Gamma2_inv_lobatto[2]+0.5*Gamma2_inv_lobatto[1]+0.5*Gamma2_inv_lobatto[0]; 
    Gamma2_inv_p1[1] = 0.2886751345948129*Gamma2_inv_lobatto[3]-0.2886751345948129*Gamma2_inv_lobatto[2]+0.2886751345948129*Gamma2_inv_lobatto[1]-0.2886751345948129*Gamma2_inv_lobatto[0]; 
    Gamma2_inv_p1[2] = 0.2886751345948129*Gamma2_inv_lobatto[3]+0.2886751345948129*Gamma2_inv_lobatto[2]-0.2886751345948129*Gamma2_inv_lobatto[1]-0.2886751345948129*Gamma2_inv_lobatto[0]; 
    Gamma2_inv_p1[3] = 0.16666666666666666*Gamma2_inv_lobatto[3]-0.16666666666666666*Gamma2_inv_lobatto[2]-0.16666666666666666*Gamma2_inv_lobatto[1]+0.16666666666666666*Gamma2_inv_lobatto[0]; 
    ser_2x_p1_sqrt(Gamma2_inv_p1, Gamma_inv_p1); 
    n[0] = 0.5*Gamma_inv_p1[3]*M0[3]+0.5*Gamma_inv_p1[2]*M0[2]+0.5*Gamma_inv_p1[1]*M0[1]+0.5*Gamma_inv_p1[0]*M0[0]; 
    n[1] = 0.44721359549995804*Gamma_inv_p1[3]*M0[6]+0.4472135954999579*Gamma_inv_p1[1]*M0[4]+0.5*Gamma_inv_p1[2]*M0[3]+0.5*M0[2]*Gamma_inv_p1[3]+0.5*Gamma_inv_p1[0]*M0[1]+0.5*M0[0]*Gamma_inv_p1[1]; 
    n[2] = 0.44721359549995804*Gamma_inv_p1[3]*M0[7]+0.4472135954999579*Gamma_inv_p1[2]*M0[5]+0.5*Gamma_inv_p1[1]*M0[3]+0.5*M0[1]*Gamma_inv_p1[3]+0.5*Gamma_inv_p1[0]*M0[2]+0.5*M0[0]*Gamma_inv_p1[2]; 
    n[3] = 0.44721359549995804*Gamma_inv_p1[2]*M0[7]+0.44721359549995804*Gamma_inv_p1[1]*M0[6]+0.4472135954999579*Gamma_inv_p1[3]*M0[5]+0.4472135954999579*Gamma_inv_p1[3]*M0[4]+0.5*Gamma_inv_p1[0]*M0[3]+0.5*M0[0]*Gamma_inv_p1[3]+0.5*Gamma_inv_p1[1]*M0[2]+0.5*M0[1]*Gamma_inv_p1[2]; 
    n[4] = 0.5000000000000001*Gamma_inv_p1[2]*M0[6]+0.5*Gamma_inv_p1[0]*M0[4]+0.4472135954999579*Gamma_inv_p1[3]*M0[3]+0.4472135954999579*Gamma_inv_p1[1]*M0[1]; 
    n[5] = 0.5000000000000001*Gamma_inv_p1[1]*M0[7]+0.5*Gamma_inv_p1[0]*M0[5]+0.4472135954999579*Gamma_inv_p1[3]*M0[3]+0.4472135954999579*Gamma_inv_p1[2]*M0[2]; 
    n[6] = 0.4*Gamma_inv_p1[3]*M0[7]+0.5*Gamma_inv_p1[0]*M0[6]+0.5000000000000001*Gamma_inv_p1[2]*M0[4]+0.44721359549995804*Gamma_inv_p1[1]*M0[3]+0.44721359549995804*M0[1]*Gamma_inv_p1[3]; 
    n[7] = 0.5*Gamma_inv_p1[0]*M0[7]+0.4*Gamma_inv_p1[3]*M0[6]+0.5000000000000001*Gamma_inv_p1[1]*M0[5]+0.44721359549995804*Gamma_inv_p1[2]*M0[3]+0.44721359549995804*M0[2]*Gamma_inv_p1[3]; 
  } 
  else { 
    double Gamma_inv[8] = {0.0}; 
    ser_2x_p2_sqrt(Gamma2_inv, Gamma_inv); 
    n[0] = 0.5*Gamma_inv[7]*M0[7]+0.5*Gamma_inv[6]*M0[6]+0.5*Gamma_inv[5]*M0[5]+0.5*Gamma_inv[4]*M0[4]+0.5*Gamma_inv[3]*M0[3]+0.5*Gamma_inv[2]*M0[2]+0.5*Gamma_inv[1]*M0[1]+0.5*Gamma_inv[0]*M0[0]; 
    n[1] = 0.5000000000000001*Gamma_inv[5]*M0[7]+0.5000000000000001*M0[5]*Gamma_inv[7]+0.44721359549995804*Gamma_inv[3]*M0[6]+0.44721359549995804*M0[3]*Gamma_inv[6]+0.4472135954999579*Gamma_inv[1]*M0[4]+0.4472135954999579*M0[1]*Gamma_inv[4]+0.5*Gamma_inv[2]*M0[3]+0.5*M0[2]*Gamma_inv[3]+0.5*Gamma_inv[0]*M0[1]+0.5*M0[0]*Gamma_inv[1]; 
    n[2] = 0.44721359549995804*Gamma_inv[3]*M0[7]+0.44721359549995804*M0[3]*Gamma_inv[7]+0.5000000000000001*Gamma_inv[4]*M0[6]+0.5000000000000001*M0[4]*Gamma_inv[6]+0.4472135954999579*Gamma_inv[2]*M0[5]+0.4472135954999579*M0[2]*Gamma_inv[5]+0.5*Gamma_inv[1]*M0[3]+0.5*M0[1]*Gamma_inv[3]+0.5*Gamma_inv[0]*M0[2]+0.5*M0[0]*Gamma_inv[2]; 
    n[3] = 0.4*Gamma_inv[6]*M0[7]+0.44721359549995804*Gamma_inv[2]*M0[7]+0.4*M0[6]*Gamma_inv[7]+0.44721359549995804*M0[2]*Gamma_inv[7]+0.44721359549995804*Gamma_inv[1]*M0[6]+0.44721359549995804*M0[1]*Gamma_inv[6]+0.4472135954999579*Gamma_inv[3]*M0[5]+0.4472135954999579*M0[3]*Gamma_inv[5]+0.4472135954999579*Gamma_inv[3]*M0[4]+0.4472135954999579*M0[3]*Gamma_inv[4]+0.5*Gamma_inv[0]*M0[3]+0.5*M0[0]*Gamma_inv[3]+0.5*Gamma_inv[1]*M0[2]+0.5*M0[1]*Gamma_inv[2]; 
    n[4] = 0.4472135954999579*Gamma_inv[7]*M0[7]+0.31943828249996997*Gamma_inv[6]*M0[6]+0.5000000000000001*Gamma_inv[2]*M0[6]+0.5000000000000001*M0[2]*Gamma_inv[6]+0.31943828249996997*Gamma_inv[4]*M0[4]+0.5*Gamma_inv[0]*M0[4]+0.5*M0[0]*Gamma_inv[4]+0.4472135954999579*Gamma_inv[3]*M0[3]+0.4472135954999579*Gamma_inv[1]*M0[1]; 
    n[5] = 0.31943828249996997*Gamma_inv[7]*M0[7]+0.5000000000000001*Gamma_inv[1]*M0[7]+0.5000000000000001*M0[1]*Gamma_inv[7]+0.4472135954999579*Gamma_inv[6]*M0[6]+0.31943828249996997*Gamma_inv[5]*M0[5]+0.5*Gamma_inv[0]*M0[5]+0.5*M0[0]*Gamma_inv[5]+0.4472135954999579*Gamma_inv[3]*M0[3]+0.4472135954999579*Gamma_inv[2]*M0[2]; 
    n[6] = 0.4*Gamma_inv[3]*M0[7]+0.4*M0[3]*Gamma_inv[7]+0.4472135954999579*Gamma_inv[5]*M0[6]+0.31943828249996997*Gamma_inv[4]*M0[6]+0.5*Gamma_inv[0]*M0[6]+0.4472135954999579*M0[5]*Gamma_inv[6]+0.31943828249996997*M0[4]*Gamma_inv[6]+0.5*M0[0]*Gamma_inv[6]+0.5000000000000001*Gamma_inv[2]*M0[4]+0.5000000000000001*M0[2]*Gamma_inv[4]+0.44721359549995804*Gamma_inv[1]*M0[3]+0.44721359549995804*M0[1]*Gamma_inv[3]; 
    n[7] = 0.31943828249996997*Gamma_inv[5]*M0[7]+0.4472135954999579*Gamma_inv[4]*M0[7]+0.5*Gamma_inv[0]*M0[7]+0.31943828249996997*M0[5]*Gamma_inv[7]+0.4472135954999579*M0[4]*Gamma_inv[7]+0.5*M0[0]*Gamma_inv[7]+0.4*Gamma_inv[3]*M0[6]+0.4*M0[3]*Gamma_inv[6]+0.5000000000000001*Gamma_inv[1]*M0[5]+0.5000000000000001*M0[1]*Gamma_inv[5]+0.44721359549995804*Gamma_inv[2]*M0[3]+0.44721359549995804*M0[2]*Gamma_inv[3]; 
  } 

} 
 
