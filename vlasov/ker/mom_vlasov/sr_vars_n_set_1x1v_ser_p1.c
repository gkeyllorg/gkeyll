#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
GKYL_CU_DH void sr_vars_n_set_1x1v_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A:     preallocated LHS matrix. 
  // rhs:   preallocated RHS vector. 
  // M0:    Lab frame density = Gamma*n.
  // M1i:   Lab frame flux = Gamma*n*V_drift_i.

  double M0_inv[2] = {0.0}; 
  ser_1x_p1_inv(M0, M0_inv); 
  double V_drift[2] = {0.0}; 
  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs0 = gkyl_nmat_get(rhs, count+0); 
  gkyl_mat_clear(&rhs0, 0.0); 
  const double *M10 = &M1i[0]; 
  binop_mul_1d_ser_p1(M0_inv, M10, V_drift); 
  gkyl_mat_set(&rhs0,0,0,V_drift[0]); 
  gkyl_mat_set(&rhs0,1,0,V_drift[1]); 

} 
