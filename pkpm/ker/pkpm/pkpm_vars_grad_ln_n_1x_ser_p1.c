#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH int pkpm_vars_grad_ln_n_1x_ser_p1(const double *rho, const double *grad_rho, double* GKYL_RESTRICT grad_ln_n) 
{ 
  // rho:       Input mass density. 
  // grad_rho:  Input gradient of [rho, 0, 0] (cdim x 3 expansions). 
  // grad_ln_n: Output [d_x ln n, ...] (cdim expansions): grad(rho)/rho at the 2 Gauss-Legendre nodes, 
  //            interpolated to the p=1 basis (the weak division at p=1). 
  // Returns the number of nodes where rho <= 0 (the quotient is zero there). 
 
  const double *grho_x = &grad_rho[0]; 
  double *gln_x = &grad_ln_n[0]; 
 
  double gln_n[1][2] = {{0.0}}; 
  int num_null = 0; 
  { 
    const double rho_q = 0.7071067811865475*rho[0]-0.7071067811865475*rho[1]; 
    if (rho_q > 0.0) { 
      const double inv_rho = 1.0/rho_q; 
      gln_n[0][0] = (0.7071067811865475*grho_x[0]-0.7071067811865475*grho_x[1])*inv_rho; 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double rho_q = 0.7071067811865475*rho[1]+0.7071067811865475*rho[0]; 
    if (rho_q > 0.0) { 
      const double inv_rho = 1.0/rho_q; 
      gln_n[0][1] = (0.7071067811865475*grho_x[1]+0.7071067811865475*grho_x[0])*inv_rho; 
    } else { 
      num_null += 1; 
    } 
  } 
 
  gln_x[0] = 0.7071067811865475*gln_n[0][1]+0.7071067811865475*gln_n[0][0]; 
  gln_x[1] = 0.7071067811865475*gln_n[0][1]-0.7071067811865475*gln_n[0][0]; 
 
  return num_null; 
} 
