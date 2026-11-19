#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH void pkpm_vars_int_2x_ser_p1(const double *fluid, const double *vlasov_pkpm_moms, const double *u, double* GKYL_RESTRICT int_pkpm_vars) 
{ 
  // fluid:            Fluid state [rho, rho ux, rho uy, rho uz, E]. 
  // vlasov_pkpm_moms: Kinetic coupling moments [p_par, p_perp, q_par, q_perp, s (x, y, z), r (x, y, z)]; p_par and p_perp are read. 
  // u:                Flow velocity [ux, uy, uz]. 
  // int_pkpm_vars:    Incremented cell integrals [rho, rho ux, rho uy, rho uz, rho ux^2, rho uy^2, rho uz^2, p_par, p_perp]. 
 
  const double *rho = &fluid[0]; 
  const double *M0 = &fluid[4]; 
  const double *u0 = &u[0]; 
  const double *M1 = &fluid[8]; 
  const double *u1 = &u[4]; 
  const double *M2 = &fluid[12]; 
  const double *u2 = &u[8]; 
  const double *ppar = &vlasov_pkpm_moms[0]; 
  const double *pperp = &vlasov_pkpm_moms[4]; 
 
  int_pkpm_vars[0] += 0.5*rho[0]; 
  int_pkpm_vars[1] += 0.5*M0[0]; 
  int_pkpm_vars[2] += 0.5*M1[0]; 
  int_pkpm_vars[3] += 0.5*M2[0]; 
  int_pkpm_vars[4] += 0.25*M0[3]*u0[3]+0.25*M0[2]*u0[2]+0.25*M0[1]*u0[1]+0.25*M0[0]*u0[0]; 
  int_pkpm_vars[5] += 0.25*M1[3]*u1[3]+0.25*M1[2]*u1[2]+0.25*M1[1]*u1[1]+0.25*M1[0]*u1[0]; 
  int_pkpm_vars[6] += 0.25*M2[3]*u2[3]+0.25*M2[2]*u2[2]+0.25*M2[1]*u2[1]+0.25*M2[0]*u2[0]; 
  int_pkpm_vars[7] += 0.5*ppar[0]; 
  int_pkpm_vars[8] += 0.5*pperp[0]; 
} 
