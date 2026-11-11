#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH void pkpm_vars_diag_1x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm,
  const double *bvar, const double *grad_b, const double *u, const double *grad_u, const double *theta, const double *alpha_0,
  const double *xi, const double *c_h, double* GKYL_RESTRICT fluid_io, double* GKYL_RESTRICT pkpm_vars_io) 
{ 
  // vlasov_pkpm_moms: Kinetic coupling moments [p_par, p_perp, q_par, q_perp]. 
  // euler_pkpm:       Fluid state [rho, rho ux, rho uy, rho uz, E]. 
  // bvar:             Magnetic field unit vector [bx, by, bz]. 
  // grad_b:           Weak gradient of b, [d_i bx, d_i by, d_i bz] per direction i. 
  // u:                Flow velocity [ux, uy, uz]. 
  // grad_u:           Gradient of u of the sources (in-cell, plus the blended face lift), [d_i ux, d_i uy, d_i uz] per direction i. 
  // theta:            T_perp/T = M0(G). 
  // alpha_0:          Pressure force a_P/xi closing the first moment. 
  // xi:               Thermal speed [xi, xi^2]. 
  // c_h:              Heating rate closing the scalar second moment. 
  // fluid_io:         Output [rho, M_i, P_ij + M_i u_j, E] (11 components). 
  // pkpm_vars_io:     Output [u_x, u_y, u_z, T_perp/m, a_P, div b, bb:grad u, T_perp/m div b, xi, c] (10 components). 
  const double *rho = &euler_pkpm[0]; 
  const double *ppar = &vlasov_pkpm_moms[0]; 
  const double *pperp = &vlasov_pkpm_moms[2]; 
  const double *bx = &bvar[0]; 
  const double *ux = &u[0]; 
  const double *Mx = &euler_pkpm[2]; 
  const double *by = &bvar[2]; 
  const double *uy = &u[2]; 
  const double *My = &euler_pkpm[4]; 
  const double *bz = &bvar[4]; 
  const double *uz = &u[4]; 
  const double *Mz = &euler_pkpm[6]; 
  const double *guxx = &grad_u[0]; 
  const double *guxy = &grad_u[2]; 
  const double *guxz = &grad_u[4]; 
  const double *gbxx = &grad_b[0]; 
  double *fluid_io_Sxx = &fluid_io[8]; 
  double *fluid_io_Sxy = &fluid_io[10]; 
  double *fluid_io_Sxz = &fluid_io[12]; 
  double *fluid_io_Syy = &fluid_io[14]; 
  double *fluid_io_Syz = &fluid_io[16]; 
  double *fluid_io_Szz = &fluid_io[18]; 
  double *pkpm_vars_io_tperp = &pkpm_vars_io[6]; 
  double *pkpm_vars_io_a_P = &pkpm_vars_io[8]; 
  double *pkpm_vars_io_div_b = &pkpm_vars_io[10]; 
  double *pkpm_vars_io_bb_grad_u = &pkpm_vars_io[12]; 
  double *pkpm_vars_io_theta_div_b = &pkpm_vars_io[14]; 
  const double *xi_1 = &xi[0]; 
  const double *xi_sq = &xi[2]; 
 
  for (int k = 0; k < 2; ++k) { 
    fluid_io[0 + k] = rho[k]; 
    fluid_io[2 + k] = Mx[k]; 
    fluid_io[4 + k] = My[k]; 
    fluid_io[6 + k] = Mz[k]; 
    fluid_io[20 + k] = euler_pkpm[8 + k]; 
    pkpm_vars_io[0 + k] = ux[k]; 
    pkpm_vars_io[2 + k] = uy[k]; 
    pkpm_vars_io[4 + k] = uz[k]; 
    pkpm_vars_io[16 + k] = xi_1[k]; 
    pkpm_vars_io[18 + k] = c_h[k]; 
  } 
 
  pkpm_vars_io_tperp[0] = 0.7071067811865475*theta[1]*xi_sq[1]+0.7071067811865475*theta[0]*xi_sq[0]; 
  pkpm_vars_io_tperp[1] = 0.7071067811865475*theta[0]*xi_sq[1]+0.7071067811865475*xi_sq[0]*theta[1]; 
  pkpm_vars_io_a_P[0] = 0.7071067811865475*alpha_0[1]*xi_1[1]+0.7071067811865475*alpha_0[0]*xi_1[0]; 
  pkpm_vars_io_a_P[1] = 0.7071067811865475*alpha_0[0]*xi_1[1]+0.7071067811865475*xi_1[0]*alpha_0[1]; 
  fluid_io_Sxx[0] = -(0.5*pperp[0]*pow(bx[1],2))+0.5*ppar[0]*pow(bx[1],2)-0.5*pperp[0]*pow(bx[0],2)+0.5*ppar[0]*pow(bx[0],2)+0.7071067811865475*Mx[1]*ux[1]-1.0*bx[0]*bx[1]*pperp[1]+bx[0]*bx[1]*ppar[1]+0.7071067811865475*Mx[0]*ux[0]+pperp[0]; 
  fluid_io_Sxx[1] = -(0.9*pperp[1]*pow(bx[1],2))+0.9*ppar[1]*pow(bx[1],2)-0.5*pperp[1]*pow(bx[0],2)+0.5*ppar[1]*pow(bx[0],2)+0.7071067811865475*Mx[0]*ux[1]+pperp[1]-1.0*bx[0]*pperp[0]*bx[1]+bx[0]*ppar[0]*bx[1]+0.7071067811865475*ux[0]*Mx[1]; 
  fluid_io_Sxy[0] = 0.7071067811865475*Mx[1]*uy[1]-0.5*bx[0]*by[1]*pperp[1]-0.5*by[0]*bx[1]*pperp[1]+0.5*bx[0]*by[1]*ppar[1]+0.5*by[0]*bx[1]*ppar[1]-0.5*pperp[0]*bx[1]*by[1]+0.5*ppar[0]*bx[1]*by[1]+0.7071067811865475*Mx[0]*uy[0]-0.5*bx[0]*by[0]*pperp[0]+0.5*bx[0]*by[0]*ppar[0]; 
  fluid_io_Sxy[1] = 0.7071067811865475*Mx[0]*uy[1]-0.9*bx[1]*by[1]*pperp[1]-0.5*bx[0]*by[0]*pperp[1]+0.9*bx[1]*by[1]*ppar[1]+0.5*bx[0]*by[0]*ppar[1]-0.5*bx[0]*pperp[0]*by[1]+0.5*bx[0]*ppar[0]*by[1]-0.5*by[0]*pperp[0]*bx[1]+0.5*by[0]*ppar[0]*bx[1]+0.7071067811865475*uy[0]*Mx[1]; 
  fluid_io_Sxz[0] = 0.7071067811865475*Mx[1]*uz[1]-0.5*bx[0]*bz[1]*pperp[1]-0.5*bz[0]*bx[1]*pperp[1]+0.5*bx[0]*bz[1]*ppar[1]+0.5*bz[0]*bx[1]*ppar[1]-0.5*pperp[0]*bx[1]*bz[1]+0.5*ppar[0]*bx[1]*bz[1]+0.7071067811865475*Mx[0]*uz[0]-0.5*bx[0]*bz[0]*pperp[0]+0.5*bx[0]*bz[0]*ppar[0]; 
  fluid_io_Sxz[1] = 0.7071067811865475*Mx[0]*uz[1]-0.9*bx[1]*bz[1]*pperp[1]-0.5*bx[0]*bz[0]*pperp[1]+0.9*bx[1]*bz[1]*ppar[1]+0.5*bx[0]*bz[0]*ppar[1]-0.5*bx[0]*pperp[0]*bz[1]+0.5*bx[0]*ppar[0]*bz[1]-0.5*bz[0]*pperp[0]*bx[1]+0.5*bz[0]*ppar[0]*bx[1]+0.7071067811865475*uz[0]*Mx[1]; 
  fluid_io_Syy[0] = -(0.5*pperp[0]*pow(by[1],2))+0.5*ppar[0]*pow(by[1],2)-0.5*pperp[0]*pow(by[0],2)+0.5*ppar[0]*pow(by[0],2)+0.7071067811865475*My[1]*uy[1]-1.0*by[0]*by[1]*pperp[1]+by[0]*by[1]*ppar[1]+0.7071067811865475*My[0]*uy[0]+pperp[0]; 
  fluid_io_Syy[1] = -(0.9*pperp[1]*pow(by[1],2))+0.9*ppar[1]*pow(by[1],2)-0.5*pperp[1]*pow(by[0],2)+0.5*ppar[1]*pow(by[0],2)+0.7071067811865475*My[0]*uy[1]+pperp[1]-1.0*by[0]*pperp[0]*by[1]+by[0]*ppar[0]*by[1]+0.7071067811865475*uy[0]*My[1]; 
  fluid_io_Syz[0] = 0.7071067811865475*My[1]*uz[1]-0.5*by[0]*bz[1]*pperp[1]-0.5*bz[0]*by[1]*pperp[1]+0.5*by[0]*bz[1]*ppar[1]+0.5*bz[0]*by[1]*ppar[1]-0.5*pperp[0]*by[1]*bz[1]+0.5*ppar[0]*by[1]*bz[1]+0.7071067811865475*My[0]*uz[0]-0.5*by[0]*bz[0]*pperp[0]+0.5*by[0]*bz[0]*ppar[0]; 
  fluid_io_Syz[1] = 0.7071067811865475*My[0]*uz[1]-0.9*by[1]*bz[1]*pperp[1]-0.5*by[0]*bz[0]*pperp[1]+0.9*by[1]*bz[1]*ppar[1]+0.5*by[0]*bz[0]*ppar[1]-0.5*by[0]*pperp[0]*bz[1]+0.5*by[0]*ppar[0]*bz[1]-0.5*bz[0]*pperp[0]*by[1]+0.5*bz[0]*ppar[0]*by[1]+0.7071067811865475*uz[0]*My[1]; 
  fluid_io_Szz[0] = -(0.5*pperp[0]*pow(bz[1],2))+0.5*ppar[0]*pow(bz[1],2)-0.5*pperp[0]*pow(bz[0],2)+0.5*ppar[0]*pow(bz[0],2)+0.7071067811865475*Mz[1]*uz[1]-1.0*bz[0]*bz[1]*pperp[1]+bz[0]*bz[1]*ppar[1]+0.7071067811865475*Mz[0]*uz[0]+pperp[0]; 
  fluid_io_Szz[1] = -(0.9*pperp[1]*pow(bz[1],2))+0.9*ppar[1]*pow(bz[1],2)-0.5*pperp[1]*pow(bz[0],2)+0.5*ppar[1]*pow(bz[0],2)+0.7071067811865475*Mz[0]*uz[1]+pperp[1]-1.0*bz[0]*pperp[0]*bz[1]+bz[0]*ppar[0]*bz[1]+0.7071067811865475*uz[0]*Mz[1]; 
 
  pkpm_vars_io_div_b[0] = gbxx[0]; 
  pkpm_vars_io_div_b[1] = gbxx[1]; 
  pkpm_vars_io_bb_grad_u[0] = 0.5*guxx[0]*pow(bx[1],2)+0.5*guxx[0]*pow(bx[0],2)+0.5*bx[0]*bz[1]*guxz[1]+0.5*bz[0]*bx[1]*guxz[1]+0.5*bx[0]*by[1]*guxy[1]+0.5*by[0]*bx[1]*guxy[1]+bx[0]*bx[1]*guxx[1]+0.5*guxz[0]*bx[1]*bz[1]+0.5*guxy[0]*bx[1]*by[1]+0.5*bx[0]*bz[0]*guxz[0]+0.5*bx[0]*by[0]*guxy[0]; 
  pkpm_vars_io_bb_grad_u[1] = 0.9*guxx[1]*pow(bx[1],2)+0.5*guxx[1]*pow(bx[0],2)+0.9*bx[1]*bz[1]*guxz[1]+0.5*bx[0]*bz[0]*guxz[1]+0.9*bx[1]*by[1]*guxy[1]+0.5*bx[0]*by[0]*guxy[1]+0.5*bx[0]*guxz[0]*bz[1]+0.5*bx[0]*guxy[0]*by[1]+0.5*bz[0]*guxz[0]*bx[1]+0.5*by[0]*guxy[0]*bx[1]+bx[0]*guxx[0]*bx[1]; 
  pkpm_vars_io_theta_div_b[0] = 0.5*gbxx[0]*theta[1]*xi_sq[1]+0.5*theta[0]*gbxx[1]*xi_sq[1]+0.5*xi_sq[0]*gbxx[1]*theta[1]+0.5*gbxx[0]*theta[0]*xi_sq[0]; 
  pkpm_vars_io_theta_div_b[1] = 0.9*gbxx[1]*theta[1]*xi_sq[1]+0.5*gbxx[0]*theta[0]*xi_sq[1]+0.5*gbxx[0]*xi_sq[0]*theta[1]+0.5*theta[0]*xi_sq[0]*gbxx[1]; 
} 
