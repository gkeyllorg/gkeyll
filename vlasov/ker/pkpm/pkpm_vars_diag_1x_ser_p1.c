#include <gkyl_pkpm_vars_kernels.h> 
GKYL_CU_DH void pkpm_vars_diag_1x_ser_p1(const double *vlasov_pkpm_moms, const double *bvar, const double *grad_b,
  const double *u, const double *grad_u, const double *theta, const double *alpha_0, const double *xi,
  const double *c_h, double* GKYL_RESTRICT pkpm_prim_io, double* GKYL_RESTRICT pkpm_vars_io) 
{ 
  // vlasov_pkpm_moms: Kinetic coupling moments [p_par, p_perp, q_par, q_perp, s (x, y, z), r (x, y, z)]: 
  //                   s the parallel-perpendicular stress, r the perpendicular heat flux. 
  // bvar:             Magnetic field unit vector [bx, by, bz]. 
  // grad_b:           Weak gradient of b, [d_i bx, d_i by, d_i bz] per direction i. 
  // u:                Flow velocity [ux, uy, uz]. 
  // grad_u:           Gradient of u of the sources (in-cell, plus the blended face lift), [d_i ux, d_i uy, d_i uz] per direction i. 
  // theta:            T_perp/T = M0(G). 
  // alpha_0:          Pressure force a_P/xi closing the first moment. 
  // xi:               Thermal speed and pressure [xi, xi^2, p]. 
  // c_h:              Heating rate closing the scalar second moment. 
  // pkpm_prim_io:     Output [P_xx, P_xy, P_xz, P_yy, P_yz, P_zz, q_x, q_y, q_z] (9 components): the pressure 
  //                   tensor and the heat flux of the fluid equations. 
  // pkpm_vars_io:     Output [u_x, u_y, u_z, T_perp/m, a_P, div b, bb:grad u, T_perp/m div b, xi, c] (10 components). 
  const double *ppar = &vlasov_pkpm_moms[0]; 
  const double *pperp = &vlasov_pkpm_moms[2]; 
  const double *qpar = &vlasov_pkpm_moms[4]; 
  const double *qperp = &vlasov_pkpm_moms[6]; 
  const double *sx = &vlasov_pkpm_moms[8]; 
  const double *rx = &vlasov_pkpm_moms[14]; 
  const double *bx = &bvar[0]; 
  const double *ux = &u[0]; 
  const double *sy = &vlasov_pkpm_moms[10]; 
  const double *ry = &vlasov_pkpm_moms[16]; 
  const double *by = &bvar[2]; 
  const double *uy = &u[2]; 
  const double *sz = &vlasov_pkpm_moms[12]; 
  const double *rz = &vlasov_pkpm_moms[18]; 
  const double *bz = &bvar[4]; 
  const double *uz = &u[4]; 
  const double *guxx = &grad_u[0]; 
  const double *guxy = &grad_u[2]; 
  const double *guxz = &grad_u[4]; 
  const double *gbxx = &grad_b[0]; 
  double *pkpm_prim_io_Pxx = &pkpm_prim_io[0]; 
  double *pkpm_prim_io_Pxy = &pkpm_prim_io[2]; 
  double *pkpm_prim_io_Pxz = &pkpm_prim_io[4]; 
  double *pkpm_prim_io_Pyy = &pkpm_prim_io[6]; 
  double *pkpm_prim_io_Pyz = &pkpm_prim_io[8]; 
  double *pkpm_prim_io_Pzz = &pkpm_prim_io[10]; 
  double *pkpm_prim_io_qx = &pkpm_prim_io[12]; 
  double *pkpm_prim_io_qy = &pkpm_prim_io[14]; 
  double *pkpm_prim_io_qz = &pkpm_prim_io[16]; 
  double *pkpm_vars_io_tperp = &pkpm_vars_io[6]; 
  double *pkpm_vars_io_a_P = &pkpm_vars_io[8]; 
  double *pkpm_vars_io_div_b = &pkpm_vars_io[10]; 
  double *pkpm_vars_io_bb_grad_u = &pkpm_vars_io[12]; 
  double *pkpm_vars_io_theta_div_b = &pkpm_vars_io[14]; 
  const double *xi_1 = &xi[0]; 
  const double *xi_sq = &xi[2]; 
 
  for (int k = 0; k < 2; ++k) { 
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
  pkpm_prim_io_Pxx[0] = -(0.5*pperp[0]*pow(bx[1],2))+0.5*ppar[0]*pow(bx[1],2)-0.5*pperp[0]*pow(bx[0],2)+0.5*ppar[0]*pow(bx[0],2)+1.4142135623730951*bx[1]*sx[1]-1.0*bx[0]*bx[1]*pperp[1]+bx[0]*bx[1]*ppar[1]+1.4142135623730951*bx[0]*sx[0]+pperp[0]; 
  pkpm_prim_io_Pxx[1] = -(0.9*pperp[1]*pow(bx[1],2))+0.9*ppar[1]*pow(bx[1],2)-0.5*pperp[1]*pow(bx[0],2)+0.5*ppar[1]*pow(bx[0],2)+1.4142135623730951*bx[0]*sx[1]+pperp[1]+1.4142135623730951*sx[0]*bx[1]-1.0*bx[0]*pperp[0]*bx[1]+bx[0]*ppar[0]*bx[1]; 
  pkpm_prim_io_Pxy[0] = 0.7071067811865475*bx[1]*sy[1]+0.7071067811865475*by[1]*sx[1]-0.5*bx[0]*by[1]*pperp[1]-0.5*by[0]*bx[1]*pperp[1]+0.5*bx[0]*by[1]*ppar[1]+0.5*by[0]*bx[1]*ppar[1]-0.5*pperp[0]*bx[1]*by[1]+0.5*ppar[0]*bx[1]*by[1]+0.7071067811865475*bx[0]*sy[0]+0.7071067811865475*by[0]*sx[0]-0.5*bx[0]*by[0]*pperp[0]+0.5*bx[0]*by[0]*ppar[0]; 
  pkpm_prim_io_Pxy[1] = 0.7071067811865475*bx[0]*sy[1]+0.7071067811865475*by[0]*sx[1]-0.9*bx[1]*by[1]*pperp[1]-0.5*bx[0]*by[0]*pperp[1]+0.9*bx[1]*by[1]*ppar[1]+0.5*bx[0]*by[0]*ppar[1]+0.7071067811865475*sx[0]*by[1]-0.5*bx[0]*pperp[0]*by[1]+0.5*bx[0]*ppar[0]*by[1]+0.7071067811865475*sy[0]*bx[1]-0.5*by[0]*pperp[0]*bx[1]+0.5*by[0]*ppar[0]*bx[1]; 
  pkpm_prim_io_Pxz[0] = 0.7071067811865475*bx[1]*sz[1]+0.7071067811865475*bz[1]*sx[1]-0.5*bx[0]*bz[1]*pperp[1]-0.5*bz[0]*bx[1]*pperp[1]+0.5*bx[0]*bz[1]*ppar[1]+0.5*bz[0]*bx[1]*ppar[1]-0.5*pperp[0]*bx[1]*bz[1]+0.5*ppar[0]*bx[1]*bz[1]+0.7071067811865475*bx[0]*sz[0]+0.7071067811865475*bz[0]*sx[0]-0.5*bx[0]*bz[0]*pperp[0]+0.5*bx[0]*bz[0]*ppar[0]; 
  pkpm_prim_io_Pxz[1] = 0.7071067811865475*bx[0]*sz[1]+0.7071067811865475*bz[0]*sx[1]-0.9*bx[1]*bz[1]*pperp[1]-0.5*bx[0]*bz[0]*pperp[1]+0.9*bx[1]*bz[1]*ppar[1]+0.5*bx[0]*bz[0]*ppar[1]+0.7071067811865475*sx[0]*bz[1]-0.5*bx[0]*pperp[0]*bz[1]+0.5*bx[0]*ppar[0]*bz[1]+0.7071067811865475*sz[0]*bx[1]-0.5*bz[0]*pperp[0]*bx[1]+0.5*bz[0]*ppar[0]*bx[1]; 
  pkpm_prim_io_Pyy[0] = -(0.5*pperp[0]*pow(by[1],2))+0.5*ppar[0]*pow(by[1],2)-0.5*pperp[0]*pow(by[0],2)+0.5*ppar[0]*pow(by[0],2)+1.4142135623730951*by[1]*sy[1]-1.0*by[0]*by[1]*pperp[1]+by[0]*by[1]*ppar[1]+1.4142135623730951*by[0]*sy[0]+pperp[0]; 
  pkpm_prim_io_Pyy[1] = -(0.9*pperp[1]*pow(by[1],2))+0.9*ppar[1]*pow(by[1],2)-0.5*pperp[1]*pow(by[0],2)+0.5*ppar[1]*pow(by[0],2)+1.4142135623730951*by[0]*sy[1]+pperp[1]+1.4142135623730951*sy[0]*by[1]-1.0*by[0]*pperp[0]*by[1]+by[0]*ppar[0]*by[1]; 
  pkpm_prim_io_Pyz[0] = 0.7071067811865475*by[1]*sz[1]+0.7071067811865475*bz[1]*sy[1]-0.5*by[0]*bz[1]*pperp[1]-0.5*bz[0]*by[1]*pperp[1]+0.5*by[0]*bz[1]*ppar[1]+0.5*bz[0]*by[1]*ppar[1]-0.5*pperp[0]*by[1]*bz[1]+0.5*ppar[0]*by[1]*bz[1]+0.7071067811865475*by[0]*sz[0]+0.7071067811865475*bz[0]*sy[0]-0.5*by[0]*bz[0]*pperp[0]+0.5*by[0]*bz[0]*ppar[0]; 
  pkpm_prim_io_Pyz[1] = 0.7071067811865475*by[0]*sz[1]+0.7071067811865475*bz[0]*sy[1]-0.9*by[1]*bz[1]*pperp[1]-0.5*by[0]*bz[0]*pperp[1]+0.9*by[1]*bz[1]*ppar[1]+0.5*by[0]*bz[0]*ppar[1]+0.7071067811865475*sy[0]*bz[1]-0.5*by[0]*pperp[0]*bz[1]+0.5*by[0]*ppar[0]*bz[1]+0.7071067811865475*sz[0]*by[1]-0.5*bz[0]*pperp[0]*by[1]+0.5*bz[0]*ppar[0]*by[1]; 
  pkpm_prim_io_Pzz[0] = -(0.5*pperp[0]*pow(bz[1],2))+0.5*ppar[0]*pow(bz[1],2)-0.5*pperp[0]*pow(bz[0],2)+0.5*ppar[0]*pow(bz[0],2)+1.4142135623730951*bz[1]*sz[1]-1.0*bz[0]*bz[1]*pperp[1]+bz[0]*bz[1]*ppar[1]+1.4142135623730951*bz[0]*sz[0]+pperp[0]; 
  pkpm_prim_io_Pzz[1] = -(0.9*pperp[1]*pow(bz[1],2))+0.9*ppar[1]*pow(bz[1],2)-0.5*pperp[1]*pow(bz[0],2)+0.5*ppar[1]*pow(bz[0],2)+1.4142135623730951*bz[0]*sz[1]+pperp[1]+1.4142135623730951*sz[0]*bz[1]-1.0*bz[0]*pperp[0]*bz[1]+bz[0]*ppar[0]*bz[1]; 
  pkpm_prim_io_qx[0] = 0.7071067811865475*bx[1]*qperp[1]+0.3535533905932737*bx[1]*qpar[1]+rx[0]+0.7071067811865475*bx[0]*qperp[0]+0.3535533905932737*bx[0]*qpar[0]; 
  pkpm_prim_io_qx[1] = rx[1]+0.7071067811865475*bx[0]*qperp[1]+0.3535533905932737*bx[0]*qpar[1]+0.7071067811865475*qperp[0]*bx[1]+0.3535533905932737*qpar[0]*bx[1]; 
  pkpm_prim_io_qy[0] = 0.7071067811865475*by[1]*qperp[1]+0.3535533905932737*by[1]*qpar[1]+ry[0]+0.7071067811865475*by[0]*qperp[0]+0.3535533905932737*by[0]*qpar[0]; 
  pkpm_prim_io_qy[1] = ry[1]+0.7071067811865475*by[0]*qperp[1]+0.3535533905932737*by[0]*qpar[1]+0.7071067811865475*qperp[0]*by[1]+0.3535533905932737*qpar[0]*by[1]; 
  pkpm_prim_io_qz[0] = 0.7071067811865475*bz[1]*qperp[1]+0.3535533905932737*bz[1]*qpar[1]+rz[0]+0.7071067811865475*bz[0]*qperp[0]+0.3535533905932737*bz[0]*qpar[0]; 
  pkpm_prim_io_qz[1] = rz[1]+0.7071067811865475*bz[0]*qperp[1]+0.3535533905932737*bz[0]*qpar[1]+0.7071067811865475*qperp[0]*bz[1]+0.3535533905932737*qpar[0]*bz[1]; 
 
  pkpm_vars_io_div_b[0] = gbxx[0]; 
  pkpm_vars_io_div_b[1] = gbxx[1]; 
  pkpm_vars_io_bb_grad_u[0] = 0.5*guxx[0]*pow(bx[1],2)+0.5*guxx[0]*pow(bx[0],2)+0.5*bx[0]*bz[1]*guxz[1]+0.5*bz[0]*bx[1]*guxz[1]+0.5*bx[0]*by[1]*guxy[1]+0.5*by[0]*bx[1]*guxy[1]+bx[0]*bx[1]*guxx[1]+0.5*guxz[0]*bx[1]*bz[1]+0.5*guxy[0]*bx[1]*by[1]+0.5*bx[0]*bz[0]*guxz[0]+0.5*bx[0]*by[0]*guxy[0]; 
  pkpm_vars_io_bb_grad_u[1] = 0.9*guxx[1]*pow(bx[1],2)+0.5*guxx[1]*pow(bx[0],2)+0.9*bx[1]*bz[1]*guxz[1]+0.5*bx[0]*bz[0]*guxz[1]+0.9*bx[1]*by[1]*guxy[1]+0.5*bx[0]*by[0]*guxy[1]+0.5*bx[0]*guxz[0]*bz[1]+0.5*bx[0]*guxy[0]*by[1]+0.5*bz[0]*guxz[0]*bx[1]+0.5*by[0]*guxy[0]*bx[1]+bx[0]*guxx[0]*bx[1]; 
  pkpm_vars_io_theta_div_b[0] = 0.5*gbxx[0]*theta[1]*xi_sq[1]+0.5*theta[0]*gbxx[1]*xi_sq[1]+0.5*xi_sq[0]*gbxx[1]*theta[1]+0.5*gbxx[0]*theta[0]*xi_sq[0]; 
  pkpm_vars_io_theta_div_b[1] = 0.9*gbxx[1]*theta[1]*xi_sq[1]+0.5*gbxx[0]*theta[0]*xi_sq[1]+0.5*gbxx[0]*xi_sq[0]*theta[1]+0.5*theta[0]*xi_sq[0]*gbxx[1]; 
} 
