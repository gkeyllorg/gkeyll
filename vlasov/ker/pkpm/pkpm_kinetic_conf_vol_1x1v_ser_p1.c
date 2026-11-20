#include <gkyl_pkpm_kinetic_kernels.h> 

// Streaming volume term of one block X: out += integral( d_x phi . (xi b w + u) X ), with xb[d] the 
// exact degree <= 2 expansion of xi b_d (k0 + 3 k1 + 9 k2). 
GKYL_CU_DH static inline void pkpm_kinetic_conf_vol_1x1v_ser_p1_block(const double *dx1, double wvpar, double dvpar, 
  const double *xb, const double *u, const double *X, double* GKYL_RESTRICT out) 
{ 
  const double dx10 = dx1[0]; 
  const double *xbx = &xb[0]; 
  const double *ux = &u[0]; 
  const double *uy = &u[2]; 
  const double *uz = &u[4]; 
  out[1] += (1.224744871391589*X[1]*xbx[1]+1.224744871391589*X[0]*xbx[0])*dx10*wvpar+(0.3535533905932737*xbx[1]*X[3]+0.3535533905932737*xbx[0]*X[2])*dvpar*dx10+(1.224744871391589*X[1]*ux[1]+1.224744871391589*X[0]*ux[0])*dx10; 
  out[3] += (1.224744871391589*xbx[1]*X[3]+1.224744871391589*xbx[0]*X[2])*dx10*wvpar+(0.31622776601683794*xbx[1]*X[5]+0.3162277660168379*xbx[0]*X[4]+0.3535533905932737*X[1]*xbx[1]+0.3535533905932737*X[0]*xbx[0])*dvpar*dx10+(1.224744871391589*ux[1]*X[3]+1.224744871391589*ux[0]*X[2])*dx10; 
  out[5] += (1.224744871391589*xbx[1]*X[5]+1.224744871391589*xbx[0]*X[4])*dx10*wvpar+(0.31622776601683794*xbx[1]*X[3]+0.31622776601683794*xbx[0]*X[2])*dvpar*dx10+(1.224744871391589*ux[1]*X[5]+1.224744871391589*ux[0]*X[4])*dx10; 
} 

GKYL_CU_DH double pkpm_kinetic_conf_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *jacob_pos,
  const double *bvar, const double *u, const double *xi, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:      Cell-center coordinates (computational). 
  // dxv[NDIM]:    Cell spacing (computational). 
  // jacob_pos:    Position-map Jacobian per configuration direction (cdim blocks of 2 equal 
  //               entries); the streaming divergence in direction d carries 1/J_d. 
  // bvar:         Magnetic field unit vector [bx, by, bz]. 
  // u:            Flow velocity [ux, uy, uz]. 
  // xi:           Thermal speed [xi, xi^2, p] (the first expansion is used). 
  // f:            Input Jacobian-weighted normalized distributions J [F, G, M0_perp, N_perp]. 
  // out:          Incremented output for J [F, G, M0_perp, N_perp]. 
  // Returns 0 (the streaming speeds are accumulated by the surface-flux updaters). 
  double dx1[1]; 
  dx1[0] = 2.0/(dxv[0]*jacob_pos[0]); 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[2]; 
  const double *bz = &bvar[4]; 
  const double *xi_1 = &xi[0]; 
  // xi b_d in the tensor basis of degree <= 2 (exact) 
  double xb[3]; 
  xb[0] = 0.7071067811865475*bx[1]*xi_1[1]+0.7071067811865475*bx[0]*xi_1[0]; 
  xb[1] = 0.7071067811865475*bx[0]*xi_1[1]+0.7071067811865475*xi_1[0]*bx[1]; 
  xb[2] = 0.6324555320336759*bx[1]*xi_1[1]; 
  for (int c = 0; c < 8; ++c) { 
    pkpm_kinetic_conf_vol_1x1v_ser_p1_block(dx1, w[1], dxv[1], xb, u, &f[c*6], &out[c*6]); 
  } 
  return 0.0; 
} 
