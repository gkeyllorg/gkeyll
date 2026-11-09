#include <gkyl_pkpm_remap_flux_kernels.h> 
GKYL_CU_DH void pkpm_remap_flux_trace_x_1x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  // jdinv_l/r: 1/J_x of the left/right cell (J_x: position-map Jacobian of direction x). 
  // f_l/f_r:   Jacobian-weighted blocks J [F, G, M0_perp, N_perp] of one velocity cell in the left/right cell. 
  // G_l/G_r:   Output modal traces of the eight blocks on the face of direction x: 
  //            G[c*3 + a*3 + m], a transverse mode, m Legendre mode in w. 
  for (int c = 0; c < 8; ++c) { 
    const double *fl = &f_l[c*6], *fr = &f_r[c*6]; 
    double *gl = &G_l[c*3], *gr = &G_r[c*3]; 
    gl[0] = jdinv_l*(1.224744871391589*fl[1]+0.7071067811865475*fl[0]); 
    gl[1] = jdinv_l*(1.224744871391589*fl[3]+0.7071067811865475*fl[2]); 
    gl[2] = jdinv_l*(1.224744871391589*fl[5]+0.7071067811865475*fl[4]); 
    gr[0] = jdinv_r*(0.7071067811865475*fr[0]-1.224744871391589*fr[1]); 
    gr[1] = jdinv_r*(0.7071067811865475*fr[2]-1.224744871391589*fr[3]); 
    gr[2] = jdinv_r*(0.7071067811865475*fr[4]-1.224744871391589*fr[5]); 
  } 
} 
