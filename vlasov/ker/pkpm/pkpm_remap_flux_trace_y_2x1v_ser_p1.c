#include <gkyl_pkpm_remap_flux_kernels.h> 
GKYL_CU_DH void pkpm_remap_flux_trace_y_2x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  // jdinv_l/r: 1/J_y of the left/right cell (J_y: position-map Jacobian of direction y). 
  // f_l/f_r:   Jacobian-weighted blocks J [F, G, M0_perp, N_perp] of one velocity cell in the left/right cell. 
  // G_l/G_r:   Output modal traces of the eight blocks on the face of direction y: 
  //            G[c*6 + a*3 + m], a transverse mode, m Legendre mode in w. 
  for (int c = 0; c < 8; ++c) { 
    const double *fl = &f_l[c*12], *fr = &f_r[c*12]; 
    double *gl = &G_l[c*6], *gr = &G_r[c*6]; 
    gl[0] = jdinv_l*(1.224744871391589*fl[2]+0.7071067811865475*fl[0]); 
    gl[1] = jdinv_l*(1.224744871391589*fl[6]+0.7071067811865475*fl[3]); 
    gl[2] = jdinv_l*(1.224744871391589*fl[10]+0.7071067811865475*fl[7]); 
    gl[3] = jdinv_l*(1.224744871391589*fl[4]+0.7071067811865475*fl[1]); 
    gl[4] = jdinv_l*(1.224744871391589*fl[8]+0.7071067811865475*fl[5]); 
    gl[5] = jdinv_l*(1.224744871391589*fl[11]+0.7071067811865475*fl[9]); 
    gr[0] = jdinv_r*(0.7071067811865475*fr[0]-1.224744871391589*fr[2]); 
    gr[1] = jdinv_r*(0.7071067811865475*fr[3]-1.224744871391589*fr[6]); 
    gr[2] = jdinv_r*(0.7071067811865475*fr[7]-1.224744871391589*fr[10]); 
    gr[3] = jdinv_r*(0.7071067811865475*fr[1]-1.224744871391589*fr[4]); 
    gr[4] = jdinv_r*(0.7071067811865475*fr[5]-1.224744871391589*fr[8]); 
    gr[5] = jdinv_r*(0.7071067811865475*fr[9]-1.224744871391589*fr[11]); 
  } 
} 
