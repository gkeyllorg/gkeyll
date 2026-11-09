#include <gkyl_pkpm_remap_flux_kernels.h> 
GKYL_CU_DH void pkpm_remap_flux_trace_x_3x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  // jdinv_l/r: 1/J_x of the left/right cell (J_x: position-map Jacobian of direction x). 
  // f_l/f_r:   Jacobian-weighted blocks J [F, G, M0_perp, N_perp] of one velocity cell in the left/right cell. 
  // G_l/G_r:   Output modal traces of the eight blocks on the face of direction x: 
  //            G[c*12 + a*3 + m], a transverse mode, m Legendre mode in w. 
  for (int c = 0; c < 8; ++c) { 
    const double *fl = &f_l[c*24], *fr = &f_r[c*24]; 
    double *gl = &G_l[c*12], *gr = &G_r[c*12]; 
    gl[0] = jdinv_l*(1.224744871391589*fl[1]+0.7071067811865475*fl[0]); 
    gl[1] = jdinv_l*(1.224744871391589*fl[8]+0.7071067811865475*fl[4]); 
    gl[2] = jdinv_l*(1.224744871391589*fl[16]+0.7071067811865475*fl[11]); 
    gl[3] = jdinv_l*(1.224744871391589*fl[5]+0.7071067811865475*fl[2]); 
    gl[4] = jdinv_l*(1.224744871391589*fl[13]+0.7071067811865475*fl[9]); 
    gl[5] = jdinv_l*(1.224744871391589*fl[20]+0.7071067811865475*fl[17]); 
    gl[6] = jdinv_l*(1.224744871391589*fl[6]+0.7071067811865475*fl[3]); 
    gl[7] = jdinv_l*(1.224744871391589*fl[14]+0.7071067811865475*fl[10]); 
    gl[8] = jdinv_l*(1.224744871391589*fl[21]+0.7071067811865475*fl[18]); 
    gl[9] = jdinv_l*(1.224744871391589*fl[12]+0.7071067811865475*fl[7]); 
    gl[10] = jdinv_l*(1.224744871391589*fl[19]+0.7071067811865475*fl[15]); 
    gl[11] = jdinv_l*(1.224744871391589*fl[23]+0.7071067811865475*fl[22]); 
    gr[0] = jdinv_r*(0.7071067811865475*fr[0]-1.224744871391589*fr[1]); 
    gr[1] = jdinv_r*(0.7071067811865475*fr[4]-1.224744871391589*fr[8]); 
    gr[2] = jdinv_r*(0.7071067811865475*fr[11]-1.224744871391589*fr[16]); 
    gr[3] = jdinv_r*(0.7071067811865475*fr[2]-1.224744871391589*fr[5]); 
    gr[4] = jdinv_r*(0.7071067811865475*fr[9]-1.224744871391589*fr[13]); 
    gr[5] = jdinv_r*(0.7071067811865475*fr[17]-1.224744871391589*fr[20]); 
    gr[6] = jdinv_r*(0.7071067811865475*fr[3]-1.224744871391589*fr[6]); 
    gr[7] = jdinv_r*(0.7071067811865475*fr[10]-1.224744871391589*fr[14]); 
    gr[8] = jdinv_r*(0.7071067811865475*fr[18]-1.224744871391589*fr[21]); 
    gr[9] = jdinv_r*(0.7071067811865475*fr[7]-1.224744871391589*fr[12]); 
    gr[10] = jdinv_r*(0.7071067811865475*fr[15]-1.224744871391589*fr[19]); 
    gr[11] = jdinv_r*(0.7071067811865475*fr[22]-1.224744871391589*fr[23]); 
  } 
} 
