#include <gkyl_pkpm_remap_flux_kernels.h> 
GKYL_CU_DH void pkpm_remap_flux_trace_z_3x1v_ser_p1(double jdinv_l, double jdinv_r, const double *f_l, const double *f_r,
  double* GKYL_RESTRICT G_l, double* GKYL_RESTRICT G_r) 
{ 
  // jdinv_l/r: 1/J_z of the left/right cell (J_z: position-map Jacobian of direction z). 
  // f_l/f_r:   Jacobian-weighted blocks J [F, G, M0_perp, N_perp] of one velocity cell in the left/right cell. 
  // G_l/G_r:   Output modal traces of the eight blocks on the face of direction z: 
  //            G[c*12 + a*3 + m], a transverse mode, m Legendre mode in w. 
  for (int c = 0; c < 8; ++c) { 
    const double *fl = &f_l[c*24], *fr = &f_r[c*24]; 
    double *gl = &G_l[c*12], *gr = &G_r[c*12]; 
    gl[0] = jdinv_l*(1.224744871391589*fl[3]+0.7071067811865475*fl[0]); 
    gl[1] = jdinv_l*(1.224744871391589*fl[10]+0.7071067811865475*fl[4]); 
    gl[2] = jdinv_l*(1.224744871391589*fl[18]+0.7071067811865475*fl[11]); 
    gl[3] = jdinv_l*(1.224744871391589*fl[6]+0.7071067811865475*fl[1]); 
    gl[4] = jdinv_l*(1.224744871391589*fl[14]+0.7071067811865475*fl[8]); 
    gl[5] = jdinv_l*(1.224744871391589*fl[21]+0.7071067811865475*fl[16]); 
    gl[6] = jdinv_l*(1.224744871391589*fl[7]+0.7071067811865475*fl[2]); 
    gl[7] = jdinv_l*(1.224744871391589*fl[15]+0.7071067811865475*fl[9]); 
    gl[8] = jdinv_l*(1.224744871391589*fl[22]+0.7071067811865475*fl[17]); 
    gl[9] = jdinv_l*(1.224744871391589*fl[12]+0.7071067811865475*fl[5]); 
    gl[10] = jdinv_l*(1.224744871391589*fl[19]+0.7071067811865475*fl[13]); 
    gl[11] = jdinv_l*(1.224744871391589*fl[23]+0.7071067811865475*fl[20]); 
    gr[0] = jdinv_r*(0.7071067811865475*fr[0]-1.224744871391589*fr[3]); 
    gr[1] = jdinv_r*(0.7071067811865475*fr[4]-1.224744871391589*fr[10]); 
    gr[2] = jdinv_r*(0.7071067811865475*fr[11]-1.224744871391589*fr[18]); 
    gr[3] = jdinv_r*(0.7071067811865475*fr[1]-1.224744871391589*fr[6]); 
    gr[4] = jdinv_r*(0.7071067811865475*fr[8]-1.224744871391589*fr[14]); 
    gr[5] = jdinv_r*(0.7071067811865475*fr[16]-1.224744871391589*fr[21]); 
    gr[6] = jdinv_r*(0.7071067811865475*fr[2]-1.224744871391589*fr[7]); 
    gr[7] = jdinv_r*(0.7071067811865475*fr[9]-1.224744871391589*fr[15]); 
    gr[8] = jdinv_r*(0.7071067811865475*fr[17]-1.224744871391589*fr[22]); 
    gr[9] = jdinv_r*(0.7071067811865475*fr[5]-1.224744871391589*fr[12]); 
    gr[10] = jdinv_r*(0.7071067811865475*fr[13]-1.224744871391589*fr[19]); 
    gr[11] = jdinv_r*(0.7071067811865475*fr[20]-1.224744871391589*fr[23]); 
  } 
} 
