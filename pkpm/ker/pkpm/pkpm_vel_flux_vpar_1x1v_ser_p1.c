#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 

// The p=1 field f (gkyl basis) at the 5 Gauss-Legendre nodes (5 per direction, row-major: the 
// last direction fastest), direction by direction. 
GKYL_CU_DH static inline void pkpm_vel_flux_vpar_1x1v_ser_p1_p1(const double *f, double* GKYL_RESTRICT fn) 
{ 
  fn[0] = 0.7071067811865475*f[0]-1.1098391188717989*f[1]; 
  fn[1] = 0.7071067811865475*f[0]-0.6594875259537023*f[1]; 
  fn[2] = 0.7071067811865475*f[0]; 
  fn[3] = 0.6594875259537023*f[1]+0.7071067811865475*f[0]; 
  fn[4] = 1.1098391188717989*f[1]+0.7071067811865475*f[0]; 
} 

// Projection of node values onto the stored face format of the p=1 surface basis (gkyl order). 
GKYL_CU_DH static inline void pkpm_vel_flux_vpar_1x1v_ser_p1_prj(const double *v, double* GKYL_RESTRICT out) 
{ 
  out[0] = 0.23692688505618908*v[4]+0.47862867049936647*v[3]+0.5688888888888889*v[2]+0.47862867049936647*v[1]+0.23692688505618908*v[0]; 
  out[1] = 0.5367459204973624*v[4]+0.6443171250014853*v[3]-0.6443171250014853*v[1]-0.5367459204973624*v[0]; 
} 

GKYL_CU_DH void pkpm_vel_flux_vpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nv,
  const double *f_l, const double *f_r, double* GKYL_RESTRICT flux) 
{ 
  // w[NDIM]:   Cell-center coordinates (current cell). 
  // dxv[NDIM]: Cell spacing (current cell). 
  // nv:        Node values [5 nodes][PKPM_HS_NODE_NV] of the w flux of the configuration cell. 
  // f_l/f_r:   Distributions [F, G, M0_perp, N_perp] in the cell below/current cell. 
  // flux:      Output surface fluxes [8 blocks] on the lower w face (stored face format). 
  // (The time-step rate of the w fluxes is the column rate of the closure kernel, plan 17.27.) 
  const double v_face = w[1] - 0.5*dxv[1]; 
  double tl[2], tr[2], ul[8][5], ur[8][5]; 
  tl[0] = 1.5811388300841895*f_l[4]+1.224744871391589*f_l[2]+0.7071067811865475*f_l[0]; 
  tl[1] = 1.5811388300841895*f_l[5]+1.224744871391589*f_l[3]+0.7071067811865475*f_l[1]; 
  tr[0] = 1.5811388300841895*f_r[4]-1.224744871391589*f_r[2]+0.7071067811865475*f_r[0]; 
  tr[1] = 1.5811388300841895*f_r[5]-1.224744871391589*f_r[3]+0.7071067811865475*f_r[1]; 
  pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tl, ul[0]); pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tr, ur[0]); 
  tl[0] = 1.5811388300841895*f_l[10]+1.224744871391589*f_l[8]+0.7071067811865475*f_l[6]; 
  tl[1] = 1.5811388300841895*f_l[11]+1.224744871391589*f_l[9]+0.7071067811865475*f_l[7]; 
  tr[0] = 1.5811388300841895*f_r[10]-1.224744871391589*f_r[8]+0.7071067811865475*f_r[6]; 
  tr[1] = 1.5811388300841895*f_r[11]-1.224744871391589*f_r[9]+0.7071067811865475*f_r[7]; 
  pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tl, ul[1]); pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tr, ur[1]); 
  tl[0] = 1.5811388300841895*f_l[16]+1.224744871391589*f_l[14]+0.7071067811865475*f_l[12]; 
  tl[1] = 1.5811388300841895*f_l[17]+1.224744871391589*f_l[15]+0.7071067811865475*f_l[13]; 
  tr[0] = 1.5811388300841895*f_r[16]-1.224744871391589*f_r[14]+0.7071067811865475*f_r[12]; 
  tr[1] = 1.5811388300841895*f_r[17]-1.224744871391589*f_r[15]+0.7071067811865475*f_r[13]; 
  pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tl, ul[2]); pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tr, ur[2]); 
  tl[0] = 1.5811388300841895*f_l[22]+1.224744871391589*f_l[20]+0.7071067811865475*f_l[18]; 
  tl[1] = 1.5811388300841895*f_l[23]+1.224744871391589*f_l[21]+0.7071067811865475*f_l[19]; 
  tr[0] = 1.5811388300841895*f_r[22]-1.224744871391589*f_r[20]+0.7071067811865475*f_r[18]; 
  tr[1] = 1.5811388300841895*f_r[23]-1.224744871391589*f_r[21]+0.7071067811865475*f_r[19]; 
  pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tl, ul[3]); pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tr, ur[3]); 
  tl[0] = 1.5811388300841895*f_l[28]+1.224744871391589*f_l[26]+0.7071067811865475*f_l[24]; 
  tl[1] = 1.5811388300841895*f_l[29]+1.224744871391589*f_l[27]+0.7071067811865475*f_l[25]; 
  tr[0] = 1.5811388300841895*f_r[28]-1.224744871391589*f_r[26]+0.7071067811865475*f_r[24]; 
  tr[1] = 1.5811388300841895*f_r[29]-1.224744871391589*f_r[27]+0.7071067811865475*f_r[25]; 
  pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tl, ul[4]); pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tr, ur[4]); 
  tl[0] = 1.5811388300841895*f_l[34]+1.224744871391589*f_l[32]+0.7071067811865475*f_l[30]; 
  tl[1] = 1.5811388300841895*f_l[35]+1.224744871391589*f_l[33]+0.7071067811865475*f_l[31]; 
  tr[0] = 1.5811388300841895*f_r[34]-1.224744871391589*f_r[32]+0.7071067811865475*f_r[30]; 
  tr[1] = 1.5811388300841895*f_r[35]-1.224744871391589*f_r[33]+0.7071067811865475*f_r[31]; 
  pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tl, ul[5]); pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tr, ur[5]); 
  tl[0] = 1.5811388300841895*f_l[40]+1.224744871391589*f_l[38]+0.7071067811865475*f_l[36]; 
  tl[1] = 1.5811388300841895*f_l[41]+1.224744871391589*f_l[39]+0.7071067811865475*f_l[37]; 
  tr[0] = 1.5811388300841895*f_r[40]-1.224744871391589*f_r[38]+0.7071067811865475*f_r[36]; 
  tr[1] = 1.5811388300841895*f_r[41]-1.224744871391589*f_r[39]+0.7071067811865475*f_r[37]; 
  pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tl, ul[6]); pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tr, ur[6]); 
  tl[0] = 1.5811388300841895*f_l[46]+1.224744871391589*f_l[44]+0.7071067811865475*f_l[42]; 
  tl[1] = 1.5811388300841895*f_l[47]+1.224744871391589*f_l[45]+0.7071067811865475*f_l[43]; 
  tr[0] = 1.5811388300841895*f_r[46]-1.224744871391589*f_r[44]+0.7071067811865475*f_r[42]; 
  tr[1] = 1.5811388300841895*f_r[47]-1.224744871391589*f_r[45]+0.7071067811865475*f_r[43]; 
  pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tl, ul[7]); pkpm_vel_flux_vpar_1x1v_ser_p1_p1(tr, ur[7]); 
  double Fhat[8][5]; 
  for (int n = 0; n < 5; ++n) { 
    double U_l[8], U_r[8], Fn[8]; 
    for (int c = 0; c < 8; ++c) { U_l[c] = ul[c][n]; U_r[c] = ur[c][n]; } 
    pkpm_harm_node_wflux(&nv[n*PKPM_HS_NODE_NV], v_face, U_l, U_r, Fn); 
    for (int c = 0; c < 8; ++c) Fhat[c][n] = Fn[c]; 
  } 
  for (int c = 0; c < 8; ++c) pkpm_vel_flux_vpar_1x1v_ser_p1_prj(Fhat[c], &flux[c*2]); 
} 
