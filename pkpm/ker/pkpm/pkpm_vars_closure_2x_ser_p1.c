#include <gkyl_pkpm_vars_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 

// The p=1 field f (gkyl basis) at the 25 Gauss-Legendre nodes (5 per direction, row-major: the 
// last direction fastest), direction by direction. 
GKYL_CU_DH static inline void pkpm_vars_closure_2x_ser_p1_p1(const double *f, double* GKYL_RESTRICT fn) 
{ 
  double ev1[10]; 
  ev1[0] = 0.7071067811865475*f[0]-1.1098391188717989*f[1]; 
  ev1[1] = 0.7071067811865475*f[2]-1.1098391188717989*f[3]; 
  ev1[2] = 0.7071067811865475*f[0]-0.6594875259537023*f[1]; 
  ev1[3] = 0.7071067811865475*f[2]-0.6594875259537023*f[3]; 
  ev1[4] = 0.7071067811865475*f[0]; 
  ev1[5] = 0.7071067811865475*f[2]; 
  ev1[6] = 0.6594875259537023*f[1]+0.7071067811865475*f[0]; 
  ev1[7] = 0.6594875259537023*f[3]+0.7071067811865475*f[2]; 
  ev1[8] = 1.1098391188717989*f[1]+0.7071067811865475*f[0]; 
  ev1[9] = 1.1098391188717989*f[3]+0.7071067811865475*f[2]; 
  fn[0] = 0.7071067811865475*ev1[0]-1.1098391188717989*ev1[1]; 
  fn[1] = 0.7071067811865475*ev1[0]-0.6594875259537023*ev1[1]; 
  fn[2] = 0.7071067811865475*ev1[0]; 
  fn[3] = 0.6594875259537023*ev1[1]+0.7071067811865475*ev1[0]; 
  fn[4] = 1.1098391188717989*ev1[1]+0.7071067811865475*ev1[0]; 
  fn[5] = 0.7071067811865475*ev1[2]-1.1098391188717989*ev1[3]; 
  fn[6] = 0.7071067811865475*ev1[2]-0.6594875259537023*ev1[3]; 
  fn[7] = 0.7071067811865475*ev1[2]; 
  fn[8] = 0.6594875259537023*ev1[3]+0.7071067811865475*ev1[2]; 
  fn[9] = 1.1098391188717989*ev1[3]+0.7071067811865475*ev1[2]; 
  fn[10] = 0.7071067811865475*ev1[4]-1.1098391188717989*ev1[5]; 
  fn[11] = 0.7071067811865475*ev1[4]-0.6594875259537023*ev1[5]; 
  fn[12] = 0.7071067811865475*ev1[4]; 
  fn[13] = 0.6594875259537023*ev1[5]+0.7071067811865475*ev1[4]; 
  fn[14] = 1.1098391188717989*ev1[5]+0.7071067811865475*ev1[4]; 
  fn[15] = 0.7071067811865475*ev1[6]-1.1098391188717989*ev1[7]; 
  fn[16] = 0.7071067811865475*ev1[6]-0.6594875259537023*ev1[7]; 
  fn[17] = 0.7071067811865475*ev1[6]; 
  fn[18] = 0.6594875259537023*ev1[7]+0.7071067811865475*ev1[6]; 
  fn[19] = 1.1098391188717989*ev1[7]+0.7071067811865475*ev1[6]; 
  fn[20] = 0.7071067811865475*ev1[8]-1.1098391188717989*ev1[9]; 
  fn[21] = 0.7071067811865475*ev1[8]-0.6594875259537023*ev1[9]; 
  fn[22] = 0.7071067811865475*ev1[8]; 
  fn[23] = 0.6594875259537023*ev1[9]+0.7071067811865475*ev1[8]; 
  fn[24] = 1.1098391188717989*ev1[9]+0.7071067811865475*ev1[8]; 
} 

// Projection of node values onto the p=1 basis (orthonormal coefficients, gkyl order). 
GKYL_CU_DH static inline void pkpm_vars_closure_2x_ser_p1_prj(const double *v, double* GKYL_RESTRICT out) 
{ 
  double pr1[10]; 
  pr1[0] = 0.16753260706863699*v[20]+0.33844157858040363*v[15]+0.4022651910750137*v[10]+0.33844157858040363*v[5]+0.16753260706863699*v[0]; 
  pr1[1] = 0.16753260706863699*v[21]+0.33844157858040363*v[16]+0.4022651910750137*v[11]+0.33844157858040363*v[6]+0.16753260706863699*v[1]; 
  pr1[2] = 0.16753260706863699*v[22]+0.33844157858040363*v[17]+0.4022651910750137*v[12]+0.33844157858040363*v[7]+0.16753260706863699*v[2]; 
  pr1[3] = 0.16753260706863699*v[23]+0.33844157858040363*v[18]+0.4022651910750137*v[13]+0.33844157858040363*v[8]+0.16753260706863699*v[3]; 
  pr1[4] = 0.16753260706863699*v[24]+0.33844157858040363*v[19]+0.4022651910750137*v[14]+0.33844157858040363*v[9]+0.16753260706863699*v[4]; 
  pr1[5] = 0.26295072534780084*v[20]+0.31564963775813704*v[15]-0.31564963775813704*v[5]-0.26295072534780084*v[0]; 
  pr1[6] = 0.26295072534780084*v[21]+0.31564963775813704*v[16]-0.31564963775813704*v[6]-0.26295072534780084*v[1]; 
  pr1[7] = 0.26295072534780084*v[22]+0.31564963775813704*v[17]-0.31564963775813704*v[7]-0.26295072534780084*v[2]; 
  pr1[8] = 0.26295072534780084*v[23]+0.31564963775813704*v[18]-0.31564963775813704*v[8]-0.26295072534780084*v[3]; 
  pr1[9] = 0.26295072534780084*v[24]+0.31564963775813704*v[19]-0.31564963775813704*v[9]-0.26295072534780084*v[4]; 
  out[0] = 0.16753260706863699*pr1[4]+0.33844157858040363*pr1[3]+0.4022651910750137*pr1[2]+0.33844157858040363*pr1[1]+0.16753260706863699*pr1[0]; 
  out[2] = 0.26295072534780084*pr1[4]+0.31564963775813704*pr1[3]-0.31564963775813704*pr1[1]-0.26295072534780084*pr1[0]; 
  out[1] = 0.16753260706863699*pr1[9]+0.33844157858040363*pr1[8]+0.4022651910750137*pr1[7]+0.33844157858040363*pr1[6]+0.16753260706863699*pr1[5]; 
  out[3] = 0.26295072534780084*pr1[9]+0.31564963775813704*pr1[8]-0.31564963775813704*pr1[6]-0.26295072534780084*pr1[5]; 
} 

GKYL_CU_DH int pkpm_vars_closure_2x_ser_p1(double w_lo, double w_up, double dv, const double *norm_moms, const double *rhs_moms,
  const double *bvar, const double *wn, const double *theta, const double *alpha_0_prev, const double *c_h_prev,
  const double *f_lo, const double *f_up,
  double* GKYL_RESTRICT c_0, double* GKYL_RESTRICT alpha_0, double* GKYL_RESTRICT c_h, double* GKYL_RESTRICT a_perp,
  double* GKYL_RESTRICT flux_lo, double* GKYL_RESTRICT flux_up, double* GKYL_RESTRICT cfl,
  double* GKYL_RESTRICT cond) 
{ 
  // w_lo, w_up:   Velocity domain edges; dv: velocity cell size. 
  // norm_moms:    Input [M0, M2, M0_G, M1, M3, M1_G, mu (x, y, z)] moments of the normalized distributions (stage state). 
  // rhs_moms:     Input moments of the stage-A right-hand sides, same layout. 
  // bvar:         Magnetic field unit vector [bx, by, bz]. 
  // wn:           Raw w-flux values at the 25 w-face nodes [node][PKPM_HS_W_END] (the build of the stored operator). 
  // theta:        T_perp/T. 
  // alpha_0_prev, c_h_prev: closures of the previous stage (may alias the outputs), for the edge speeds. 
  // f_lo, f_up:   [F, G, M0_perp, N_perp] of the lowest and uppermost velocity cells of this configuration cell. 
  // c_0, alpha_0, c_h, a_perp: Output closures at the 4 Gauss-Legendre nodes, interpolated to p=1. 
  // flux_lo, flux_up: Output edge fluxes [8 blocks] in the stored face format. 
  // cfl:          Output [0]: the time-step rate of the w fluxes and closure sources of the whole column 
  //               (absolute coefficients at the larger |w| edge, plan 17.27); [1] = 0. 
  // cond:         Output [1]: the smallest determinant of the 3x3 system over the nodes relative to 
  //               its edge-free value (1 when the edge terms vanish; small or negative when the 
  //               velocity-domain edges carry a non-negligible trace: too small velocity extents). 
  // Returns the number of nodes where the closures are undefined (everything is zero in the cell then). 
 
  const double dv1 = 2.0/dv; 
  const double *n0 = &norm_moms[0]; 
  const double *n2 = &norm_moms[4]; 
  const double *nG = &norm_moms[8]; 
  const double *n1 = &norm_moms[12]; 
  const double *r0 = &rhs_moms[0]; 
  const double *r2 = &rhs_moms[4]; 
  const double *rG = &rhs_moms[8]; 
  const double *r1 = &rhs_moms[12]; 
  const double *mu0 = &norm_moms[24]; 
  const double *mu1 = &norm_moms[28]; 
  const double *mu2 = &norm_moms[32]; 
  const double *rM0 = &rhs_moms[24]; 
  const double *rM1 = &rhs_moms[28]; 
  const double *rM2 = &rhs_moms[32]; 
  const double *b0 = &bvar[0]; 
  const double *b1 = &bvar[4]; 
  const double *b2 = &bvar[8]; 
 
  // The edge fluxes at the 25 w-face nodes, constructed as on the interior w faces (traces, raw 
  // node values, HLL against the empty exterior with the previous closures): 
  // Phi_e = K_e + (alpha_0 + c w_e) g_e, K_e = wf (k_e U + T(U)) + wu U, g_e = wf U, e = lower, upper. 
  double tr[4], ulo[8][25], uup[8][25]; 
  tr[0] = 1.5811388300841895*f_lo[7]-1.224744871391589*f_lo[3]+0.7071067811865475*f_lo[0]; 
  tr[1] = 1.5811388300841895*f_lo[9]-1.224744871391589*f_lo[5]+0.7071067811865475*f_lo[1]; 
  tr[2] = 1.5811388300841895*f_lo[10]-1.224744871391589*f_lo[6]+0.7071067811865475*f_lo[2]; 
  tr[3] = 1.5811388300841895*f_lo[11]-1.224744871391589*f_lo[8]+0.7071067811865475*f_lo[4]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, ulo[0]); 
  tr[0] = 1.5811388300841895*f_up[7]+1.224744871391589*f_up[3]+0.7071067811865475*f_up[0]; 
  tr[1] = 1.5811388300841895*f_up[9]+1.224744871391589*f_up[5]+0.7071067811865475*f_up[1]; 
  tr[2] = 1.5811388300841895*f_up[10]+1.224744871391589*f_up[6]+0.7071067811865475*f_up[2]; 
  tr[3] = 1.5811388300841895*f_up[11]+1.224744871391589*f_up[8]+0.7071067811865475*f_up[4]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, uup[0]); 
  tr[0] = 1.5811388300841895*f_lo[19]-1.224744871391589*f_lo[15]+0.7071067811865475*f_lo[12]; 
  tr[1] = 1.5811388300841895*f_lo[21]-1.224744871391589*f_lo[17]+0.7071067811865475*f_lo[13]; 
  tr[2] = 1.5811388300841895*f_lo[22]-1.224744871391589*f_lo[18]+0.7071067811865475*f_lo[14]; 
  tr[3] = 1.5811388300841895*f_lo[23]-1.224744871391589*f_lo[20]+0.7071067811865475*f_lo[16]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, ulo[1]); 
  tr[0] = 1.5811388300841895*f_up[19]+1.224744871391589*f_up[15]+0.7071067811865475*f_up[12]; 
  tr[1] = 1.5811388300841895*f_up[21]+1.224744871391589*f_up[17]+0.7071067811865475*f_up[13]; 
  tr[2] = 1.5811388300841895*f_up[22]+1.224744871391589*f_up[18]+0.7071067811865475*f_up[14]; 
  tr[3] = 1.5811388300841895*f_up[23]+1.224744871391589*f_up[20]+0.7071067811865475*f_up[16]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, uup[1]); 
  tr[0] = 1.5811388300841895*f_lo[31]-1.224744871391589*f_lo[27]+0.7071067811865475*f_lo[24]; 
  tr[1] = 1.5811388300841895*f_lo[33]-1.224744871391589*f_lo[29]+0.7071067811865475*f_lo[25]; 
  tr[2] = 1.5811388300841895*f_lo[34]-1.224744871391589*f_lo[30]+0.7071067811865475*f_lo[26]; 
  tr[3] = 1.5811388300841895*f_lo[35]-1.224744871391589*f_lo[32]+0.7071067811865475*f_lo[28]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, ulo[2]); 
  tr[0] = 1.5811388300841895*f_up[31]+1.224744871391589*f_up[27]+0.7071067811865475*f_up[24]; 
  tr[1] = 1.5811388300841895*f_up[33]+1.224744871391589*f_up[29]+0.7071067811865475*f_up[25]; 
  tr[2] = 1.5811388300841895*f_up[34]+1.224744871391589*f_up[30]+0.7071067811865475*f_up[26]; 
  tr[3] = 1.5811388300841895*f_up[35]+1.224744871391589*f_up[32]+0.7071067811865475*f_up[28]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, uup[2]); 
  tr[0] = 1.5811388300841895*f_lo[43]-1.224744871391589*f_lo[39]+0.7071067811865475*f_lo[36]; 
  tr[1] = 1.5811388300841895*f_lo[45]-1.224744871391589*f_lo[41]+0.7071067811865475*f_lo[37]; 
  tr[2] = 1.5811388300841895*f_lo[46]-1.224744871391589*f_lo[42]+0.7071067811865475*f_lo[38]; 
  tr[3] = 1.5811388300841895*f_lo[47]-1.224744871391589*f_lo[44]+0.7071067811865475*f_lo[40]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, ulo[3]); 
  tr[0] = 1.5811388300841895*f_up[43]+1.224744871391589*f_up[39]+0.7071067811865475*f_up[36]; 
  tr[1] = 1.5811388300841895*f_up[45]+1.224744871391589*f_up[41]+0.7071067811865475*f_up[37]; 
  tr[2] = 1.5811388300841895*f_up[46]+1.224744871391589*f_up[42]+0.7071067811865475*f_up[38]; 
  tr[3] = 1.5811388300841895*f_up[47]+1.224744871391589*f_up[44]+0.7071067811865475*f_up[40]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, uup[3]); 
  tr[0] = 1.5811388300841895*f_lo[55]-1.224744871391589*f_lo[51]+0.7071067811865475*f_lo[48]; 
  tr[1] = 1.5811388300841895*f_lo[57]-1.224744871391589*f_lo[53]+0.7071067811865475*f_lo[49]; 
  tr[2] = 1.5811388300841895*f_lo[58]-1.224744871391589*f_lo[54]+0.7071067811865475*f_lo[50]; 
  tr[3] = 1.5811388300841895*f_lo[59]-1.224744871391589*f_lo[56]+0.7071067811865475*f_lo[52]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, ulo[4]); 
  tr[0] = 1.5811388300841895*f_up[55]+1.224744871391589*f_up[51]+0.7071067811865475*f_up[48]; 
  tr[1] = 1.5811388300841895*f_up[57]+1.224744871391589*f_up[53]+0.7071067811865475*f_up[49]; 
  tr[2] = 1.5811388300841895*f_up[58]+1.224744871391589*f_up[54]+0.7071067811865475*f_up[50]; 
  tr[3] = 1.5811388300841895*f_up[59]+1.224744871391589*f_up[56]+0.7071067811865475*f_up[52]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, uup[4]); 
  tr[0] = 1.5811388300841895*f_lo[67]-1.224744871391589*f_lo[63]+0.7071067811865475*f_lo[60]; 
  tr[1] = 1.5811388300841895*f_lo[69]-1.224744871391589*f_lo[65]+0.7071067811865475*f_lo[61]; 
  tr[2] = 1.5811388300841895*f_lo[70]-1.224744871391589*f_lo[66]+0.7071067811865475*f_lo[62]; 
  tr[3] = 1.5811388300841895*f_lo[71]-1.224744871391589*f_lo[68]+0.7071067811865475*f_lo[64]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, ulo[5]); 
  tr[0] = 1.5811388300841895*f_up[67]+1.224744871391589*f_up[63]+0.7071067811865475*f_up[60]; 
  tr[1] = 1.5811388300841895*f_up[69]+1.224744871391589*f_up[65]+0.7071067811865475*f_up[61]; 
  tr[2] = 1.5811388300841895*f_up[70]+1.224744871391589*f_up[66]+0.7071067811865475*f_up[62]; 
  tr[3] = 1.5811388300841895*f_up[71]+1.224744871391589*f_up[68]+0.7071067811865475*f_up[64]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, uup[5]); 
  tr[0] = 1.5811388300841895*f_lo[79]-1.224744871391589*f_lo[75]+0.7071067811865475*f_lo[72]; 
  tr[1] = 1.5811388300841895*f_lo[81]-1.224744871391589*f_lo[77]+0.7071067811865475*f_lo[73]; 
  tr[2] = 1.5811388300841895*f_lo[82]-1.224744871391589*f_lo[78]+0.7071067811865475*f_lo[74]; 
  tr[3] = 1.5811388300841895*f_lo[83]-1.224744871391589*f_lo[80]+0.7071067811865475*f_lo[76]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, ulo[6]); 
  tr[0] = 1.5811388300841895*f_up[79]+1.224744871391589*f_up[75]+0.7071067811865475*f_up[72]; 
  tr[1] = 1.5811388300841895*f_up[81]+1.224744871391589*f_up[77]+0.7071067811865475*f_up[73]; 
  tr[2] = 1.5811388300841895*f_up[82]+1.224744871391589*f_up[78]+0.7071067811865475*f_up[74]; 
  tr[3] = 1.5811388300841895*f_up[83]+1.224744871391589*f_up[80]+0.7071067811865475*f_up[76]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, uup[6]); 
  tr[0] = 1.5811388300841895*f_lo[91]-1.224744871391589*f_lo[87]+0.7071067811865475*f_lo[84]; 
  tr[1] = 1.5811388300841895*f_lo[93]-1.224744871391589*f_lo[89]+0.7071067811865475*f_lo[85]; 
  tr[2] = 1.5811388300841895*f_lo[94]-1.224744871391589*f_lo[90]+0.7071067811865475*f_lo[86]; 
  tr[3] = 1.5811388300841895*f_lo[95]-1.224744871391589*f_lo[92]+0.7071067811865475*f_lo[88]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, ulo[7]); 
  tr[0] = 1.5811388300841895*f_up[91]+1.224744871391589*f_up[87]+0.7071067811865475*f_up[84]; 
  tr[1] = 1.5811388300841895*f_up[93]+1.224744871391589*f_up[89]+0.7071067811865475*f_up[85]; 
  tr[2] = 1.5811388300841895*f_up[94]+1.224744871391589*f_up[90]+0.7071067811865475*f_up[86]; 
  tr[3] = 1.5811388300841895*f_up[95]+1.224744871391589*f_up[92]+0.7071067811865475*f_up[88]; 
  pkpm_vars_closure_2x_ser_p1_p1(tr, uup[7]); 
  double th_n[25]; 
  pkpm_vars_closure_2x_ser_p1_p1(theta, th_n); 
  double b0_n[25]; 
  pkpm_vars_closure_2x_ser_p1_p1(&bvar[0], b0_n); 
  double b1_n[25]; 
  pkpm_vars_closure_2x_ser_p1_p1(&bvar[4], b1_n); 
  double b2_n[25]; 
  pkpm_vars_closure_2x_ser_p1_p1(&bvar[8], b2_n); 
  double a0_n[25]; 
  pkpm_vars_closure_2x_ser_p1_p1(alpha_0_prev, a0_n); 
  double ch_n[25]; 
  pkpm_vars_closure_2x_ser_p1_p1(c_h_prev, ch_n); 
  double Kn[2][8][25], gn[2][8][25], kw[25], rest[25]; 
  const double wm = fmax(fabs(w_lo), fabs(w_up)); 
  for (int n = 0; n < 25; ++n) { 
    double nv[PKPM_HS_NODE_NV]; 
    for (int s = 0; s < PKPM_HS_W_END; ++s) nv[s] = wn[n*PKPM_HS_W_END + s]; 
    nv[PKPM_HS_NODE_THETA] = th_n[n]; 
    nv[PKPM_HS_NODE_B] = b0_n[n]; nv[PKPM_HS_NODE_B + 1] = b1_n[n]; nv[PKPM_HS_NODE_B + 2] = b2_n[n]; 
    nv[PKPM_HS_NODE_A0] = a0_n[n]; nv[PKPM_HS_NODE_CH] = ch_n[n]; 
    pkpm_harm_node_spec(nv); 
    kw[n] = nv[PKPM_HS_KW]; 
    rest[n] = pkpm_harm_node_speed_rest(nv, wm); 
    for (int e = 0; e < 2; ++e) { 
      const double w = e == 0 ? w_lo : w_up; 
      double U[8], T[8], lo, hi; 
      for (int c = 0; c < 8; ++c) U[c] = e == 0 ? ulo[c][n] : uup[c][n]; 
      pkpm_harm_node_coupling(nv, w, U, T); 
      pkpm_harm_node_bounds(nv, w, &lo, &hi); 
      const double k = pkpm_harm_node_alpha_k(nv, w); 
      const double ap = nv[PKPM_HS_NODE_A0] + nv[PKPM_HS_NODE_CH]*w + k; 
      const double sl = ap + lo, sr = ap + hi; 
      // the HLL flux against the empty exterior (below the lower edge, above the upper one): 
      // wf (alpha U + T(U)) + wu U 
      double wf = 0.0, wu = 0.0; 
      if (e == 0) { 
        if (sr < 0.0) { 
          wf = 1.0; 
        } else if (sl < 0.0) { 
          wf = -sl/(sr - sl); wu = sl*sr/(sr - sl); 
        } 
      } else { 
        if (sl > 0.0) { 
          wf = 1.0; 
        } else if (sr > 0.0) { 
          wf = sr/(sr - sl); wu = -sl*sr/(sr - sl); 
        } 
      } 
      for (int c = 0; c < 8; ++c) { 
        Kn[e][c][n] = wf*(k*U[c] + T[c]) + wu*U[c]; 
        gn[e][c][n] = wf*U[c]; 
      } 
    } 
  } 
  double Kh[2][32], gh[2][32]; 
  for (int e = 0; e < 2; ++e) { 
    for (int c = 0; c < 8; ++c) { 
      pkpm_vars_closure_2x_ser_p1_prj(Kn[e][c], &Kh[e][c*4]); 
      pkpm_vars_closure_2x_ser_p1_prj(gn[e][c], &gh[e][c*4]); 
    } 
  } 
 
  // The closures at the 4 Gauss-Legendre nodes, with the edge fluxes' p=1 projections there. 
  double c0_s[4] = {0.0}, a0_s[4] = {0.0}, c_s[4] = {0.0}; 
  double ap_s[3][4] = {{0.0}}, flo_s[8][4] = {{0.0}}, fup_s[8][4] = {{0.0}}; 
  int num_null = 0; 
  double cond_min = 1.0; 
 
  { 
    const double n0_q = 0.5*n0[3]-0.5*n0[2]-0.5*n0[1]+0.5*n0[0]; 
    const double n1_q = 0.5*n1[3]-0.5*n1[2]-0.5*n1[1]+0.5*n1[0]; 
    const double n2_q = 0.5*n2[3]-0.5*n2[2]-0.5*n2[1]+0.5*n2[0]; 
    const double nG_q = 0.5*nG[3]-0.5*nG[2]-0.5*nG[1]+0.5*nG[0]; 
    const double r0_q = 0.5*r0[3]-0.5*r0[2]-0.5*r0[1]+0.5*r0[0]; 
    const double r1_q = 0.5*r1[3]-0.5*r1[2]-0.5*r1[1]+0.5*r1[0]; 
    const double r2_q = 0.5*r2[3]-0.5*r2[2]-0.5*r2[1]+0.5*r2[0] + 2.0*(0.5*rG[3]-0.5*rG[2]-0.5*rG[1]+0.5*rG[0]); 
    const double mu_q[3] = {0.5*mu0[3]-0.5*mu0[2]-0.5*mu0[1]+0.5*mu0[0], 0.5*mu1[3]-0.5*mu1[2]-0.5*mu1[1]+0.5*mu1[0], 0.5*mu2[3]-0.5*mu2[2]-0.5*mu2[1]+0.5*mu2[0]}; 
    const double rM_q[3] = {0.5*rM0[3]-0.5*rM0[2]-0.5*rM0[1]+0.5*rM0[0], 0.5*rM1[3]-0.5*rM1[2]-0.5*rM1[1]+0.5*rM1[0], 0.5*rM2[3]-0.5*rM2[2]-0.5*rM2[1]+0.5*rM2[0]}; 
    const double b_q[3] = {0.5*b0[3]-0.5*b0[2]-0.5*b0[1]+0.5*b0[0], 0.5*b1[3]-0.5*b1[2]-0.5*b1[1]+0.5*b1[0], 0.5*b2[3]-0.5*b2[2]-0.5*b2[1]+0.5*b2[0]}; 
    // the p=1 projections of the edge fluxes' parts at the node 
    double Klo[8], Kup[8], glo[8], gup[8]; 
    for (int c = 0; c < 8; ++c) { 
      const double *h = &Kh[0][c*4]; 
      Klo[c] = 0.5*h[3]-0.5*h[2]-0.5*h[1]+0.5*h[0]; 
      h = &Kh[1][c*4]; 
      Kup[c] = 0.5*h[3]-0.5*h[2]-0.5*h[1]+0.5*h[0]; 
      h = &gh[0][c*4]; 
      glo[c] = 0.5*h[3]-0.5*h[2]-0.5*h[1]+0.5*h[0]; 
      h = &gh[1][c*4]; 
      gup[c] = 0.5*h[3]-0.5*h[2]-0.5*h[1]+0.5*h[0]; 
    } 
    // edge contributions to the moment rates: the k-th moment gains sgn_e w_e^k Phi_e per edge, 
    // sgn_e = +1 lo, -1 up: closure part (alpha_0 + c w_e) g_e and known part K_e 
    const double eF0 = glo[0] - gup[0]; 
    const double eF1 = w_lo*glo[0] - w_up*gup[0]; 
    const double eF2 = w_lo*w_lo*glo[0] - w_up*w_up*gup[0]; 
    const double eF3 = w_lo*w_lo*w_lo*glo[0] - w_up*w_up*w_up*gup[0]; 
    const double eG0 = glo[1] - gup[1]; 
    const double eG1 = w_lo*glo[1] - w_up*gup[1]; 
    const double KF0 = Klo[0] - Kup[0]; 
    const double KF1 = w_lo*Klo[0] - w_up*Kup[0]; 
    const double KF2 = w_lo*w_lo*Klo[0] - w_up*w_up*Kup[0]; 
    const double KG0 = Klo[1] - Kup[1]; 
    // perpendicular rows M0 a_perp + (c_0 + c) Q + alpha_0 A_a + c A_c = -P (b is a unit vector here) 
    double P[3], Q[3], Aa[3], Ac[3], v[3]; 
    for (int a = 0; a < 3; ++a) v[a] = rM_q[a] + Klo[2 + a] - Kup[2 + a]; 
    pkpm_harm_perp(b_q, v, P); 
    pkpm_harm_perp(b_q, mu_q, Q); 
    for (int a = 0; a < 3; ++a) v[a] = glo[2 + a] - gup[2 + a]; 
    pkpm_harm_perp(b_q, v, Aa); 
    for (int a = 0; a < 3; ++a) v[a] = w_lo*glo[2 + a] - w_up*gup[2 + a]; 
    pkpm_harm_perp(b_q, v, Ac); 
    // eliminated into the energy row, which carries 2 mu . a_perp (Schur complement) 
    const double fac = n0_q > 0.0 ? 2.0/n0_q : 0.0; 
    const double sQQ = pkpm_harm_dot(Q, Q), sQP = pkpm_harm_dot(Q, P), sQA = pkpm_harm_dot(Q, Aa), sQC = pkpm_harm_dot(Q, Ac); 
    // 3x3 system A [c_0, alpha_0, c] = rhs: stage-A rates + interior stage-B rates + edge rates = 0 
    const double A00 = n0_q, A01 = eF0, A02 = eF1; 
    const double A10 = n1_q, A11 = n0_q + eF1, A12 = n1_q + eF2; 
    const double A20 = n2_q + 2.0*nG_q - fac*sQQ, A21 = 2.0*n1_q + eF2 + 2.0*eG0 - fac*sQA, A22 = 2.0*(n2_q + 2.0*nG_q) + eF3 + 2.0*eG1 - fac*(sQQ + sQC); 
    const double b0 = -r0_q - KF0, b1 = -r1_q - KF1, b2 = -r2_q - KF2 - 2.0*KG0 + fac*sQP; 
    const double det = A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20); 
    // conditioning proxy: the determinant relative to its edge-free value 
    // n0 (2 n0 (M2 + 2 M0G) - 2 M1^2) (positive for admissible moments); a small or negative 
    // ratio means the traces at the velocity-domain edges are not small 
    const double det0 = n0_q*(2.0*n0_q*(n2_q + 2.0*nG_q) - 2.0*n1_q*n1_q); 
    if (det0 > 0.0) cond_min = fmin(cond_min, det/det0); 
    if (n0_q > 0.0 && det > 0.0) { 
      c0_s[0] = (b0*(A11*A22 - A12*A21) - A01*(b1*A22 - A12*b2) + A02*(b1*A21 - A11*b2))/det; 
      a0_s[0] = (A00*(b1*A22 - A12*b2) - b0*(A10*A22 - A12*A20) + A02*(A10*b2 - b1*A20))/det; 
      c_s[0] = (A00*(A11*b2 - b1*A21) - A01*(A10*b2 - b1*A20) + b0*(A10*A21 - A11*A20))/det; 
      for (int a = 0; a < 3; ++a) 
        ap_s[a][0] = -(P[a] + (c0_s[0] + c_s[0])*Q[a] + Aa[a]*a0_s[0] + Ac[a]*c_s[0])/n0_q; 
      for (int c = 0; c < 8; ++c) { 
        flo_s[c][0] = (a0_s[0] + c_s[0]*w_lo)*glo[c] + Klo[c]; 
        fup_s[c][0] = (a0_s[0] + c_s[0]*w_up)*gup[c] + Kup[c]; 
      } 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double n0_q = -(0.5*n0[3])+0.5*n0[2]-0.5*n0[1]+0.5*n0[0]; 
    const double n1_q = -(0.5*n1[3])+0.5*n1[2]-0.5*n1[1]+0.5*n1[0]; 
    const double n2_q = -(0.5*n2[3])+0.5*n2[2]-0.5*n2[1]+0.5*n2[0]; 
    const double nG_q = -(0.5*nG[3])+0.5*nG[2]-0.5*nG[1]+0.5*nG[0]; 
    const double r0_q = -(0.5*r0[3])+0.5*r0[2]-0.5*r0[1]+0.5*r0[0]; 
    const double r1_q = -(0.5*r1[3])+0.5*r1[2]-0.5*r1[1]+0.5*r1[0]; 
    const double r2_q = -(0.5*r2[3])+0.5*r2[2]-0.5*r2[1]+0.5*r2[0] + 2.0*(-(0.5*rG[3])+0.5*rG[2]-0.5*rG[1]+0.5*rG[0]); 
    const double mu_q[3] = {-(0.5*mu0[3])+0.5*mu0[2]-0.5*mu0[1]+0.5*mu0[0], -(0.5*mu1[3])+0.5*mu1[2]-0.5*mu1[1]+0.5*mu1[0], -(0.5*mu2[3])+0.5*mu2[2]-0.5*mu2[1]+0.5*mu2[0]}; 
    const double rM_q[3] = {-(0.5*rM0[3])+0.5*rM0[2]-0.5*rM0[1]+0.5*rM0[0], -(0.5*rM1[3])+0.5*rM1[2]-0.5*rM1[1]+0.5*rM1[0], -(0.5*rM2[3])+0.5*rM2[2]-0.5*rM2[1]+0.5*rM2[0]}; 
    const double b_q[3] = {-(0.5*b0[3])+0.5*b0[2]-0.5*b0[1]+0.5*b0[0], -(0.5*b1[3])+0.5*b1[2]-0.5*b1[1]+0.5*b1[0], -(0.5*b2[3])+0.5*b2[2]-0.5*b2[1]+0.5*b2[0]}; 
    // the p=1 projections of the edge fluxes' parts at the node 
    double Klo[8], Kup[8], glo[8], gup[8]; 
    for (int c = 0; c < 8; ++c) { 
      const double *h = &Kh[0][c*4]; 
      Klo[c] = -(0.5*h[3])+0.5*h[2]-0.5*h[1]+0.5*h[0]; 
      h = &Kh[1][c*4]; 
      Kup[c] = -(0.5*h[3])+0.5*h[2]-0.5*h[1]+0.5*h[0]; 
      h = &gh[0][c*4]; 
      glo[c] = -(0.5*h[3])+0.5*h[2]-0.5*h[1]+0.5*h[0]; 
      h = &gh[1][c*4]; 
      gup[c] = -(0.5*h[3])+0.5*h[2]-0.5*h[1]+0.5*h[0]; 
    } 
    // edge contributions to the moment rates: the k-th moment gains sgn_e w_e^k Phi_e per edge, 
    // sgn_e = +1 lo, -1 up: closure part (alpha_0 + c w_e) g_e and known part K_e 
    const double eF0 = glo[0] - gup[0]; 
    const double eF1 = w_lo*glo[0] - w_up*gup[0]; 
    const double eF2 = w_lo*w_lo*glo[0] - w_up*w_up*gup[0]; 
    const double eF3 = w_lo*w_lo*w_lo*glo[0] - w_up*w_up*w_up*gup[0]; 
    const double eG0 = glo[1] - gup[1]; 
    const double eG1 = w_lo*glo[1] - w_up*gup[1]; 
    const double KF0 = Klo[0] - Kup[0]; 
    const double KF1 = w_lo*Klo[0] - w_up*Kup[0]; 
    const double KF2 = w_lo*w_lo*Klo[0] - w_up*w_up*Kup[0]; 
    const double KG0 = Klo[1] - Kup[1]; 
    // perpendicular rows M0 a_perp + (c_0 + c) Q + alpha_0 A_a + c A_c = -P (b is a unit vector here) 
    double P[3], Q[3], Aa[3], Ac[3], v[3]; 
    for (int a = 0; a < 3; ++a) v[a] = rM_q[a] + Klo[2 + a] - Kup[2 + a]; 
    pkpm_harm_perp(b_q, v, P); 
    pkpm_harm_perp(b_q, mu_q, Q); 
    for (int a = 0; a < 3; ++a) v[a] = glo[2 + a] - gup[2 + a]; 
    pkpm_harm_perp(b_q, v, Aa); 
    for (int a = 0; a < 3; ++a) v[a] = w_lo*glo[2 + a] - w_up*gup[2 + a]; 
    pkpm_harm_perp(b_q, v, Ac); 
    // eliminated into the energy row, which carries 2 mu . a_perp (Schur complement) 
    const double fac = n0_q > 0.0 ? 2.0/n0_q : 0.0; 
    const double sQQ = pkpm_harm_dot(Q, Q), sQP = pkpm_harm_dot(Q, P), sQA = pkpm_harm_dot(Q, Aa), sQC = pkpm_harm_dot(Q, Ac); 
    // 3x3 system A [c_0, alpha_0, c] = rhs: stage-A rates + interior stage-B rates + edge rates = 0 
    const double A00 = n0_q, A01 = eF0, A02 = eF1; 
    const double A10 = n1_q, A11 = n0_q + eF1, A12 = n1_q + eF2; 
    const double A20 = n2_q + 2.0*nG_q - fac*sQQ, A21 = 2.0*n1_q + eF2 + 2.0*eG0 - fac*sQA, A22 = 2.0*(n2_q + 2.0*nG_q) + eF3 + 2.0*eG1 - fac*(sQQ + sQC); 
    const double b0 = -r0_q - KF0, b1 = -r1_q - KF1, b2 = -r2_q - KF2 - 2.0*KG0 + fac*sQP; 
    const double det = A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20); 
    // conditioning proxy: the determinant relative to its edge-free value 
    // n0 (2 n0 (M2 + 2 M0G) - 2 M1^2) (positive for admissible moments); a small or negative 
    // ratio means the traces at the velocity-domain edges are not small 
    const double det0 = n0_q*(2.0*n0_q*(n2_q + 2.0*nG_q) - 2.0*n1_q*n1_q); 
    if (det0 > 0.0) cond_min = fmin(cond_min, det/det0); 
    if (n0_q > 0.0 && det > 0.0) { 
      c0_s[1] = (b0*(A11*A22 - A12*A21) - A01*(b1*A22 - A12*b2) + A02*(b1*A21 - A11*b2))/det; 
      a0_s[1] = (A00*(b1*A22 - A12*b2) - b0*(A10*A22 - A12*A20) + A02*(A10*b2 - b1*A20))/det; 
      c_s[1] = (A00*(A11*b2 - b1*A21) - A01*(A10*b2 - b1*A20) + b0*(A10*A21 - A11*A20))/det; 
      for (int a = 0; a < 3; ++a) 
        ap_s[a][1] = -(P[a] + (c0_s[1] + c_s[1])*Q[a] + Aa[a]*a0_s[1] + Ac[a]*c_s[1])/n0_q; 
      for (int c = 0; c < 8; ++c) { 
        flo_s[c][1] = (a0_s[1] + c_s[1]*w_lo)*glo[c] + Klo[c]; 
        fup_s[c][1] = (a0_s[1] + c_s[1]*w_up)*gup[c] + Kup[c]; 
      } 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double n0_q = -(0.5*n0[3])-0.5*n0[2]+0.5*n0[1]+0.5*n0[0]; 
    const double n1_q = -(0.5*n1[3])-0.5*n1[2]+0.5*n1[1]+0.5*n1[0]; 
    const double n2_q = -(0.5*n2[3])-0.5*n2[2]+0.5*n2[1]+0.5*n2[0]; 
    const double nG_q = -(0.5*nG[3])-0.5*nG[2]+0.5*nG[1]+0.5*nG[0]; 
    const double r0_q = -(0.5*r0[3])-0.5*r0[2]+0.5*r0[1]+0.5*r0[0]; 
    const double r1_q = -(0.5*r1[3])-0.5*r1[2]+0.5*r1[1]+0.5*r1[0]; 
    const double r2_q = -(0.5*r2[3])-0.5*r2[2]+0.5*r2[1]+0.5*r2[0] + 2.0*(-(0.5*rG[3])-0.5*rG[2]+0.5*rG[1]+0.5*rG[0]); 
    const double mu_q[3] = {-(0.5*mu0[3])-0.5*mu0[2]+0.5*mu0[1]+0.5*mu0[0], -(0.5*mu1[3])-0.5*mu1[2]+0.5*mu1[1]+0.5*mu1[0], -(0.5*mu2[3])-0.5*mu2[2]+0.5*mu2[1]+0.5*mu2[0]}; 
    const double rM_q[3] = {-(0.5*rM0[3])-0.5*rM0[2]+0.5*rM0[1]+0.5*rM0[0], -(0.5*rM1[3])-0.5*rM1[2]+0.5*rM1[1]+0.5*rM1[0], -(0.5*rM2[3])-0.5*rM2[2]+0.5*rM2[1]+0.5*rM2[0]}; 
    const double b_q[3] = {-(0.5*b0[3])-0.5*b0[2]+0.5*b0[1]+0.5*b0[0], -(0.5*b1[3])-0.5*b1[2]+0.5*b1[1]+0.5*b1[0], -(0.5*b2[3])-0.5*b2[2]+0.5*b2[1]+0.5*b2[0]}; 
    // the p=1 projections of the edge fluxes' parts at the node 
    double Klo[8], Kup[8], glo[8], gup[8]; 
    for (int c = 0; c < 8; ++c) { 
      const double *h = &Kh[0][c*4]; 
      Klo[c] = -(0.5*h[3])-0.5*h[2]+0.5*h[1]+0.5*h[0]; 
      h = &Kh[1][c*4]; 
      Kup[c] = -(0.5*h[3])-0.5*h[2]+0.5*h[1]+0.5*h[0]; 
      h = &gh[0][c*4]; 
      glo[c] = -(0.5*h[3])-0.5*h[2]+0.5*h[1]+0.5*h[0]; 
      h = &gh[1][c*4]; 
      gup[c] = -(0.5*h[3])-0.5*h[2]+0.5*h[1]+0.5*h[0]; 
    } 
    // edge contributions to the moment rates: the k-th moment gains sgn_e w_e^k Phi_e per edge, 
    // sgn_e = +1 lo, -1 up: closure part (alpha_0 + c w_e) g_e and known part K_e 
    const double eF0 = glo[0] - gup[0]; 
    const double eF1 = w_lo*glo[0] - w_up*gup[0]; 
    const double eF2 = w_lo*w_lo*glo[0] - w_up*w_up*gup[0]; 
    const double eF3 = w_lo*w_lo*w_lo*glo[0] - w_up*w_up*w_up*gup[0]; 
    const double eG0 = glo[1] - gup[1]; 
    const double eG1 = w_lo*glo[1] - w_up*gup[1]; 
    const double KF0 = Klo[0] - Kup[0]; 
    const double KF1 = w_lo*Klo[0] - w_up*Kup[0]; 
    const double KF2 = w_lo*w_lo*Klo[0] - w_up*w_up*Kup[0]; 
    const double KG0 = Klo[1] - Kup[1]; 
    // perpendicular rows M0 a_perp + (c_0 + c) Q + alpha_0 A_a + c A_c = -P (b is a unit vector here) 
    double P[3], Q[3], Aa[3], Ac[3], v[3]; 
    for (int a = 0; a < 3; ++a) v[a] = rM_q[a] + Klo[2 + a] - Kup[2 + a]; 
    pkpm_harm_perp(b_q, v, P); 
    pkpm_harm_perp(b_q, mu_q, Q); 
    for (int a = 0; a < 3; ++a) v[a] = glo[2 + a] - gup[2 + a]; 
    pkpm_harm_perp(b_q, v, Aa); 
    for (int a = 0; a < 3; ++a) v[a] = w_lo*glo[2 + a] - w_up*gup[2 + a]; 
    pkpm_harm_perp(b_q, v, Ac); 
    // eliminated into the energy row, which carries 2 mu . a_perp (Schur complement) 
    const double fac = n0_q > 0.0 ? 2.0/n0_q : 0.0; 
    const double sQQ = pkpm_harm_dot(Q, Q), sQP = pkpm_harm_dot(Q, P), sQA = pkpm_harm_dot(Q, Aa), sQC = pkpm_harm_dot(Q, Ac); 
    // 3x3 system A [c_0, alpha_0, c] = rhs: stage-A rates + interior stage-B rates + edge rates = 0 
    const double A00 = n0_q, A01 = eF0, A02 = eF1; 
    const double A10 = n1_q, A11 = n0_q + eF1, A12 = n1_q + eF2; 
    const double A20 = n2_q + 2.0*nG_q - fac*sQQ, A21 = 2.0*n1_q + eF2 + 2.0*eG0 - fac*sQA, A22 = 2.0*(n2_q + 2.0*nG_q) + eF3 + 2.0*eG1 - fac*(sQQ + sQC); 
    const double b0 = -r0_q - KF0, b1 = -r1_q - KF1, b2 = -r2_q - KF2 - 2.0*KG0 + fac*sQP; 
    const double det = A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20); 
    // conditioning proxy: the determinant relative to its edge-free value 
    // n0 (2 n0 (M2 + 2 M0G) - 2 M1^2) (positive for admissible moments); a small or negative 
    // ratio means the traces at the velocity-domain edges are not small 
    const double det0 = n0_q*(2.0*n0_q*(n2_q + 2.0*nG_q) - 2.0*n1_q*n1_q); 
    if (det0 > 0.0) cond_min = fmin(cond_min, det/det0); 
    if (n0_q > 0.0 && det > 0.0) { 
      c0_s[2] = (b0*(A11*A22 - A12*A21) - A01*(b1*A22 - A12*b2) + A02*(b1*A21 - A11*b2))/det; 
      a0_s[2] = (A00*(b1*A22 - A12*b2) - b0*(A10*A22 - A12*A20) + A02*(A10*b2 - b1*A20))/det; 
      c_s[2] = (A00*(A11*b2 - b1*A21) - A01*(A10*b2 - b1*A20) + b0*(A10*A21 - A11*A20))/det; 
      for (int a = 0; a < 3; ++a) 
        ap_s[a][2] = -(P[a] + (c0_s[2] + c_s[2])*Q[a] + Aa[a]*a0_s[2] + Ac[a]*c_s[2])/n0_q; 
      for (int c = 0; c < 8; ++c) { 
        flo_s[c][2] = (a0_s[2] + c_s[2]*w_lo)*glo[c] + Klo[c]; 
        fup_s[c][2] = (a0_s[2] + c_s[2]*w_up)*gup[c] + Kup[c]; 
      } 
    } else { 
      num_null += 1; 
    } 
  } 
  { 
    const double n0_q = 0.5*n0[3]+0.5*n0[2]+0.5*n0[1]+0.5*n0[0]; 
    const double n1_q = 0.5*n1[3]+0.5*n1[2]+0.5*n1[1]+0.5*n1[0]; 
    const double n2_q = 0.5*n2[3]+0.5*n2[2]+0.5*n2[1]+0.5*n2[0]; 
    const double nG_q = 0.5*nG[3]+0.5*nG[2]+0.5*nG[1]+0.5*nG[0]; 
    const double r0_q = 0.5*r0[3]+0.5*r0[2]+0.5*r0[1]+0.5*r0[0]; 
    const double r1_q = 0.5*r1[3]+0.5*r1[2]+0.5*r1[1]+0.5*r1[0]; 
    const double r2_q = 0.5*r2[3]+0.5*r2[2]+0.5*r2[1]+0.5*r2[0] + 2.0*(0.5*rG[3]+0.5*rG[2]+0.5*rG[1]+0.5*rG[0]); 
    const double mu_q[3] = {0.5*mu0[3]+0.5*mu0[2]+0.5*mu0[1]+0.5*mu0[0], 0.5*mu1[3]+0.5*mu1[2]+0.5*mu1[1]+0.5*mu1[0], 0.5*mu2[3]+0.5*mu2[2]+0.5*mu2[1]+0.5*mu2[0]}; 
    const double rM_q[3] = {0.5*rM0[3]+0.5*rM0[2]+0.5*rM0[1]+0.5*rM0[0], 0.5*rM1[3]+0.5*rM1[2]+0.5*rM1[1]+0.5*rM1[0], 0.5*rM2[3]+0.5*rM2[2]+0.5*rM2[1]+0.5*rM2[0]}; 
    const double b_q[3] = {0.5*b0[3]+0.5*b0[2]+0.5*b0[1]+0.5*b0[0], 0.5*b1[3]+0.5*b1[2]+0.5*b1[1]+0.5*b1[0], 0.5*b2[3]+0.5*b2[2]+0.5*b2[1]+0.5*b2[0]}; 
    // the p=1 projections of the edge fluxes' parts at the node 
    double Klo[8], Kup[8], glo[8], gup[8]; 
    for (int c = 0; c < 8; ++c) { 
      const double *h = &Kh[0][c*4]; 
      Klo[c] = 0.5*h[3]+0.5*h[2]+0.5*h[1]+0.5*h[0]; 
      h = &Kh[1][c*4]; 
      Kup[c] = 0.5*h[3]+0.5*h[2]+0.5*h[1]+0.5*h[0]; 
      h = &gh[0][c*4]; 
      glo[c] = 0.5*h[3]+0.5*h[2]+0.5*h[1]+0.5*h[0]; 
      h = &gh[1][c*4]; 
      gup[c] = 0.5*h[3]+0.5*h[2]+0.5*h[1]+0.5*h[0]; 
    } 
    // edge contributions to the moment rates: the k-th moment gains sgn_e w_e^k Phi_e per edge, 
    // sgn_e = +1 lo, -1 up: closure part (alpha_0 + c w_e) g_e and known part K_e 
    const double eF0 = glo[0] - gup[0]; 
    const double eF1 = w_lo*glo[0] - w_up*gup[0]; 
    const double eF2 = w_lo*w_lo*glo[0] - w_up*w_up*gup[0]; 
    const double eF3 = w_lo*w_lo*w_lo*glo[0] - w_up*w_up*w_up*gup[0]; 
    const double eG0 = glo[1] - gup[1]; 
    const double eG1 = w_lo*glo[1] - w_up*gup[1]; 
    const double KF0 = Klo[0] - Kup[0]; 
    const double KF1 = w_lo*Klo[0] - w_up*Kup[0]; 
    const double KF2 = w_lo*w_lo*Klo[0] - w_up*w_up*Kup[0]; 
    const double KG0 = Klo[1] - Kup[1]; 
    // perpendicular rows M0 a_perp + (c_0 + c) Q + alpha_0 A_a + c A_c = -P (b is a unit vector here) 
    double P[3], Q[3], Aa[3], Ac[3], v[3]; 
    for (int a = 0; a < 3; ++a) v[a] = rM_q[a] + Klo[2 + a] - Kup[2 + a]; 
    pkpm_harm_perp(b_q, v, P); 
    pkpm_harm_perp(b_q, mu_q, Q); 
    for (int a = 0; a < 3; ++a) v[a] = glo[2 + a] - gup[2 + a]; 
    pkpm_harm_perp(b_q, v, Aa); 
    for (int a = 0; a < 3; ++a) v[a] = w_lo*glo[2 + a] - w_up*gup[2 + a]; 
    pkpm_harm_perp(b_q, v, Ac); 
    // eliminated into the energy row, which carries 2 mu . a_perp (Schur complement) 
    const double fac = n0_q > 0.0 ? 2.0/n0_q : 0.0; 
    const double sQQ = pkpm_harm_dot(Q, Q), sQP = pkpm_harm_dot(Q, P), sQA = pkpm_harm_dot(Q, Aa), sQC = pkpm_harm_dot(Q, Ac); 
    // 3x3 system A [c_0, alpha_0, c] = rhs: stage-A rates + interior stage-B rates + edge rates = 0 
    const double A00 = n0_q, A01 = eF0, A02 = eF1; 
    const double A10 = n1_q, A11 = n0_q + eF1, A12 = n1_q + eF2; 
    const double A20 = n2_q + 2.0*nG_q - fac*sQQ, A21 = 2.0*n1_q + eF2 + 2.0*eG0 - fac*sQA, A22 = 2.0*(n2_q + 2.0*nG_q) + eF3 + 2.0*eG1 - fac*(sQQ + sQC); 
    const double b0 = -r0_q - KF0, b1 = -r1_q - KF1, b2 = -r2_q - KF2 - 2.0*KG0 + fac*sQP; 
    const double det = A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20); 
    // conditioning proxy: the determinant relative to its edge-free value 
    // n0 (2 n0 (M2 + 2 M0G) - 2 M1^2) (positive for admissible moments); a small or negative 
    // ratio means the traces at the velocity-domain edges are not small 
    const double det0 = n0_q*(2.0*n0_q*(n2_q + 2.0*nG_q) - 2.0*n1_q*n1_q); 
    if (det0 > 0.0) cond_min = fmin(cond_min, det/det0); 
    if (n0_q > 0.0 && det > 0.0) { 
      c0_s[3] = (b0*(A11*A22 - A12*A21) - A01*(b1*A22 - A12*b2) + A02*(b1*A21 - A11*b2))/det; 
      a0_s[3] = (A00*(b1*A22 - A12*b2) - b0*(A10*A22 - A12*A20) + A02*(A10*b2 - b1*A20))/det; 
      c_s[3] = (A00*(A11*b2 - b1*A21) - A01*(A10*b2 - b1*A20) + b0*(A10*A21 - A11*A20))/det; 
      for (int a = 0; a < 3; ++a) 
        ap_s[a][3] = -(P[a] + (c0_s[3] + c_s[3])*Q[a] + Aa[a]*a0_s[3] + Ac[a]*c_s[3])/n0_q; 
      for (int c = 0; c < 8; ++c) { 
        flo_s[c][3] = (a0_s[3] + c_s[3]*w_lo)*glo[c] + Klo[c]; 
        fup_s[c][3] = (a0_s[3] + c_s[3]*w_up)*gup[c] + Kup[c]; 
      } 
    } else { 
      num_null += 1; 
    } 
  } 
 
  cond[0] = cond_min; 
  if (num_null > 0) { 
    // a node with M0 <= 0 or a non-positive determinant zeroes everything in the cell (the 
    // closures it would have needed are not applied); the rate is that of the known w fluxes 
    double rate = 0.0; 
    for (int n = 0; n < 25; ++n) { 
      rate = fmax(rate, fabs(kw[n])*wm + rest[n]); 
    } 
    cfl[0] = 2.5*dv1*rate; cfl[1] = 0.0; 
    for (int k = 0; k < 4; ++k) { 
      c_0[k] = 0.0; alpha_0[k] = 0.0; c_h[k] = 0.0; 
    } 
    for (int k = 0; k < 12; ++k) { 
      a_perp[k] = 0.0; 
    } 
    for (int k = 0; k < 32; ++k) { 
      flux_lo[k] = 0.0; flux_up[k] = 0.0; 
    } 
    return num_null; 
  } 
 
  c_0[0] = 0.5*c0_s[3]+0.5*c0_s[2]+0.5*c0_s[1]+0.5*c0_s[0]; 
  c_0[1] = 0.5*c0_s[3]+0.5*c0_s[2]-0.5*c0_s[1]-0.5*c0_s[0]; 
  c_0[2] = 0.5*c0_s[3]-0.5*c0_s[2]+0.5*c0_s[1]-0.5*c0_s[0]; 
  c_0[3] = 0.5*c0_s[3]-0.5*c0_s[2]-0.5*c0_s[1]+0.5*c0_s[0]; 
  alpha_0[0] = 0.5*a0_s[3]+0.5*a0_s[2]+0.5*a0_s[1]+0.5*a0_s[0]; 
  alpha_0[1] = 0.5*a0_s[3]+0.5*a0_s[2]-0.5*a0_s[1]-0.5*a0_s[0]; 
  alpha_0[2] = 0.5*a0_s[3]-0.5*a0_s[2]+0.5*a0_s[1]-0.5*a0_s[0]; 
  alpha_0[3] = 0.5*a0_s[3]-0.5*a0_s[2]-0.5*a0_s[1]+0.5*a0_s[0]; 
  c_h[0] = 0.5*c_s[3]+0.5*c_s[2]+0.5*c_s[1]+0.5*c_s[0]; 
  c_h[1] = 0.5*c_s[3]+0.5*c_s[2]-0.5*c_s[1]-0.5*c_s[0]; 
  c_h[2] = 0.5*c_s[3]-0.5*c_s[2]+0.5*c_s[1]-0.5*c_s[0]; 
  c_h[3] = 0.5*c_s[3]-0.5*c_s[2]-0.5*c_s[1]+0.5*c_s[0]; 
 
  for (int a = 0; a < 3; ++a) { 
    a_perp[a*4 + 0] = (0.5)*ap_s[a][0] + (0.5)*ap_s[a][1] + (0.5)*ap_s[a][2] + (0.5)*ap_s[a][3]; 
    a_perp[a*4 + 1] = (-0.5)*ap_s[a][0] + (-0.5)*ap_s[a][1] + (0.5)*ap_s[a][2] + (0.5)*ap_s[a][3]; 
    a_perp[a*4 + 2] = (-0.5)*ap_s[a][0] + (0.5)*ap_s[a][1] + (-0.5)*ap_s[a][2] + (0.5)*ap_s[a][3]; 
    a_perp[a*4 + 3] = (0.5)*ap_s[a][0] + (-0.5)*ap_s[a][1] + (-0.5)*ap_s[a][2] + (0.5)*ap_s[a][3]; 
  } 
  for (int c = 0; c < 8; ++c) { 
    flux_lo[c*4 + 0] = (1.0)*flo_s[c][0] + (1.0)*flo_s[c][1] + (1.0)*flo_s[c][2] + (1.0)*flo_s[c][3]; 
    flux_up[c*4 + 0] = (1.0)*fup_s[c][0] + (1.0)*fup_s[c][1] + (1.0)*fup_s[c][2] + (1.0)*fup_s[c][3]; 
    flux_lo[c*4 + 1] = (-1.4433756729740645)*flo_s[c][0] + (-1.4433756729740645)*flo_s[c][1] + (1.4433756729740645)*flo_s[c][2] + (1.4433756729740645)*flo_s[c][3]; 
    flux_up[c*4 + 1] = (-1.4433756729740645)*fup_s[c][0] + (-1.4433756729740645)*fup_s[c][1] + (1.4433756729740645)*fup_s[c][2] + (1.4433756729740645)*fup_s[c][3]; 
    flux_lo[c*4 + 2] = (-1.4433756729740645)*flo_s[c][0] + (1.4433756729740645)*flo_s[c][1] + (-1.4433756729740645)*flo_s[c][2] + (1.4433756729740645)*flo_s[c][3]; 
    flux_up[c*4 + 2] = (-1.4433756729740645)*fup_s[c][0] + (1.4433756729740645)*fup_s[c][1] + (-1.4433756729740645)*fup_s[c][2] + (1.4433756729740645)*fup_s[c][3]; 
    flux_lo[c*4 + 3] = (2.0833333333333335)*flo_s[c][0] + (-2.0833333333333335)*flo_s[c][1] + (-2.0833333333333335)*flo_s[c][2] + (2.0833333333333335)*flo_s[c][3]; 
    flux_up[c*4 + 3] = (2.0833333333333335)*fup_s[c][0] + (-2.0833333333333335)*fup_s[c][1] + (-2.0833333333333335)*fup_s[c][2] + (2.0833333333333335)*fup_s[c][3]; 
  } 
 
  // the column time-step rate with the new closures (plan 17.27) 
  double a0f[25], cf[25]; 
  pkpm_vars_closure_2x_ser_p1_p1(alpha_0, a0f); 
  pkpm_vars_closure_2x_ser_p1_p1(c_h, cf); 
  double rate = 0.0; 
  for (int n = 0; n < 25; ++n) { 
    rate = fmax(rate, fabs(a0f[n]) + fabs(cf[n] + kw[n])*wm + rest[n]); 
  } 
  double src = 0.0; 
  src += 0.5*(fabs(c_0[0]) + 3.0*fabs(c_h[0])); 
  src += 0.8660254037844386*(fabs(c_0[1]) + 3.0*fabs(c_h[1])); 
  src += 0.8660254037844386*(fabs(c_0[2]) + 3.0*fabs(c_h[2])); 
  src += 1.5*(fabs(c_0[3]) + 3.0*fabs(c_h[3])); 
  cfl[0] = 2.5*dv1*rate + src; 
  cfl[1] = 0.0; 
 
  return 0; 
} 
