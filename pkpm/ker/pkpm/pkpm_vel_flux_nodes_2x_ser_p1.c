#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 

// The p=1 field f (gkyl basis) at the 25 Gauss-Legendre nodes (5 per direction, row-major: the 
// last direction fastest), direction by direction. 
GKYL_CU_DH static inline void pkpm_vel_flux_nodes_2x_ser_p1_p1(const double *f, double* GKYL_RESTRICT fn) 
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

GKYL_CU_DH void pkpm_vel_flux_nodes_2x_ser_p1(const double *wn, const double *theta, const double *bvar, const double *alpha_0,
  const double *c_h, double* GKYL_RESTRICT nv) 
{ 
  // wn:      Raw w-flux values at the 25 nodes [node][PKPM_HS_W_END] (the build of the stored operator). 
  // theta:   T_perp/T. 
  // bvar:    Magnetic field unit vector [bx, by, bz]. 
  // alpha_0: Pressure force a_P/xi. 
  // c_h:     Heating rate. 
  // nv:      Output node values [25 nodes][PKPM_HS_NODE_NV]. 
  double th_n[25]; 
  pkpm_vel_flux_nodes_2x_ser_p1_p1(theta, th_n); 
  double b0_n[25]; 
  pkpm_vel_flux_nodes_2x_ser_p1_p1(&bvar[0], b0_n); 
  double b1_n[25]; 
  pkpm_vel_flux_nodes_2x_ser_p1_p1(&bvar[4], b1_n); 
  double b2_n[25]; 
  pkpm_vel_flux_nodes_2x_ser_p1_p1(&bvar[8], b2_n); 
  double a0_n[25]; 
  pkpm_vel_flux_nodes_2x_ser_p1_p1(alpha_0, a0_n); 
  double ch_n[25]; 
  pkpm_vel_flux_nodes_2x_ser_p1_p1(c_h, ch_n); 
  for (int n = 0; n < 25; ++n) { 
    double *v = &nv[n*PKPM_HS_NODE_NV]; 
    for (int s = 0; s < PKPM_HS_W_END; ++s) v[s] = wn[n*PKPM_HS_W_END + s]; 
    v[PKPM_HS_NODE_THETA] = th_n[n]; 
    v[PKPM_HS_NODE_B] = b0_n[n]; v[PKPM_HS_NODE_B + 1] = b1_n[n]; v[PKPM_HS_NODE_B + 2] = b2_n[n]; 
    v[PKPM_HS_NODE_A0] = a0_n[n]; v[PKPM_HS_NODE_CH] = ch_n[n]; 
    pkpm_harm_node_spec(v); 
  } 
} 
