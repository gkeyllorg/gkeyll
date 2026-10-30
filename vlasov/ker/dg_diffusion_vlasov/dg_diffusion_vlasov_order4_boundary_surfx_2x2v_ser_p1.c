#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order4_boundary_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[36] = {0.0}; 

  double edgeSurf_incr[36] = {0.0}; 
  double boundSurf_incr[36] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.6237976320958223*coeff[0]*fSkin[1]+1.6237976320958223*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]+1.6237976320958223*coeff[0]*fSkin[0]-1.6237976320958223*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 1.6237976320958223*coeff[0]*fSkin[5]+1.6237976320958223*coeff[0]*fEdge[5]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 1.6237976320958223*coeff[0]*fSkin[6]+1.6237976320958223*coeff[0]*fEdge[6]+0.9375*coeff[0]*fSkin[3]-0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 1.6237976320958223*coeff[0]*fSkin[8]+1.6237976320958223*coeff[0]*fEdge[8]+0.9375*coeff[0]*fSkin[4]-0.9375*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*fSkin[5]+2.0625*coeff[0]*fEdge[5]+1.6237976320958223*coeff[0]*fSkin[2]-1.6237976320958223*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = 3.5625*coeff[0]*fSkin[6]+2.0625*coeff[0]*fEdge[6]+1.6237976320958223*coeff[0]*fSkin[3]-1.6237976320958223*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 1.6237976320958223*coeff[0]*fSkin[13]+1.6237976320958223*coeff[0]*fEdge[13]+0.9375*coeff[0]*fSkin[7]-0.9375*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = 3.5625*coeff[0]*fSkin[8]+2.0625*coeff[0]*fEdge[8]+1.6237976320958223*coeff[0]*fSkin[4]-1.6237976320958223*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 1.6237976320958223*coeff[0]*fSkin[14]+1.6237976320958223*coeff[0]*fEdge[14]+0.9375*coeff[0]*fSkin[9]-0.9375*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 1.6237976320958223*coeff[0]*fSkin[15]+1.6237976320958223*coeff[0]*fEdge[15]+0.9375*coeff[0]*fSkin[10]-0.9375*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 1.6237976320958225*coeff[0]*fSkin[17]+1.6237976320958225*coeff[0]*fEdge[17]+0.9375*coeff[0]*fSkin[11]-0.9375*coeff[0]*fEdge[11]; 
  edgeSurf_incr[12] = 1.6237976320958225*coeff[0]*fSkin[20]+1.6237976320958225*coeff[0]*fEdge[20]+0.9375*coeff[0]*fSkin[12]-0.9375*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = 3.5625*coeff[0]*fSkin[13]+2.0625*coeff[0]*fEdge[13]+1.6237976320958223*coeff[0]*fSkin[7]-1.6237976320958223*coeff[0]*fEdge[7]; 
  edgeSurf_incr[14] = 3.5625*coeff[0]*fSkin[14]+2.0625*coeff[0]*fEdge[14]+1.6237976320958223*coeff[0]*fSkin[9]-1.6237976320958223*coeff[0]*fEdge[9]; 
  edgeSurf_incr[15] = 3.5625*coeff[0]*fSkin[15]+2.0625*coeff[0]*fEdge[15]+1.6237976320958223*coeff[0]*fSkin[10]-1.6237976320958223*coeff[0]*fEdge[10]; 
  edgeSurf_incr[16] = 1.6237976320958223*coeff[0]*fSkin[23]+1.6237976320958223*coeff[0]*fEdge[23]+0.9375*coeff[0]*fSkin[16]-0.9375*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 3.5625*coeff[0]*fSkin[17]+2.0625*coeff[0]*fEdge[17]+1.6237976320958225*coeff[0]*fSkin[11]-1.6237976320958225*coeff[0]*fEdge[11]; 
  edgeSurf_incr[18] = 1.6237976320958225*coeff[0]*fSkin[24]+1.6237976320958225*coeff[0]*fEdge[24]+0.9375*coeff[0]*fSkin[18]-0.9375*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = 1.6237976320958225*coeff[0]*fSkin[25]+1.6237976320958225*coeff[0]*fEdge[25]+0.9375*coeff[0]*fSkin[19]-0.9375*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = 3.5625*coeff[0]*fSkin[20]+2.0625*coeff[0]*fEdge[20]+1.6237976320958225*coeff[0]*fSkin[12]-1.6237976320958225*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = 1.6237976320958225*coeff[0]*fSkin[27]+1.6237976320958225*coeff[0]*fEdge[27]+0.9375*coeff[0]*fSkin[21]-0.9375*coeff[0]*fEdge[21]; 
  edgeSurf_incr[22] = 1.6237976320958225*coeff[0]*fSkin[28]+1.6237976320958225*coeff[0]*fEdge[28]+0.9375*coeff[0]*fSkin[22]-0.9375*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 3.5625*coeff[0]*fSkin[23]+2.0625*coeff[0]*fEdge[23]+1.6237976320958223*coeff[0]*fSkin[16]-1.6237976320958223*coeff[0]*fEdge[16]; 
  edgeSurf_incr[24] = 3.5625*coeff[0]*fSkin[24]+2.0625*coeff[0]*fEdge[24]+1.6237976320958225*coeff[0]*fSkin[18]-1.6237976320958225*coeff[0]*fEdge[18]; 
  edgeSurf_incr[25] = 3.5625*coeff[0]*fSkin[25]+2.0625*coeff[0]*fEdge[25]+1.6237976320958225*coeff[0]*fSkin[19]-1.6237976320958225*coeff[0]*fEdge[19]; 
  edgeSurf_incr[26] = 1.6237976320958225*coeff[0]*fSkin[31]+1.6237976320958225*coeff[0]*fEdge[31]+0.9375*coeff[0]*fSkin[26]-0.9375*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = 3.5625*coeff[0]*fSkin[27]+2.0625*coeff[0]*fEdge[27]+1.6237976320958225*coeff[0]*fSkin[21]-1.6237976320958225*coeff[0]*fEdge[21]; 
  edgeSurf_incr[28] = 3.5625*coeff[0]*fSkin[28]+2.0625*coeff[0]*fEdge[28]+1.6237976320958225*coeff[0]*fSkin[22]-1.6237976320958225*coeff[0]*fEdge[22]; 
  edgeSurf_incr[29] = 1.6237976320958225*coeff[0]*fSkin[32]+1.6237976320958225*coeff[0]*fEdge[32]+0.9375*coeff[0]*fSkin[29]-0.9375*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = 1.6237976320958223*coeff[0]*fSkin[33]+1.6237976320958223*coeff[0]*fEdge[33]+0.9375*coeff[0]*fSkin[30]-0.9375*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = 3.5625*coeff[0]*fSkin[31]+2.0625*coeff[0]*fEdge[31]+1.6237976320958225*coeff[0]*fSkin[26]-1.6237976320958225*coeff[0]*fEdge[26]; 
  edgeSurf_incr[32] = 3.5625*coeff[0]*fSkin[32]+2.0625*coeff[0]*fEdge[32]+1.6237976320958225*coeff[0]*fSkin[29]-1.6237976320958225*coeff[0]*fEdge[29]; 
  edgeSurf_incr[33] = 3.5625*coeff[0]*fSkin[33]+2.0625*coeff[0]*fEdge[33]+1.6237976320958223*coeff[0]*fSkin[30]-1.6237976320958223*coeff[0]*fEdge[30]; 
  edgeSurf_incr[34] = 1.6237976320958223*coeff[0]*fSkin[35]+1.6237976320958223*coeff[0]*fEdge[35]+0.9375*coeff[0]*fSkin[34]-0.9375*coeff[0]*fEdge[34]; 
  edgeSurf_incr[35] = 3.5625*coeff[0]*fSkin[35]+2.0625*coeff[0]*fEdge[35]+1.6237976320958223*coeff[0]*fSkin[34]-1.6237976320958223*coeff[0]*fEdge[34]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = 1.5*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 1.5*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = 1.5*coeff[0]*fSkin[8]; 
  boundSurf_incr[13] = 1.5*coeff[0]*fSkin[13]; 
  boundSurf_incr[14] = 1.5*coeff[0]*fSkin[14]; 
  boundSurf_incr[15] = 1.5*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 1.5*coeff[0]*fSkin[17]; 
  boundSurf_incr[20] = 1.5*coeff[0]*fSkin[20]; 
  boundSurf_incr[23] = 1.5*coeff[0]*fSkin[23]; 
  boundSurf_incr[24] = 1.5*coeff[0]*fSkin[24]; 
  boundSurf_incr[25] = 1.5*coeff[0]*fSkin[25]; 
  boundSurf_incr[27] = 1.5*coeff[0]*fSkin[27]; 
  boundSurf_incr[28] = 1.5*coeff[0]*fSkin[28]; 
  boundSurf_incr[31] = 1.5*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = 1.5*coeff[0]*fSkin[32]; 
  boundSurf_incr[33] = 1.5*coeff[0]*fSkin[33]; 
  boundSurf_incr[35] = 1.5*coeff[0]*fSkin[35]; 

  } else { 

  edgeSurf_incr[0] = -(1.6237976320958223*coeff[0]*fSkin[1])-1.6237976320958223*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]-1.6237976320958223*coeff[0]*fSkin[0]+1.6237976320958223*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(1.6237976320958223*coeff[0]*fSkin[5])-1.6237976320958223*coeff[0]*fEdge[5]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(1.6237976320958223*coeff[0]*fSkin[6])-1.6237976320958223*coeff[0]*fEdge[6]+0.9375*coeff[0]*fSkin[3]-0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(1.6237976320958223*coeff[0]*fSkin[8])-1.6237976320958223*coeff[0]*fEdge[8]+0.9375*coeff[0]*fSkin[4]-0.9375*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*fSkin[5]+2.0625*coeff[0]*fEdge[5]-1.6237976320958223*coeff[0]*fSkin[2]+1.6237976320958223*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = 3.5625*coeff[0]*fSkin[6]+2.0625*coeff[0]*fEdge[6]-1.6237976320958223*coeff[0]*fSkin[3]+1.6237976320958223*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = -(1.6237976320958223*coeff[0]*fSkin[13])-1.6237976320958223*coeff[0]*fEdge[13]+0.9375*coeff[0]*fSkin[7]-0.9375*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = 3.5625*coeff[0]*fSkin[8]+2.0625*coeff[0]*fEdge[8]-1.6237976320958223*coeff[0]*fSkin[4]+1.6237976320958223*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = -(1.6237976320958223*coeff[0]*fSkin[14])-1.6237976320958223*coeff[0]*fEdge[14]+0.9375*coeff[0]*fSkin[9]-0.9375*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = -(1.6237976320958223*coeff[0]*fSkin[15])-1.6237976320958223*coeff[0]*fEdge[15]+0.9375*coeff[0]*fSkin[10]-0.9375*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = -(1.6237976320958225*coeff[0]*fSkin[17])-1.6237976320958225*coeff[0]*fEdge[17]+0.9375*coeff[0]*fSkin[11]-0.9375*coeff[0]*fEdge[11]; 
  edgeSurf_incr[12] = -(1.6237976320958225*coeff[0]*fSkin[20])-1.6237976320958225*coeff[0]*fEdge[20]+0.9375*coeff[0]*fSkin[12]-0.9375*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = 3.5625*coeff[0]*fSkin[13]+2.0625*coeff[0]*fEdge[13]-1.6237976320958223*coeff[0]*fSkin[7]+1.6237976320958223*coeff[0]*fEdge[7]; 
  edgeSurf_incr[14] = 3.5625*coeff[0]*fSkin[14]+2.0625*coeff[0]*fEdge[14]-1.6237976320958223*coeff[0]*fSkin[9]+1.6237976320958223*coeff[0]*fEdge[9]; 
  edgeSurf_incr[15] = 3.5625*coeff[0]*fSkin[15]+2.0625*coeff[0]*fEdge[15]-1.6237976320958223*coeff[0]*fSkin[10]+1.6237976320958223*coeff[0]*fEdge[10]; 
  edgeSurf_incr[16] = -(1.6237976320958223*coeff[0]*fSkin[23])-1.6237976320958223*coeff[0]*fEdge[23]+0.9375*coeff[0]*fSkin[16]-0.9375*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 3.5625*coeff[0]*fSkin[17]+2.0625*coeff[0]*fEdge[17]-1.6237976320958225*coeff[0]*fSkin[11]+1.6237976320958225*coeff[0]*fEdge[11]; 
  edgeSurf_incr[18] = -(1.6237976320958225*coeff[0]*fSkin[24])-1.6237976320958225*coeff[0]*fEdge[24]+0.9375*coeff[0]*fSkin[18]-0.9375*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = -(1.6237976320958225*coeff[0]*fSkin[25])-1.6237976320958225*coeff[0]*fEdge[25]+0.9375*coeff[0]*fSkin[19]-0.9375*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = 3.5625*coeff[0]*fSkin[20]+2.0625*coeff[0]*fEdge[20]-1.6237976320958225*coeff[0]*fSkin[12]+1.6237976320958225*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = -(1.6237976320958225*coeff[0]*fSkin[27])-1.6237976320958225*coeff[0]*fEdge[27]+0.9375*coeff[0]*fSkin[21]-0.9375*coeff[0]*fEdge[21]; 
  edgeSurf_incr[22] = -(1.6237976320958225*coeff[0]*fSkin[28])-1.6237976320958225*coeff[0]*fEdge[28]+0.9375*coeff[0]*fSkin[22]-0.9375*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 3.5625*coeff[0]*fSkin[23]+2.0625*coeff[0]*fEdge[23]-1.6237976320958223*coeff[0]*fSkin[16]+1.6237976320958223*coeff[0]*fEdge[16]; 
  edgeSurf_incr[24] = 3.5625*coeff[0]*fSkin[24]+2.0625*coeff[0]*fEdge[24]-1.6237976320958225*coeff[0]*fSkin[18]+1.6237976320958225*coeff[0]*fEdge[18]; 
  edgeSurf_incr[25] = 3.5625*coeff[0]*fSkin[25]+2.0625*coeff[0]*fEdge[25]-1.6237976320958225*coeff[0]*fSkin[19]+1.6237976320958225*coeff[0]*fEdge[19]; 
  edgeSurf_incr[26] = -(1.6237976320958225*coeff[0]*fSkin[31])-1.6237976320958225*coeff[0]*fEdge[31]+0.9375*coeff[0]*fSkin[26]-0.9375*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = 3.5625*coeff[0]*fSkin[27]+2.0625*coeff[0]*fEdge[27]-1.6237976320958225*coeff[0]*fSkin[21]+1.6237976320958225*coeff[0]*fEdge[21]; 
  edgeSurf_incr[28] = 3.5625*coeff[0]*fSkin[28]+2.0625*coeff[0]*fEdge[28]-1.6237976320958225*coeff[0]*fSkin[22]+1.6237976320958225*coeff[0]*fEdge[22]; 
  edgeSurf_incr[29] = -(1.6237976320958225*coeff[0]*fSkin[32])-1.6237976320958225*coeff[0]*fEdge[32]+0.9375*coeff[0]*fSkin[29]-0.9375*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = -(1.6237976320958223*coeff[0]*fSkin[33])-1.6237976320958223*coeff[0]*fEdge[33]+0.9375*coeff[0]*fSkin[30]-0.9375*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = 3.5625*coeff[0]*fSkin[31]+2.0625*coeff[0]*fEdge[31]-1.6237976320958225*coeff[0]*fSkin[26]+1.6237976320958225*coeff[0]*fEdge[26]; 
  edgeSurf_incr[32] = 3.5625*coeff[0]*fSkin[32]+2.0625*coeff[0]*fEdge[32]-1.6237976320958225*coeff[0]*fSkin[29]+1.6237976320958225*coeff[0]*fEdge[29]; 
  edgeSurf_incr[33] = 3.5625*coeff[0]*fSkin[33]+2.0625*coeff[0]*fEdge[33]-1.6237976320958223*coeff[0]*fSkin[30]+1.6237976320958223*coeff[0]*fEdge[30]; 
  edgeSurf_incr[34] = -(1.6237976320958223*coeff[0]*fSkin[35])-1.6237976320958223*coeff[0]*fEdge[35]+0.9375*coeff[0]*fSkin[34]-0.9375*coeff[0]*fEdge[34]; 
  edgeSurf_incr[35] = 3.5625*coeff[0]*fSkin[35]+2.0625*coeff[0]*fEdge[35]-1.6237976320958223*coeff[0]*fSkin[34]+1.6237976320958223*coeff[0]*fEdge[34]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = 1.5*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 1.5*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = 1.5*coeff[0]*fSkin[8]; 
  boundSurf_incr[13] = 1.5*coeff[0]*fSkin[13]; 
  boundSurf_incr[14] = 1.5*coeff[0]*fSkin[14]; 
  boundSurf_incr[15] = 1.5*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 1.5*coeff[0]*fSkin[17]; 
  boundSurf_incr[20] = 1.5*coeff[0]*fSkin[20]; 
  boundSurf_incr[23] = 1.5*coeff[0]*fSkin[23]; 
  boundSurf_incr[24] = 1.5*coeff[0]*fSkin[24]; 
  boundSurf_incr[25] = 1.5*coeff[0]*fSkin[25]; 
  boundSurf_incr[27] = 1.5*coeff[0]*fSkin[27]; 
  boundSurf_incr[28] = 1.5*coeff[0]*fSkin[28]; 
  boundSurf_incr[31] = 1.5*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = 1.5*coeff[0]*fSkin[32]; 
  boundSurf_incr[33] = 1.5*coeff[0]*fSkin[33]; 
  boundSurf_incr[35] = 1.5*coeff[0]*fSkin[35]; 

  }

  out[0] += -(1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac); 
  out[1] += -(1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac); 
  out[2] += -(1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac); 
  out[3] += -(1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac); 
  out[4] += -(1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac); 
  out[5] += -(1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac); 
  out[6] += -(1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac); 
  out[7] += -(1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac); 
  out[8] += -(1.0*(vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac); 
  out[9] += -(1.0*(vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac); 
  out[10] += -(1.0*(vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac); 
  out[11] += -(1.0*(vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac); 
  out[12] += -(1.0*(vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac); 
  out[13] += -(1.0*(vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac); 
  out[14] += -(1.0*(vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac); 
  out[15] += -(1.0*(vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac); 
  out[16] += -(1.0*(vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac); 
  out[17] += -(1.0*(vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac); 
  out[18] += -(1.0*(vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac); 
  out[19] += -(1.0*(vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac); 
  out[20] += -(1.0*(vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac); 
  out[21] += -(1.0*(vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac); 
  out[22] += -(1.0*(vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac); 
  out[23] += -(1.0*(vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac); 
  out[24] += -(1.0*(vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac); 
  out[25] += -(1.0*(vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac); 
  out[26] += -(1.0*(vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac); 
  out[27] += -(1.0*(vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*Jfac); 
  out[28] += -(1.0*(vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*Jfac); 
  out[29] += -(1.0*(vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*Jfac); 
  out[30] += -(1.0*(vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*Jfac); 
  out[31] += -(1.0*(vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*Jfac); 
  out[32] += -(1.0*(vol_incr[32]+edgeSurf_incr[32]+boundSurf_incr[32])*Jfac); 
  out[33] += -(1.0*(vol_incr[33]+edgeSurf_incr[33]+boundSurf_incr[33])*Jfac); 
  out[34] += -(1.0*(vol_incr[34]+edgeSurf_incr[34]+boundSurf_incr[34])*Jfac); 
  out[35] += -(1.0*(vol_incr[35]+edgeSurf_incr[35]+boundSurf_incr[35])*Jfac); 

  return 0.;
}

