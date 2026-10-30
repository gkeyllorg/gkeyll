#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order2_boundary_surfx_1x3v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[54] = {0.0}; 

  double edgeSurf_incr[54] = {0.0}; 
  double boundSurf_incr[54] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.5412658773652741*coeff[0]*fSkin[1])-0.5412658773652741*coeff[0]*fEdge[1]-0.5625*coeff[0]*fSkin[0]+0.5625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.4375*coeff[0]*fSkin[1])-0.4375*coeff[0]*fEdge[1]-1.4072912811497125*coeff[0]*fSkin[0]+0.5412658773652739*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(0.5412658773652741*coeff[0]*fSkin[5])-0.5412658773652741*coeff[0]*fEdge[5]-0.5625*coeff[0]*fSkin[2]+0.5625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.5412658773652741*coeff[0]*fSkin[6])-0.5412658773652741*coeff[0]*fEdge[6]-0.5625*coeff[0]*fSkin[3]+0.5625*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(0.5412658773652741*coeff[0]*fSkin[8])-0.5412658773652741*coeff[0]*fEdge[8]-0.5625*coeff[0]*fSkin[4]+0.5625*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = -(1.4375*coeff[0]*fSkin[5])-0.4375*coeff[0]*fEdge[5]-1.4072912811497125*coeff[0]*fSkin[2]+0.5412658773652739*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = -(1.4375*coeff[0]*fSkin[6])-0.4375*coeff[0]*fEdge[6]-1.4072912811497125*coeff[0]*fSkin[3]+0.5412658773652739*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = -(0.5412658773652741*coeff[0]*fSkin[14])-0.5412658773652741*coeff[0]*fEdge[14]-0.5625*coeff[0]*fSkin[7]+0.5625*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = -(1.4375*coeff[0]*fSkin[8])-0.4375*coeff[0]*fEdge[8]-1.4072912811497125*coeff[0]*fSkin[4]+0.5412658773652739*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = -(0.5412658773652741*coeff[0]*fSkin[15])-0.5412658773652741*coeff[0]*fEdge[15]-0.5625*coeff[0]*fSkin[9]+0.5625*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = -(0.5412658773652741*coeff[0]*fSkin[16])-0.5412658773652741*coeff[0]*fEdge[16]-0.5625*coeff[0]*fSkin[10]+0.5625*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = -(0.5412658773652742*coeff[0]*fSkin[18])-0.5412658773652742*coeff[0]*fEdge[18]-0.5625*coeff[0]*fSkin[11]+0.5625*coeff[0]*fEdge[11]; 
  edgeSurf_incr[12] = -(0.5412658773652742*coeff[0]*fSkin[20])-0.5412658773652742*coeff[0]*fEdge[20]-0.5625*coeff[0]*fSkin[12]+0.5625*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = -(0.5412658773652742*coeff[0]*fSkin[24])-0.5412658773652742*coeff[0]*fEdge[24]-0.5625*coeff[0]*fSkin[13]+0.5625*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = -(1.4375*coeff[0]*fSkin[14])-0.4375*coeff[0]*fEdge[14]-1.4072912811497125*coeff[0]*fSkin[7]+0.5412658773652739*coeff[0]*fEdge[7]; 
  edgeSurf_incr[15] = -(1.4375*coeff[0]*fSkin[15])-0.4375*coeff[0]*fEdge[15]-1.4072912811497125*coeff[0]*fSkin[9]+0.5412658773652739*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = -(1.4375*coeff[0]*fSkin[16])-0.4375*coeff[0]*fEdge[16]-1.4072912811497125*coeff[0]*fSkin[10]+0.5412658773652739*coeff[0]*fEdge[10]; 
  edgeSurf_incr[17] = -(0.5412658773652741*coeff[0]*fSkin[27])-0.5412658773652741*coeff[0]*fEdge[27]-0.5625*coeff[0]*fSkin[17]+0.5625*coeff[0]*fEdge[17]; 
  edgeSurf_incr[18] = -(1.4375*coeff[0]*fSkin[18])-0.4375*coeff[0]*fEdge[18]-1.4072912811497127*coeff[0]*fSkin[11]+0.5412658773652742*coeff[0]*fEdge[11]; 
  edgeSurf_incr[19] = -(0.5412658773652742*coeff[0]*fSkin[28])-0.5412658773652742*coeff[0]*fEdge[28]-0.5625*coeff[0]*fSkin[19]+0.5625*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = -(1.4375*coeff[0]*fSkin[20])-0.4375*coeff[0]*fEdge[20]-1.4072912811497127*coeff[0]*fSkin[12]+0.5412658773652742*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = -(0.5412658773652742*coeff[0]*fSkin[29])-0.5412658773652742*coeff[0]*fEdge[29]-0.5625*coeff[0]*fSkin[21]+0.5625*coeff[0]*fEdge[21]; 
  edgeSurf_incr[22] = -(0.5412658773652742*coeff[0]*fSkin[30])-0.5412658773652742*coeff[0]*fEdge[30]-0.5625*coeff[0]*fSkin[22]+0.5625*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = -(0.5412658773652742*coeff[0]*fSkin[32])-0.5412658773652742*coeff[0]*fEdge[32]-0.5625*coeff[0]*fSkin[23]+0.5625*coeff[0]*fEdge[23]; 
  edgeSurf_incr[24] = -(1.4375*coeff[0]*fSkin[24])-0.4375*coeff[0]*fEdge[24]-1.4072912811497127*coeff[0]*fSkin[13]+0.5412658773652742*coeff[0]*fEdge[13]; 
  edgeSurf_incr[25] = -(0.5412658773652742*coeff[0]*fSkin[34])-0.5412658773652742*coeff[0]*fEdge[34]-0.5625*coeff[0]*fSkin[25]+0.5625*coeff[0]*fEdge[25]; 
  edgeSurf_incr[26] = -(0.5412658773652742*coeff[0]*fSkin[35])-0.5412658773652742*coeff[0]*fEdge[35]-0.5625*coeff[0]*fSkin[26]+0.5625*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = -(1.4375*coeff[0]*fSkin[27])-0.4375*coeff[0]*fEdge[27]-1.4072912811497125*coeff[0]*fSkin[17]+0.5412658773652739*coeff[0]*fEdge[17]; 
  edgeSurf_incr[28] = -(1.4375*coeff[0]*fSkin[28])-0.4375*coeff[0]*fEdge[28]-1.4072912811497127*coeff[0]*fSkin[19]+0.5412658773652742*coeff[0]*fEdge[19]; 
  edgeSurf_incr[29] = -(1.4375*coeff[0]*fSkin[29])-0.4375*coeff[0]*fEdge[29]-1.4072912811497127*coeff[0]*fSkin[21]+0.5412658773652742*coeff[0]*fEdge[21]; 
  edgeSurf_incr[30] = -(1.4375*coeff[0]*fSkin[30])-0.4375*coeff[0]*fEdge[30]-1.4072912811497127*coeff[0]*fSkin[22]+0.5412658773652742*coeff[0]*fEdge[22]; 
  edgeSurf_incr[31] = -(0.5412658773652742*coeff[0]*fSkin[40])-0.5412658773652742*coeff[0]*fEdge[40]-0.5625*coeff[0]*fSkin[31]+0.5625*coeff[0]*fEdge[31]; 
  edgeSurf_incr[32] = -(1.4375*coeff[0]*fSkin[32])-0.4375*coeff[0]*fEdge[32]-1.4072912811497127*coeff[0]*fSkin[23]+0.5412658773652742*coeff[0]*fEdge[23]; 
  edgeSurf_incr[33] = -(0.5412658773652742*coeff[0]*fSkin[41])-0.5412658773652742*coeff[0]*fEdge[41]-0.5625*coeff[0]*fSkin[33]+0.5625*coeff[0]*fEdge[33]; 
  edgeSurf_incr[34] = -(1.4375*coeff[0]*fSkin[34])-0.4375*coeff[0]*fEdge[34]-1.4072912811497127*coeff[0]*fSkin[25]+0.5412658773652742*coeff[0]*fEdge[25]; 
  edgeSurf_incr[35] = -(1.4375*coeff[0]*fSkin[35])-0.4375*coeff[0]*fEdge[35]-1.4072912811497127*coeff[0]*fSkin[26]+0.5412658773652742*coeff[0]*fEdge[26]; 
  edgeSurf_incr[36] = -(0.5412658773652742*coeff[0]*fSkin[42])-0.5412658773652742*coeff[0]*fEdge[42]-0.5625*coeff[0]*fSkin[36]+0.5625*coeff[0]*fEdge[36]; 
  edgeSurf_incr[37] = -(0.5412658773652741*coeff[0]*fSkin[43])-0.5412658773652741*coeff[0]*fEdge[43]-0.5625*coeff[0]*fSkin[37]+0.5625*coeff[0]*fEdge[37]; 
  edgeSurf_incr[38] = -(0.5412658773652741*coeff[0]*fSkin[45])-0.5412658773652741*coeff[0]*fEdge[45]-0.5625*coeff[0]*fSkin[38]+0.5625*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = -(0.5412658773652741*coeff[0]*fSkin[47])-0.5412658773652741*coeff[0]*fEdge[47]-0.5625*coeff[0]*fSkin[39]+0.5625*coeff[0]*fEdge[39]; 
  edgeSurf_incr[40] = -(1.4375*coeff[0]*fSkin[40])-0.4375*coeff[0]*fEdge[40]-1.4072912811497127*coeff[0]*fSkin[31]+0.5412658773652742*coeff[0]*fEdge[31]; 
  edgeSurf_incr[41] = -(1.4375*coeff[0]*fSkin[41])-0.4375*coeff[0]*fEdge[41]-1.4072912811497127*coeff[0]*fSkin[33]+0.5412658773652742*coeff[0]*fEdge[33]; 
  edgeSurf_incr[42] = -(1.4375*coeff[0]*fSkin[42])-0.4375*coeff[0]*fEdge[42]-1.4072912811497127*coeff[0]*fSkin[36]+0.5412658773652742*coeff[0]*fEdge[36]; 
  edgeSurf_incr[43] = -(1.4375*coeff[0]*fSkin[43])-0.4375*coeff[0]*fEdge[43]-1.4072912811497125*coeff[0]*fSkin[37]+0.5412658773652739*coeff[0]*fEdge[37]; 
  edgeSurf_incr[44] = -(0.5412658773652741*coeff[0]*fSkin[49])-0.5412658773652741*coeff[0]*fEdge[49]-0.5625*coeff[0]*fSkin[44]+0.5625*coeff[0]*fEdge[44]; 
  edgeSurf_incr[45] = -(1.4375*coeff[0]*fSkin[45])-0.4375*coeff[0]*fEdge[45]-1.4072912811497125*coeff[0]*fSkin[38]+0.5412658773652739*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = -(0.5412658773652741*coeff[0]*fSkin[50])-0.5412658773652741*coeff[0]*fEdge[50]-0.5625*coeff[0]*fSkin[46]+0.5625*coeff[0]*fEdge[46]; 
  edgeSurf_incr[47] = -(1.4375*coeff[0]*fSkin[47])-0.4375*coeff[0]*fEdge[47]-1.4072912811497125*coeff[0]*fSkin[39]+0.5412658773652739*coeff[0]*fEdge[39]; 
  edgeSurf_incr[48] = -(0.5412658773652741*coeff[0]*fSkin[51])-0.5412658773652741*coeff[0]*fEdge[51]-0.5625*coeff[0]*fSkin[48]+0.5625*coeff[0]*fEdge[48]; 
  edgeSurf_incr[49] = -(1.4375*coeff[0]*fSkin[49])-0.4375*coeff[0]*fEdge[49]-1.4072912811497125*coeff[0]*fSkin[44]+0.5412658773652739*coeff[0]*fEdge[44]; 
  edgeSurf_incr[50] = -(1.4375*coeff[0]*fSkin[50])-0.4375*coeff[0]*fEdge[50]-1.4072912811497125*coeff[0]*fSkin[46]+0.5412658773652739*coeff[0]*fEdge[46]; 
  edgeSurf_incr[51] = -(1.4375*coeff[0]*fSkin[51])-0.4375*coeff[0]*fEdge[51]-1.4072912811497125*coeff[0]*fSkin[48]+0.5412658773652739*coeff[0]*fEdge[48]; 
  edgeSurf_incr[52] = -(0.5412658773652742*coeff[0]*fSkin[53])-0.5412658773652742*coeff[0]*fEdge[53]-0.5625*coeff[0]*fSkin[52]+0.5625*coeff[0]*fEdge[52]; 
  edgeSurf_incr[53] = -(1.4375*coeff[0]*fSkin[53])-0.4375*coeff[0]*fEdge[53]-1.4072912811497127*coeff[0]*fSkin[52]+0.5412658773652742*coeff[0]*fEdge[52]; 

  boundSurf_incr[1] = 0.8660254037844386*coeff[0]*fSkin[0]-1.0*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = 0.8660254037844386*coeff[0]*fSkin[2]-1.0*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 0.8660254037844386*coeff[0]*fSkin[3]-1.0*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = 0.8660254037844386*coeff[0]*fSkin[4]-1.0*coeff[0]*fSkin[8]; 
  boundSurf_incr[14] = 0.8660254037844386*coeff[0]*fSkin[7]-1.0*coeff[0]*fSkin[14]; 
  boundSurf_incr[15] = 0.8660254037844386*coeff[0]*fSkin[9]-1.0*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = 0.8660254037844386*coeff[0]*fSkin[10]-1.0*coeff[0]*fSkin[16]; 
  boundSurf_incr[18] = 0.8660254037844387*coeff[0]*fSkin[11]-1.0*coeff[0]*fSkin[18]; 
  boundSurf_incr[20] = 0.8660254037844387*coeff[0]*fSkin[12]-1.0*coeff[0]*fSkin[20]; 
  boundSurf_incr[24] = 0.8660254037844387*coeff[0]*fSkin[13]-1.0*coeff[0]*fSkin[24]; 
  boundSurf_incr[27] = 0.8660254037844386*coeff[0]*fSkin[17]-1.0*coeff[0]*fSkin[27]; 
  boundSurf_incr[28] = 0.8660254037844387*coeff[0]*fSkin[19]-1.0*coeff[0]*fSkin[28]; 
  boundSurf_incr[29] = 0.8660254037844387*coeff[0]*fSkin[21]-1.0*coeff[0]*fSkin[29]; 
  boundSurf_incr[30] = 0.8660254037844387*coeff[0]*fSkin[22]-1.0*coeff[0]*fSkin[30]; 
  boundSurf_incr[32] = 0.8660254037844387*coeff[0]*fSkin[23]-1.0*coeff[0]*fSkin[32]; 
  boundSurf_incr[34] = 0.8660254037844387*coeff[0]*fSkin[25]-1.0*coeff[0]*fSkin[34]; 
  boundSurf_incr[35] = 0.8660254037844387*coeff[0]*fSkin[26]-1.0*coeff[0]*fSkin[35]; 
  boundSurf_incr[40] = 0.8660254037844387*coeff[0]*fSkin[31]-1.0*coeff[0]*fSkin[40]; 
  boundSurf_incr[41] = 0.8660254037844387*coeff[0]*fSkin[33]-1.0*coeff[0]*fSkin[41]; 
  boundSurf_incr[42] = 0.8660254037844387*coeff[0]*fSkin[36]-1.0*coeff[0]*fSkin[42]; 
  boundSurf_incr[43] = 0.8660254037844386*coeff[0]*fSkin[37]-1.0*coeff[0]*fSkin[43]; 
  boundSurf_incr[45] = 0.8660254037844386*coeff[0]*fSkin[38]-1.0*coeff[0]*fSkin[45]; 
  boundSurf_incr[47] = 0.8660254037844386*coeff[0]*fSkin[39]-1.0*coeff[0]*fSkin[47]; 
  boundSurf_incr[49] = 0.8660254037844386*coeff[0]*fSkin[44]-1.0*coeff[0]*fSkin[49]; 
  boundSurf_incr[50] = 0.8660254037844386*coeff[0]*fSkin[46]-1.0*coeff[0]*fSkin[50]; 
  boundSurf_incr[51] = 0.8660254037844386*coeff[0]*fSkin[48]-1.0*coeff[0]*fSkin[51]; 
  boundSurf_incr[53] = 0.8660254037844387*coeff[0]*fSkin[52]-1.0*coeff[0]*fSkin[53]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[0]*fSkin[1]+0.5412658773652741*coeff[0]*fEdge[1]-0.5625*coeff[0]*fSkin[0]+0.5625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.4375*coeff[0]*fSkin[1])-0.4375*coeff[0]*fEdge[1]+1.4072912811497125*coeff[0]*fSkin[0]-0.5412658773652739*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.5412658773652741*coeff[0]*fSkin[5]+0.5412658773652741*coeff[0]*fEdge[5]-0.5625*coeff[0]*fSkin[2]+0.5625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.5412658773652741*coeff[0]*fSkin[6]+0.5412658773652741*coeff[0]*fEdge[6]-0.5625*coeff[0]*fSkin[3]+0.5625*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 0.5412658773652741*coeff[0]*fSkin[8]+0.5412658773652741*coeff[0]*fEdge[8]-0.5625*coeff[0]*fSkin[4]+0.5625*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = -(1.4375*coeff[0]*fSkin[5])-0.4375*coeff[0]*fEdge[5]+1.4072912811497125*coeff[0]*fSkin[2]-0.5412658773652739*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = -(1.4375*coeff[0]*fSkin[6])-0.4375*coeff[0]*fEdge[6]+1.4072912811497125*coeff[0]*fSkin[3]-0.5412658773652739*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 0.5412658773652741*coeff[0]*fSkin[14]+0.5412658773652741*coeff[0]*fEdge[14]-0.5625*coeff[0]*fSkin[7]+0.5625*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = -(1.4375*coeff[0]*fSkin[8])-0.4375*coeff[0]*fEdge[8]+1.4072912811497125*coeff[0]*fSkin[4]-0.5412658773652739*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 0.5412658773652741*coeff[0]*fSkin[15]+0.5412658773652741*coeff[0]*fEdge[15]-0.5625*coeff[0]*fSkin[9]+0.5625*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 0.5412658773652741*coeff[0]*fSkin[16]+0.5412658773652741*coeff[0]*fEdge[16]-0.5625*coeff[0]*fSkin[10]+0.5625*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 0.5412658773652742*coeff[0]*fSkin[18]+0.5412658773652742*coeff[0]*fEdge[18]-0.5625*coeff[0]*fSkin[11]+0.5625*coeff[0]*fEdge[11]; 
  edgeSurf_incr[12] = 0.5412658773652742*coeff[0]*fSkin[20]+0.5412658773652742*coeff[0]*fEdge[20]-0.5625*coeff[0]*fSkin[12]+0.5625*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = 0.5412658773652742*coeff[0]*fSkin[24]+0.5412658773652742*coeff[0]*fEdge[24]-0.5625*coeff[0]*fSkin[13]+0.5625*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = -(1.4375*coeff[0]*fSkin[14])-0.4375*coeff[0]*fEdge[14]+1.4072912811497125*coeff[0]*fSkin[7]-0.5412658773652739*coeff[0]*fEdge[7]; 
  edgeSurf_incr[15] = -(1.4375*coeff[0]*fSkin[15])-0.4375*coeff[0]*fEdge[15]+1.4072912811497125*coeff[0]*fSkin[9]-0.5412658773652739*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = -(1.4375*coeff[0]*fSkin[16])-0.4375*coeff[0]*fEdge[16]+1.4072912811497125*coeff[0]*fSkin[10]-0.5412658773652739*coeff[0]*fEdge[10]; 
  edgeSurf_incr[17] = 0.5412658773652741*coeff[0]*fSkin[27]+0.5412658773652741*coeff[0]*fEdge[27]-0.5625*coeff[0]*fSkin[17]+0.5625*coeff[0]*fEdge[17]; 
  edgeSurf_incr[18] = -(1.4375*coeff[0]*fSkin[18])-0.4375*coeff[0]*fEdge[18]+1.4072912811497127*coeff[0]*fSkin[11]-0.5412658773652742*coeff[0]*fEdge[11]; 
  edgeSurf_incr[19] = 0.5412658773652742*coeff[0]*fSkin[28]+0.5412658773652742*coeff[0]*fEdge[28]-0.5625*coeff[0]*fSkin[19]+0.5625*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = -(1.4375*coeff[0]*fSkin[20])-0.4375*coeff[0]*fEdge[20]+1.4072912811497127*coeff[0]*fSkin[12]-0.5412658773652742*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = 0.5412658773652742*coeff[0]*fSkin[29]+0.5412658773652742*coeff[0]*fEdge[29]-0.5625*coeff[0]*fSkin[21]+0.5625*coeff[0]*fEdge[21]; 
  edgeSurf_incr[22] = 0.5412658773652742*coeff[0]*fSkin[30]+0.5412658773652742*coeff[0]*fEdge[30]-0.5625*coeff[0]*fSkin[22]+0.5625*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 0.5412658773652742*coeff[0]*fSkin[32]+0.5412658773652742*coeff[0]*fEdge[32]-0.5625*coeff[0]*fSkin[23]+0.5625*coeff[0]*fEdge[23]; 
  edgeSurf_incr[24] = -(1.4375*coeff[0]*fSkin[24])-0.4375*coeff[0]*fEdge[24]+1.4072912811497127*coeff[0]*fSkin[13]-0.5412658773652742*coeff[0]*fEdge[13]; 
  edgeSurf_incr[25] = 0.5412658773652742*coeff[0]*fSkin[34]+0.5412658773652742*coeff[0]*fEdge[34]-0.5625*coeff[0]*fSkin[25]+0.5625*coeff[0]*fEdge[25]; 
  edgeSurf_incr[26] = 0.5412658773652742*coeff[0]*fSkin[35]+0.5412658773652742*coeff[0]*fEdge[35]-0.5625*coeff[0]*fSkin[26]+0.5625*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = -(1.4375*coeff[0]*fSkin[27])-0.4375*coeff[0]*fEdge[27]+1.4072912811497125*coeff[0]*fSkin[17]-0.5412658773652739*coeff[0]*fEdge[17]; 
  edgeSurf_incr[28] = -(1.4375*coeff[0]*fSkin[28])-0.4375*coeff[0]*fEdge[28]+1.4072912811497127*coeff[0]*fSkin[19]-0.5412658773652742*coeff[0]*fEdge[19]; 
  edgeSurf_incr[29] = -(1.4375*coeff[0]*fSkin[29])-0.4375*coeff[0]*fEdge[29]+1.4072912811497127*coeff[0]*fSkin[21]-0.5412658773652742*coeff[0]*fEdge[21]; 
  edgeSurf_incr[30] = -(1.4375*coeff[0]*fSkin[30])-0.4375*coeff[0]*fEdge[30]+1.4072912811497127*coeff[0]*fSkin[22]-0.5412658773652742*coeff[0]*fEdge[22]; 
  edgeSurf_incr[31] = 0.5412658773652742*coeff[0]*fSkin[40]+0.5412658773652742*coeff[0]*fEdge[40]-0.5625*coeff[0]*fSkin[31]+0.5625*coeff[0]*fEdge[31]; 
  edgeSurf_incr[32] = -(1.4375*coeff[0]*fSkin[32])-0.4375*coeff[0]*fEdge[32]+1.4072912811497127*coeff[0]*fSkin[23]-0.5412658773652742*coeff[0]*fEdge[23]; 
  edgeSurf_incr[33] = 0.5412658773652742*coeff[0]*fSkin[41]+0.5412658773652742*coeff[0]*fEdge[41]-0.5625*coeff[0]*fSkin[33]+0.5625*coeff[0]*fEdge[33]; 
  edgeSurf_incr[34] = -(1.4375*coeff[0]*fSkin[34])-0.4375*coeff[0]*fEdge[34]+1.4072912811497127*coeff[0]*fSkin[25]-0.5412658773652742*coeff[0]*fEdge[25]; 
  edgeSurf_incr[35] = -(1.4375*coeff[0]*fSkin[35])-0.4375*coeff[0]*fEdge[35]+1.4072912811497127*coeff[0]*fSkin[26]-0.5412658773652742*coeff[0]*fEdge[26]; 
  edgeSurf_incr[36] = 0.5412658773652742*coeff[0]*fSkin[42]+0.5412658773652742*coeff[0]*fEdge[42]-0.5625*coeff[0]*fSkin[36]+0.5625*coeff[0]*fEdge[36]; 
  edgeSurf_incr[37] = 0.5412658773652741*coeff[0]*fSkin[43]+0.5412658773652741*coeff[0]*fEdge[43]-0.5625*coeff[0]*fSkin[37]+0.5625*coeff[0]*fEdge[37]; 
  edgeSurf_incr[38] = 0.5412658773652741*coeff[0]*fSkin[45]+0.5412658773652741*coeff[0]*fEdge[45]-0.5625*coeff[0]*fSkin[38]+0.5625*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = 0.5412658773652741*coeff[0]*fSkin[47]+0.5412658773652741*coeff[0]*fEdge[47]-0.5625*coeff[0]*fSkin[39]+0.5625*coeff[0]*fEdge[39]; 
  edgeSurf_incr[40] = -(1.4375*coeff[0]*fSkin[40])-0.4375*coeff[0]*fEdge[40]+1.4072912811497127*coeff[0]*fSkin[31]-0.5412658773652742*coeff[0]*fEdge[31]; 
  edgeSurf_incr[41] = -(1.4375*coeff[0]*fSkin[41])-0.4375*coeff[0]*fEdge[41]+1.4072912811497127*coeff[0]*fSkin[33]-0.5412658773652742*coeff[0]*fEdge[33]; 
  edgeSurf_incr[42] = -(1.4375*coeff[0]*fSkin[42])-0.4375*coeff[0]*fEdge[42]+1.4072912811497127*coeff[0]*fSkin[36]-0.5412658773652742*coeff[0]*fEdge[36]; 
  edgeSurf_incr[43] = -(1.4375*coeff[0]*fSkin[43])-0.4375*coeff[0]*fEdge[43]+1.4072912811497125*coeff[0]*fSkin[37]-0.5412658773652739*coeff[0]*fEdge[37]; 
  edgeSurf_incr[44] = 0.5412658773652741*coeff[0]*fSkin[49]+0.5412658773652741*coeff[0]*fEdge[49]-0.5625*coeff[0]*fSkin[44]+0.5625*coeff[0]*fEdge[44]; 
  edgeSurf_incr[45] = -(1.4375*coeff[0]*fSkin[45])-0.4375*coeff[0]*fEdge[45]+1.4072912811497125*coeff[0]*fSkin[38]-0.5412658773652739*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = 0.5412658773652741*coeff[0]*fSkin[50]+0.5412658773652741*coeff[0]*fEdge[50]-0.5625*coeff[0]*fSkin[46]+0.5625*coeff[0]*fEdge[46]; 
  edgeSurf_incr[47] = -(1.4375*coeff[0]*fSkin[47])-0.4375*coeff[0]*fEdge[47]+1.4072912811497125*coeff[0]*fSkin[39]-0.5412658773652739*coeff[0]*fEdge[39]; 
  edgeSurf_incr[48] = 0.5412658773652741*coeff[0]*fSkin[51]+0.5412658773652741*coeff[0]*fEdge[51]-0.5625*coeff[0]*fSkin[48]+0.5625*coeff[0]*fEdge[48]; 
  edgeSurf_incr[49] = -(1.4375*coeff[0]*fSkin[49])-0.4375*coeff[0]*fEdge[49]+1.4072912811497125*coeff[0]*fSkin[44]-0.5412658773652739*coeff[0]*fEdge[44]; 
  edgeSurf_incr[50] = -(1.4375*coeff[0]*fSkin[50])-0.4375*coeff[0]*fEdge[50]+1.4072912811497125*coeff[0]*fSkin[46]-0.5412658773652739*coeff[0]*fEdge[46]; 
  edgeSurf_incr[51] = -(1.4375*coeff[0]*fSkin[51])-0.4375*coeff[0]*fEdge[51]+1.4072912811497125*coeff[0]*fSkin[48]-0.5412658773652739*coeff[0]*fEdge[48]; 
  edgeSurf_incr[52] = 0.5412658773652742*coeff[0]*fSkin[53]+0.5412658773652742*coeff[0]*fEdge[53]-0.5625*coeff[0]*fSkin[52]+0.5625*coeff[0]*fEdge[52]; 
  edgeSurf_incr[53] = -(1.4375*coeff[0]*fSkin[53])-0.4375*coeff[0]*fEdge[53]+1.4072912811497127*coeff[0]*fSkin[52]-0.5412658773652742*coeff[0]*fEdge[52]; 

  boundSurf_incr[1] = -(1.0*coeff[0]*fSkin[1])-0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[5] = -(1.0*coeff[0]*fSkin[5])-0.8660254037844386*coeff[0]*fSkin[2]; 
  boundSurf_incr[6] = -(1.0*coeff[0]*fSkin[6])-0.8660254037844386*coeff[0]*fSkin[3]; 
  boundSurf_incr[8] = -(1.0*coeff[0]*fSkin[8])-0.8660254037844386*coeff[0]*fSkin[4]; 
  boundSurf_incr[14] = -(1.0*coeff[0]*fSkin[14])-0.8660254037844386*coeff[0]*fSkin[7]; 
  boundSurf_incr[15] = -(1.0*coeff[0]*fSkin[15])-0.8660254037844386*coeff[0]*fSkin[9]; 
  boundSurf_incr[16] = -(1.0*coeff[0]*fSkin[16])-0.8660254037844386*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = -(1.0*coeff[0]*fSkin[18])-0.8660254037844387*coeff[0]*fSkin[11]; 
  boundSurf_incr[20] = -(1.0*coeff[0]*fSkin[20])-0.8660254037844387*coeff[0]*fSkin[12]; 
  boundSurf_incr[24] = -(1.0*coeff[0]*fSkin[24])-0.8660254037844387*coeff[0]*fSkin[13]; 
  boundSurf_incr[27] = -(1.0*coeff[0]*fSkin[27])-0.8660254037844386*coeff[0]*fSkin[17]; 
  boundSurf_incr[28] = -(1.0*coeff[0]*fSkin[28])-0.8660254037844387*coeff[0]*fSkin[19]; 
  boundSurf_incr[29] = -(1.0*coeff[0]*fSkin[29])-0.8660254037844387*coeff[0]*fSkin[21]; 
  boundSurf_incr[30] = -(1.0*coeff[0]*fSkin[30])-0.8660254037844387*coeff[0]*fSkin[22]; 
  boundSurf_incr[32] = -(1.0*coeff[0]*fSkin[32])-0.8660254037844387*coeff[0]*fSkin[23]; 
  boundSurf_incr[34] = -(1.0*coeff[0]*fSkin[34])-0.8660254037844387*coeff[0]*fSkin[25]; 
  boundSurf_incr[35] = -(1.0*coeff[0]*fSkin[35])-0.8660254037844387*coeff[0]*fSkin[26]; 
  boundSurf_incr[40] = -(1.0*coeff[0]*fSkin[40])-0.8660254037844387*coeff[0]*fSkin[31]; 
  boundSurf_incr[41] = -(1.0*coeff[0]*fSkin[41])-0.8660254037844387*coeff[0]*fSkin[33]; 
  boundSurf_incr[42] = -(1.0*coeff[0]*fSkin[42])-0.8660254037844387*coeff[0]*fSkin[36]; 
  boundSurf_incr[43] = -(1.0*coeff[0]*fSkin[43])-0.8660254037844386*coeff[0]*fSkin[37]; 
  boundSurf_incr[45] = -(1.0*coeff[0]*fSkin[45])-0.8660254037844386*coeff[0]*fSkin[38]; 
  boundSurf_incr[47] = -(1.0*coeff[0]*fSkin[47])-0.8660254037844386*coeff[0]*fSkin[39]; 
  boundSurf_incr[49] = -(1.0*coeff[0]*fSkin[49])-0.8660254037844386*coeff[0]*fSkin[44]; 
  boundSurf_incr[50] = -(1.0*coeff[0]*fSkin[50])-0.8660254037844386*coeff[0]*fSkin[46]; 
  boundSurf_incr[51] = -(1.0*coeff[0]*fSkin[51])-0.8660254037844386*coeff[0]*fSkin[48]; 
  boundSurf_incr[53] = -(1.0*coeff[0]*fSkin[53])-0.8660254037844387*coeff[0]*fSkin[52]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac; 
  out[24] += (vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac; 
  out[25] += (vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac; 
  out[26] += (vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac; 
  out[27] += (vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*Jfac; 
  out[28] += (vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*Jfac; 
  out[29] += (vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*Jfac; 
  out[30] += (vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*Jfac; 
  out[31] += (vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*Jfac; 
  out[32] += (vol_incr[32]+edgeSurf_incr[32]+boundSurf_incr[32])*Jfac; 
  out[33] += (vol_incr[33]+edgeSurf_incr[33]+boundSurf_incr[33])*Jfac; 
  out[34] += (vol_incr[34]+edgeSurf_incr[34]+boundSurf_incr[34])*Jfac; 
  out[35] += (vol_incr[35]+edgeSurf_incr[35]+boundSurf_incr[35])*Jfac; 
  out[36] += (vol_incr[36]+edgeSurf_incr[36]+boundSurf_incr[36])*Jfac; 
  out[37] += (vol_incr[37]+edgeSurf_incr[37]+boundSurf_incr[37])*Jfac; 
  out[38] += (vol_incr[38]+edgeSurf_incr[38]+boundSurf_incr[38])*Jfac; 
  out[39] += (vol_incr[39]+edgeSurf_incr[39]+boundSurf_incr[39])*Jfac; 
  out[40] += (vol_incr[40]+edgeSurf_incr[40]+boundSurf_incr[40])*Jfac; 
  out[41] += (vol_incr[41]+edgeSurf_incr[41]+boundSurf_incr[41])*Jfac; 
  out[42] += (vol_incr[42]+edgeSurf_incr[42]+boundSurf_incr[42])*Jfac; 
  out[43] += (vol_incr[43]+edgeSurf_incr[43]+boundSurf_incr[43])*Jfac; 
  out[44] += (vol_incr[44]+edgeSurf_incr[44]+boundSurf_incr[44])*Jfac; 
  out[45] += (vol_incr[45]+edgeSurf_incr[45]+boundSurf_incr[45])*Jfac; 
  out[46] += (vol_incr[46]+edgeSurf_incr[46]+boundSurf_incr[46])*Jfac; 
  out[47] += (vol_incr[47]+edgeSurf_incr[47]+boundSurf_incr[47])*Jfac; 
  out[48] += (vol_incr[48]+edgeSurf_incr[48]+boundSurf_incr[48])*Jfac; 
  out[49] += (vol_incr[49]+edgeSurf_incr[49]+boundSurf_incr[49])*Jfac; 
  out[50] += (vol_incr[50]+edgeSurf_incr[50]+boundSurf_incr[50])*Jfac; 
  out[51] += (vol_incr[51]+edgeSurf_incr[51]+boundSurf_incr[51])*Jfac; 
  out[52] += (vol_incr[52]+edgeSurf_incr[52]+boundSurf_incr[52])*Jfac; 
  out[53] += (vol_incr[53]+edgeSurf_incr[53]+boundSurf_incr[53])*Jfac; 

  return 0.;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_boundary_surfx_1x3v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[54] = {0.0}; 
  vol_incr[1] = 2.1213203435596424*fSkin[0]*coeff[1]; 
  vol_incr[5] = 2.1213203435596424*coeff[1]*fSkin[2]; 
  vol_incr[6] = 2.1213203435596424*coeff[1]*fSkin[3]; 
  vol_incr[8] = 2.1213203435596424*coeff[1]*fSkin[4]; 
  vol_incr[14] = 2.1213203435596424*coeff[1]*fSkin[7]; 
  vol_incr[15] = 2.1213203435596424*coeff[1]*fSkin[9]; 
  vol_incr[16] = 2.1213203435596424*coeff[1]*fSkin[10]; 
  vol_incr[18] = 2.1213203435596424*coeff[1]*fSkin[11]; 
  vol_incr[20] = 2.1213203435596424*coeff[1]*fSkin[12]; 
  vol_incr[24] = 2.1213203435596424*coeff[1]*fSkin[13]; 
  vol_incr[27] = 2.1213203435596424*coeff[1]*fSkin[17]; 
  vol_incr[28] = 2.1213203435596424*coeff[1]*fSkin[19]; 
  vol_incr[29] = 2.1213203435596424*coeff[1]*fSkin[21]; 
  vol_incr[30] = 2.1213203435596424*coeff[1]*fSkin[22]; 
  vol_incr[32] = 2.1213203435596424*coeff[1]*fSkin[23]; 
  vol_incr[34] = 2.1213203435596424*coeff[1]*fSkin[25]; 
  vol_incr[35] = 2.1213203435596424*coeff[1]*fSkin[26]; 
  vol_incr[40] = 2.1213203435596424*coeff[1]*fSkin[31]; 
  vol_incr[41] = 2.1213203435596424*coeff[1]*fSkin[33]; 
  vol_incr[42] = 2.1213203435596424*coeff[1]*fSkin[36]; 
  vol_incr[43] = 2.1213203435596424*coeff[1]*fSkin[37]; 
  vol_incr[45] = 2.1213203435596424*coeff[1]*fSkin[38]; 
  vol_incr[47] = 2.1213203435596424*coeff[1]*fSkin[39]; 
  vol_incr[49] = 2.1213203435596424*coeff[1]*fSkin[44]; 
  vol_incr[50] = 2.1213203435596424*coeff[1]*fSkin[46]; 
  vol_incr[51] = 2.1213203435596424*coeff[1]*fSkin[48]; 
  vol_incr[53] = 2.1213203435596424*coeff[1]*fSkin[52]; 

  double edgeSurf_incr[54] = {0.0}; 
  double boundSurf_incr[54] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.38273277230987135*coeff[0]*fSkin[1])-0.38273277230987135*coeff[0]*fEdge[1]-0.39774756441743275*coeff[0]*fSkin[0]+0.39774756441743275*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(0.6123724356957944*coeff[1]*fSkin[1])-1.0164659979556616*coeff[0]*fSkin[1]+0.6123724356957944*coeff[1]*fEdge[1]-0.3093592167691142*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]-0.9951052080056654*coeff[0]*fSkin[0]+0.3827327723098711*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(0.38273277230987135*coeff[0]*fSkin[5])-0.38273277230987135*coeff[0]*fEdge[5]-0.39774756441743275*coeff[0]*fSkin[2]+0.39774756441743275*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(0.38273277230987135*coeff[0]*fSkin[6])-0.38273277230987135*coeff[0]*fEdge[6]-0.39774756441743275*coeff[0]*fSkin[3]+0.39774756441743275*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = -(0.38273277230987135*coeff[0]*fSkin[8])-0.38273277230987135*coeff[0]*fEdge[8]-0.39774756441743275*coeff[0]*fSkin[4]+0.39774756441743275*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = -(0.6123724356957944*coeff[1]*fSkin[5])-1.0164659979556616*coeff[0]*fSkin[5]+0.6123724356957944*coeff[1]*fEdge[5]-0.3093592167691142*coeff[0]*fEdge[5]-0.5303300858899105*coeff[1]*fSkin[2]-0.9951052080056654*coeff[0]*fSkin[2]-0.5303300858899105*coeff[1]*fEdge[2]+0.3827327723098711*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = -(0.6123724356957944*coeff[1]*fSkin[6])-1.0164659979556616*coeff[0]*fSkin[6]+0.6123724356957944*coeff[1]*fEdge[6]-0.3093592167691142*coeff[0]*fEdge[6]-0.5303300858899105*coeff[1]*fSkin[3]-0.9951052080056654*coeff[0]*fSkin[3]-0.5303300858899105*coeff[1]*fEdge[3]+0.3827327723098711*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = -(0.38273277230987135*coeff[0]*fSkin[14])-0.38273277230987135*coeff[0]*fEdge[14]-0.39774756441743275*coeff[0]*fSkin[7]+0.39774756441743275*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = -(0.6123724356957944*coeff[1]*fSkin[8])-1.0164659979556616*coeff[0]*fSkin[8]+0.6123724356957944*coeff[1]*fEdge[8]-0.3093592167691142*coeff[0]*fEdge[8]-0.5303300858899105*coeff[1]*fSkin[4]-0.9951052080056654*coeff[0]*fSkin[4]-0.5303300858899105*coeff[1]*fEdge[4]+0.3827327723098711*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = -(0.38273277230987135*coeff[0]*fSkin[15])-0.38273277230987135*coeff[0]*fEdge[15]-0.39774756441743275*coeff[0]*fSkin[9]+0.39774756441743275*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = -(0.38273277230987135*coeff[0]*fSkin[16])-0.38273277230987135*coeff[0]*fEdge[16]-0.39774756441743275*coeff[0]*fSkin[10]+0.39774756441743275*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = -(0.3827327723098714*coeff[0]*fSkin[18])-0.3827327723098714*coeff[0]*fEdge[18]-0.39774756441743275*coeff[0]*fSkin[11]+0.39774756441743275*coeff[0]*fEdge[11]; 
  edgeSurf_incr[12] = -(0.3827327723098714*coeff[0]*fSkin[20])-0.3827327723098714*coeff[0]*fEdge[20]-0.39774756441743275*coeff[0]*fSkin[12]+0.39774756441743275*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = -(0.3827327723098714*coeff[0]*fSkin[24])-0.3827327723098714*coeff[0]*fEdge[24]-0.39774756441743275*coeff[0]*fSkin[13]+0.39774756441743275*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = -(0.6123724356957944*coeff[1]*fSkin[14])-1.0164659979556616*coeff[0]*fSkin[14]+0.6123724356957944*coeff[1]*fEdge[14]-0.3093592167691142*coeff[0]*fEdge[14]-0.5303300858899105*coeff[1]*fSkin[7]-0.9951052080056654*coeff[0]*fSkin[7]-0.5303300858899105*coeff[1]*fEdge[7]+0.3827327723098711*coeff[0]*fEdge[7]; 
  edgeSurf_incr[15] = -(0.6123724356957944*coeff[1]*fSkin[15])-1.0164659979556616*coeff[0]*fSkin[15]+0.6123724356957944*coeff[1]*fEdge[15]-0.3093592167691142*coeff[0]*fEdge[15]-0.5303300858899105*coeff[1]*fSkin[9]-0.9951052080056654*coeff[0]*fSkin[9]-0.5303300858899105*coeff[1]*fEdge[9]+0.3827327723098711*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = -(0.6123724356957944*coeff[1]*fSkin[16])-1.0164659979556616*coeff[0]*fSkin[16]+0.6123724356957944*coeff[1]*fEdge[16]-0.3093592167691142*coeff[0]*fEdge[16]-0.5303300858899105*coeff[1]*fSkin[10]-0.9951052080056654*coeff[0]*fSkin[10]-0.5303300858899105*coeff[1]*fEdge[10]+0.3827327723098711*coeff[0]*fEdge[10]; 
  edgeSurf_incr[17] = -(0.38273277230987135*coeff[0]*fSkin[27])-0.38273277230987135*coeff[0]*fEdge[27]-0.39774756441743275*coeff[0]*fSkin[17]+0.39774756441743275*coeff[0]*fEdge[17]; 
  edgeSurf_incr[18] = -(0.6123724356957944*coeff[1]*fSkin[18])-1.0164659979556616*coeff[0]*fSkin[18]+0.6123724356957944*coeff[1]*fEdge[18]-0.3093592167691142*coeff[0]*fEdge[18]-0.5303300858899104*coeff[1]*fSkin[11]-0.9951052080056656*coeff[0]*fSkin[11]-0.5303300858899104*coeff[1]*fEdge[11]+0.38273277230987135*coeff[0]*fEdge[11]; 
  edgeSurf_incr[19] = -(0.3827327723098714*coeff[0]*fSkin[28])-0.3827327723098714*coeff[0]*fEdge[28]-0.39774756441743275*coeff[0]*fSkin[19]+0.39774756441743275*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = -(0.6123724356957944*coeff[1]*fSkin[20])-1.0164659979556616*coeff[0]*fSkin[20]+0.6123724356957944*coeff[1]*fEdge[20]-0.3093592167691142*coeff[0]*fEdge[20]-0.5303300858899104*coeff[1]*fSkin[12]-0.9951052080056656*coeff[0]*fSkin[12]-0.5303300858899104*coeff[1]*fEdge[12]+0.38273277230987135*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = -(0.3827327723098714*coeff[0]*fSkin[29])-0.3827327723098714*coeff[0]*fEdge[29]-0.39774756441743275*coeff[0]*fSkin[21]+0.39774756441743275*coeff[0]*fEdge[21]; 
  edgeSurf_incr[22] = -(0.3827327723098714*coeff[0]*fSkin[30])-0.3827327723098714*coeff[0]*fEdge[30]-0.39774756441743275*coeff[0]*fSkin[22]+0.39774756441743275*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = -(0.3827327723098714*coeff[0]*fSkin[32])-0.3827327723098714*coeff[0]*fEdge[32]-0.39774756441743275*coeff[0]*fSkin[23]+0.39774756441743275*coeff[0]*fEdge[23]; 
  edgeSurf_incr[24] = -(0.6123724356957944*coeff[1]*fSkin[24])-1.0164659979556616*coeff[0]*fSkin[24]+0.6123724356957944*coeff[1]*fEdge[24]-0.3093592167691142*coeff[0]*fEdge[24]-0.5303300858899104*coeff[1]*fSkin[13]-0.9951052080056656*coeff[0]*fSkin[13]-0.5303300858899104*coeff[1]*fEdge[13]+0.38273277230987135*coeff[0]*fEdge[13]; 
  edgeSurf_incr[25] = -(0.3827327723098714*coeff[0]*fSkin[34])-0.3827327723098714*coeff[0]*fEdge[34]-0.39774756441743275*coeff[0]*fSkin[25]+0.39774756441743275*coeff[0]*fEdge[25]; 
  edgeSurf_incr[26] = -(0.3827327723098714*coeff[0]*fSkin[35])-0.3827327723098714*coeff[0]*fEdge[35]-0.39774756441743275*coeff[0]*fSkin[26]+0.39774756441743275*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = -(0.6123724356957944*coeff[1]*fSkin[27])-1.0164659979556616*coeff[0]*fSkin[27]+0.6123724356957944*coeff[1]*fEdge[27]-0.3093592167691142*coeff[0]*fEdge[27]-0.5303300858899105*coeff[1]*fSkin[17]-0.9951052080056654*coeff[0]*fSkin[17]-0.5303300858899105*coeff[1]*fEdge[17]+0.3827327723098711*coeff[0]*fEdge[17]; 
  edgeSurf_incr[28] = -(0.6123724356957944*coeff[1]*fSkin[28])-1.0164659979556616*coeff[0]*fSkin[28]+0.6123724356957944*coeff[1]*fEdge[28]-0.3093592167691142*coeff[0]*fEdge[28]-0.5303300858899104*coeff[1]*fSkin[19]-0.9951052080056656*coeff[0]*fSkin[19]-0.5303300858899104*coeff[1]*fEdge[19]+0.38273277230987135*coeff[0]*fEdge[19]; 
  edgeSurf_incr[29] = -(0.6123724356957944*coeff[1]*fSkin[29])-1.0164659979556616*coeff[0]*fSkin[29]+0.6123724356957944*coeff[1]*fEdge[29]-0.3093592167691142*coeff[0]*fEdge[29]-0.5303300858899104*coeff[1]*fSkin[21]-0.9951052080056656*coeff[0]*fSkin[21]-0.5303300858899104*coeff[1]*fEdge[21]+0.38273277230987135*coeff[0]*fEdge[21]; 
  edgeSurf_incr[30] = -(0.6123724356957944*coeff[1]*fSkin[30])-1.0164659979556616*coeff[0]*fSkin[30]+0.6123724356957944*coeff[1]*fEdge[30]-0.3093592167691142*coeff[0]*fEdge[30]-0.5303300858899104*coeff[1]*fSkin[22]-0.9951052080056656*coeff[0]*fSkin[22]-0.5303300858899104*coeff[1]*fEdge[22]+0.38273277230987135*coeff[0]*fEdge[22]; 
  edgeSurf_incr[31] = -(0.3827327723098714*coeff[0]*fSkin[40])-0.3827327723098714*coeff[0]*fEdge[40]-0.39774756441743275*coeff[0]*fSkin[31]+0.39774756441743275*coeff[0]*fEdge[31]; 
  edgeSurf_incr[32] = -(0.6123724356957944*coeff[1]*fSkin[32])-1.0164659979556616*coeff[0]*fSkin[32]+0.6123724356957944*coeff[1]*fEdge[32]-0.3093592167691142*coeff[0]*fEdge[32]-0.5303300858899104*coeff[1]*fSkin[23]-0.9951052080056656*coeff[0]*fSkin[23]-0.5303300858899104*coeff[1]*fEdge[23]+0.38273277230987135*coeff[0]*fEdge[23]; 
  edgeSurf_incr[33] = -(0.3827327723098714*coeff[0]*fSkin[41])-0.3827327723098714*coeff[0]*fEdge[41]-0.39774756441743275*coeff[0]*fSkin[33]+0.39774756441743275*coeff[0]*fEdge[33]; 
  edgeSurf_incr[34] = -(0.6123724356957944*coeff[1]*fSkin[34])-1.0164659979556616*coeff[0]*fSkin[34]+0.6123724356957944*coeff[1]*fEdge[34]-0.3093592167691142*coeff[0]*fEdge[34]-0.5303300858899104*coeff[1]*fSkin[25]-0.9951052080056656*coeff[0]*fSkin[25]-0.5303300858899104*coeff[1]*fEdge[25]+0.38273277230987135*coeff[0]*fEdge[25]; 
  edgeSurf_incr[35] = -(0.6123724356957944*coeff[1]*fSkin[35])-1.0164659979556616*coeff[0]*fSkin[35]+0.6123724356957944*coeff[1]*fEdge[35]-0.3093592167691142*coeff[0]*fEdge[35]-0.5303300858899104*coeff[1]*fSkin[26]-0.9951052080056656*coeff[0]*fSkin[26]-0.5303300858899104*coeff[1]*fEdge[26]+0.38273277230987135*coeff[0]*fEdge[26]; 
  edgeSurf_incr[36] = -(0.3827327723098714*coeff[0]*fSkin[42])-0.3827327723098714*coeff[0]*fEdge[42]-0.39774756441743275*coeff[0]*fSkin[36]+0.39774756441743275*coeff[0]*fEdge[36]; 
  edgeSurf_incr[37] = -(0.38273277230987135*coeff[0]*fSkin[43])-0.38273277230987135*coeff[0]*fEdge[43]-0.39774756441743275*coeff[0]*fSkin[37]+0.39774756441743275*coeff[0]*fEdge[37]; 
  edgeSurf_incr[38] = -(0.38273277230987135*coeff[0]*fSkin[45])-0.38273277230987135*coeff[0]*fEdge[45]-0.39774756441743275*coeff[0]*fSkin[38]+0.39774756441743275*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = -(0.38273277230987135*coeff[0]*fSkin[47])-0.38273277230987135*coeff[0]*fEdge[47]-0.39774756441743275*coeff[0]*fSkin[39]+0.39774756441743275*coeff[0]*fEdge[39]; 
  edgeSurf_incr[40] = -(0.6123724356957944*coeff[1]*fSkin[40])-1.0164659979556616*coeff[0]*fSkin[40]+0.6123724356957944*coeff[1]*fEdge[40]-0.3093592167691142*coeff[0]*fEdge[40]-0.5303300858899104*coeff[1]*fSkin[31]-0.9951052080056656*coeff[0]*fSkin[31]-0.5303300858899104*coeff[1]*fEdge[31]+0.38273277230987135*coeff[0]*fEdge[31]; 
  edgeSurf_incr[41] = -(0.6123724356957944*coeff[1]*fSkin[41])-1.0164659979556616*coeff[0]*fSkin[41]+0.6123724356957944*coeff[1]*fEdge[41]-0.3093592167691142*coeff[0]*fEdge[41]-0.5303300858899104*coeff[1]*fSkin[33]-0.9951052080056656*coeff[0]*fSkin[33]-0.5303300858899104*coeff[1]*fEdge[33]+0.38273277230987135*coeff[0]*fEdge[33]; 
  edgeSurf_incr[42] = -(0.6123724356957944*coeff[1]*fSkin[42])-1.0164659979556616*coeff[0]*fSkin[42]+0.6123724356957944*coeff[1]*fEdge[42]-0.3093592167691142*coeff[0]*fEdge[42]-0.5303300858899104*coeff[1]*fSkin[36]-0.9951052080056656*coeff[0]*fSkin[36]-0.5303300858899104*coeff[1]*fEdge[36]+0.38273277230987135*coeff[0]*fEdge[36]; 
  edgeSurf_incr[43] = -(0.6123724356957944*coeff[1]*fSkin[43])-1.0164659979556616*coeff[0]*fSkin[43]+0.6123724356957944*coeff[1]*fEdge[43]-0.3093592167691142*coeff[0]*fEdge[43]-0.5303300858899105*coeff[1]*fSkin[37]-0.9951052080056654*coeff[0]*fSkin[37]-0.5303300858899105*coeff[1]*fEdge[37]+0.3827327723098711*coeff[0]*fEdge[37]; 
  edgeSurf_incr[44] = -(0.38273277230987135*coeff[0]*fSkin[49])-0.38273277230987135*coeff[0]*fEdge[49]-0.39774756441743275*coeff[0]*fSkin[44]+0.39774756441743275*coeff[0]*fEdge[44]; 
  edgeSurf_incr[45] = -(0.6123724356957944*coeff[1]*fSkin[45])-1.0164659979556616*coeff[0]*fSkin[45]+0.6123724356957944*coeff[1]*fEdge[45]-0.3093592167691142*coeff[0]*fEdge[45]-0.5303300858899105*coeff[1]*fSkin[38]-0.9951052080056654*coeff[0]*fSkin[38]-0.5303300858899105*coeff[1]*fEdge[38]+0.3827327723098711*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = -(0.38273277230987135*coeff[0]*fSkin[50])-0.38273277230987135*coeff[0]*fEdge[50]-0.39774756441743275*coeff[0]*fSkin[46]+0.39774756441743275*coeff[0]*fEdge[46]; 
  edgeSurf_incr[47] = -(0.6123724356957944*coeff[1]*fSkin[47])-1.0164659979556616*coeff[0]*fSkin[47]+0.6123724356957944*coeff[1]*fEdge[47]-0.3093592167691142*coeff[0]*fEdge[47]-0.5303300858899105*coeff[1]*fSkin[39]-0.9951052080056654*coeff[0]*fSkin[39]-0.5303300858899105*coeff[1]*fEdge[39]+0.3827327723098711*coeff[0]*fEdge[39]; 
  edgeSurf_incr[48] = -(0.38273277230987135*coeff[0]*fSkin[51])-0.38273277230987135*coeff[0]*fEdge[51]-0.39774756441743275*coeff[0]*fSkin[48]+0.39774756441743275*coeff[0]*fEdge[48]; 
  edgeSurf_incr[49] = -(0.6123724356957944*coeff[1]*fSkin[49])-1.0164659979556616*coeff[0]*fSkin[49]+0.6123724356957944*coeff[1]*fEdge[49]-0.3093592167691142*coeff[0]*fEdge[49]-0.5303300858899105*coeff[1]*fSkin[44]-0.9951052080056654*coeff[0]*fSkin[44]-0.5303300858899105*coeff[1]*fEdge[44]+0.3827327723098711*coeff[0]*fEdge[44]; 
  edgeSurf_incr[50] = -(0.6123724356957944*coeff[1]*fSkin[50])-1.0164659979556616*coeff[0]*fSkin[50]+0.6123724356957944*coeff[1]*fEdge[50]-0.3093592167691142*coeff[0]*fEdge[50]-0.5303300858899105*coeff[1]*fSkin[46]-0.9951052080056654*coeff[0]*fSkin[46]-0.5303300858899105*coeff[1]*fEdge[46]+0.3827327723098711*coeff[0]*fEdge[46]; 
  edgeSurf_incr[51] = -(0.6123724356957944*coeff[1]*fSkin[51])-1.0164659979556616*coeff[0]*fSkin[51]+0.6123724356957944*coeff[1]*fEdge[51]-0.3093592167691142*coeff[0]*fEdge[51]-0.5303300858899105*coeff[1]*fSkin[48]-0.9951052080056654*coeff[0]*fSkin[48]-0.5303300858899105*coeff[1]*fEdge[48]+0.3827327723098711*coeff[0]*fEdge[48]; 
  edgeSurf_incr[52] = -(0.3827327723098714*coeff[0]*fSkin[53])-0.3827327723098714*coeff[0]*fEdge[53]-0.39774756441743275*coeff[0]*fSkin[52]+0.39774756441743275*coeff[0]*fEdge[52]; 
  edgeSurf_incr[53] = -(0.6123724356957944*coeff[1]*fSkin[53])-1.0164659979556616*coeff[0]*fSkin[53]+0.6123724356957944*coeff[1]*fEdge[53]-0.3093592167691142*coeff[0]*fEdge[53]-0.5303300858899104*coeff[1]*fSkin[52]-0.9951052080056656*coeff[0]*fSkin[52]-0.5303300858899104*coeff[1]*fEdge[52]+0.38273277230987135*coeff[0]*fEdge[52]; 

  boundSurf_incr[1] = 1.224744871391589*coeff[1]*fSkin[1]-0.7071067811865475*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]+0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[5] = 1.224744871391589*coeff[1]*fSkin[5]-0.7071067811865475*coeff[0]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[2]+0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[6] = 1.224744871391589*coeff[1]*fSkin[6]-0.7071067811865475*coeff[0]*fSkin[6]-1.060660171779821*coeff[1]*fSkin[3]+0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[8] = 1.224744871391589*coeff[1]*fSkin[8]-0.7071067811865475*coeff[0]*fSkin[8]-1.060660171779821*coeff[1]*fSkin[4]+0.6123724356957944*coeff[0]*fSkin[4]; 
  boundSurf_incr[14] = 1.224744871391589*coeff[1]*fSkin[14]-0.7071067811865475*coeff[0]*fSkin[14]-1.060660171779821*coeff[1]*fSkin[7]+0.6123724356957944*coeff[0]*fSkin[7]; 
  boundSurf_incr[15] = 1.224744871391589*coeff[1]*fSkin[15]-0.7071067811865475*coeff[0]*fSkin[15]-1.060660171779821*coeff[1]*fSkin[9]+0.6123724356957944*coeff[0]*fSkin[9]; 
  boundSurf_incr[16] = 1.224744871391589*coeff[1]*fSkin[16]-0.7071067811865475*coeff[0]*fSkin[16]-1.060660171779821*coeff[1]*fSkin[10]+0.6123724356957944*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = 1.224744871391589*coeff[1]*fSkin[18]-0.7071067811865475*coeff[0]*fSkin[18]-1.060660171779821*coeff[1]*fSkin[11]+0.6123724356957944*coeff[0]*fSkin[11]; 
  boundSurf_incr[20] = 1.224744871391589*coeff[1]*fSkin[20]-0.7071067811865475*coeff[0]*fSkin[20]-1.060660171779821*coeff[1]*fSkin[12]+0.6123724356957944*coeff[0]*fSkin[12]; 
  boundSurf_incr[24] = 1.224744871391589*coeff[1]*fSkin[24]-0.7071067811865475*coeff[0]*fSkin[24]-1.060660171779821*coeff[1]*fSkin[13]+0.6123724356957944*coeff[0]*fSkin[13]; 
  boundSurf_incr[27] = 1.224744871391589*coeff[1]*fSkin[27]-0.7071067811865475*coeff[0]*fSkin[27]-1.060660171779821*coeff[1]*fSkin[17]+0.6123724356957944*coeff[0]*fSkin[17]; 
  boundSurf_incr[28] = 1.224744871391589*coeff[1]*fSkin[28]-0.7071067811865475*coeff[0]*fSkin[28]-1.060660171779821*coeff[1]*fSkin[19]+0.6123724356957944*coeff[0]*fSkin[19]; 
  boundSurf_incr[29] = 1.224744871391589*coeff[1]*fSkin[29]-0.7071067811865475*coeff[0]*fSkin[29]-1.060660171779821*coeff[1]*fSkin[21]+0.6123724356957944*coeff[0]*fSkin[21]; 
  boundSurf_incr[30] = 1.224744871391589*coeff[1]*fSkin[30]-0.7071067811865475*coeff[0]*fSkin[30]-1.060660171779821*coeff[1]*fSkin[22]+0.6123724356957944*coeff[0]*fSkin[22]; 
  boundSurf_incr[32] = 1.224744871391589*coeff[1]*fSkin[32]-0.7071067811865475*coeff[0]*fSkin[32]-1.060660171779821*coeff[1]*fSkin[23]+0.6123724356957944*coeff[0]*fSkin[23]; 
  boundSurf_incr[34] = 1.224744871391589*coeff[1]*fSkin[34]-0.7071067811865475*coeff[0]*fSkin[34]-1.060660171779821*coeff[1]*fSkin[25]+0.6123724356957944*coeff[0]*fSkin[25]; 
  boundSurf_incr[35] = 1.224744871391589*coeff[1]*fSkin[35]-0.7071067811865475*coeff[0]*fSkin[35]-1.060660171779821*coeff[1]*fSkin[26]+0.6123724356957944*coeff[0]*fSkin[26]; 
  boundSurf_incr[40] = 1.224744871391589*coeff[1]*fSkin[40]-0.7071067811865475*coeff[0]*fSkin[40]-1.060660171779821*coeff[1]*fSkin[31]+0.6123724356957944*coeff[0]*fSkin[31]; 
  boundSurf_incr[41] = 1.224744871391589*coeff[1]*fSkin[41]-0.7071067811865475*coeff[0]*fSkin[41]-1.060660171779821*coeff[1]*fSkin[33]+0.6123724356957944*coeff[0]*fSkin[33]; 
  boundSurf_incr[42] = 1.224744871391589*coeff[1]*fSkin[42]-0.7071067811865475*coeff[0]*fSkin[42]-1.060660171779821*coeff[1]*fSkin[36]+0.6123724356957944*coeff[0]*fSkin[36]; 
  boundSurf_incr[43] = 1.224744871391589*coeff[1]*fSkin[43]-0.7071067811865475*coeff[0]*fSkin[43]-1.060660171779821*coeff[1]*fSkin[37]+0.6123724356957944*coeff[0]*fSkin[37]; 
  boundSurf_incr[45] = 1.224744871391589*coeff[1]*fSkin[45]-0.7071067811865475*coeff[0]*fSkin[45]-1.060660171779821*coeff[1]*fSkin[38]+0.6123724356957944*coeff[0]*fSkin[38]; 
  boundSurf_incr[47] = 1.224744871391589*coeff[1]*fSkin[47]-0.7071067811865475*coeff[0]*fSkin[47]-1.060660171779821*coeff[1]*fSkin[39]+0.6123724356957944*coeff[0]*fSkin[39]; 
  boundSurf_incr[49] = 1.224744871391589*coeff[1]*fSkin[49]-0.7071067811865475*coeff[0]*fSkin[49]-1.060660171779821*coeff[1]*fSkin[44]+0.6123724356957944*coeff[0]*fSkin[44]; 
  boundSurf_incr[50] = 1.224744871391589*coeff[1]*fSkin[50]-0.7071067811865475*coeff[0]*fSkin[50]-1.060660171779821*coeff[1]*fSkin[46]+0.6123724356957944*coeff[0]*fSkin[46]; 
  boundSurf_incr[51] = 1.224744871391589*coeff[1]*fSkin[51]-0.7071067811865475*coeff[0]*fSkin[51]-1.060660171779821*coeff[1]*fSkin[48]+0.6123724356957944*coeff[0]*fSkin[48]; 
  boundSurf_incr[53] = 1.224744871391589*coeff[1]*fSkin[53]-0.7071067811865475*coeff[0]*fSkin[53]-1.060660171779821*coeff[1]*fSkin[52]+0.6123724356957944*coeff[0]*fSkin[52]; 

  } else { 

  edgeSurf_incr[0] = 0.38273277230987135*coeff[0]*fSkin[1]+0.38273277230987135*coeff[0]*fEdge[1]-0.39774756441743275*coeff[0]*fSkin[0]+0.39774756441743275*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.6123724356957944*coeff[1]*fSkin[1]-1.0164659979556616*coeff[0]*fSkin[1]-0.6123724356957944*coeff[1]*fEdge[1]-0.3093592167691142*coeff[0]*fEdge[1]-0.5303300858899105*fSkin[0]*coeff[1]-0.5303300858899105*fEdge[0]*coeff[1]+0.9951052080056654*coeff[0]*fSkin[0]-0.3827327723098711*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.38273277230987135*coeff[0]*fSkin[5]+0.38273277230987135*coeff[0]*fEdge[5]-0.39774756441743275*coeff[0]*fSkin[2]+0.39774756441743275*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.38273277230987135*coeff[0]*fSkin[6]+0.38273277230987135*coeff[0]*fEdge[6]-0.39774756441743275*coeff[0]*fSkin[3]+0.39774756441743275*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 0.38273277230987135*coeff[0]*fSkin[8]+0.38273277230987135*coeff[0]*fEdge[8]-0.39774756441743275*coeff[0]*fSkin[4]+0.39774756441743275*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 0.6123724356957944*coeff[1]*fSkin[5]-1.0164659979556616*coeff[0]*fSkin[5]-0.6123724356957944*coeff[1]*fEdge[5]-0.3093592167691142*coeff[0]*fEdge[5]-0.5303300858899105*coeff[1]*fSkin[2]+0.9951052080056654*coeff[0]*fSkin[2]-0.5303300858899105*coeff[1]*fEdge[2]-0.3827327723098711*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = 0.6123724356957944*coeff[1]*fSkin[6]-1.0164659979556616*coeff[0]*fSkin[6]-0.6123724356957944*coeff[1]*fEdge[6]-0.3093592167691142*coeff[0]*fEdge[6]-0.5303300858899105*coeff[1]*fSkin[3]+0.9951052080056654*coeff[0]*fSkin[3]-0.5303300858899105*coeff[1]*fEdge[3]-0.3827327723098711*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 0.38273277230987135*coeff[0]*fSkin[14]+0.38273277230987135*coeff[0]*fEdge[14]-0.39774756441743275*coeff[0]*fSkin[7]+0.39774756441743275*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = 0.6123724356957944*coeff[1]*fSkin[8]-1.0164659979556616*coeff[0]*fSkin[8]-0.6123724356957944*coeff[1]*fEdge[8]-0.3093592167691142*coeff[0]*fEdge[8]-0.5303300858899105*coeff[1]*fSkin[4]+0.9951052080056654*coeff[0]*fSkin[4]-0.5303300858899105*coeff[1]*fEdge[4]-0.3827327723098711*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 0.38273277230987135*coeff[0]*fSkin[15]+0.38273277230987135*coeff[0]*fEdge[15]-0.39774756441743275*coeff[0]*fSkin[9]+0.39774756441743275*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 0.38273277230987135*coeff[0]*fSkin[16]+0.38273277230987135*coeff[0]*fEdge[16]-0.39774756441743275*coeff[0]*fSkin[10]+0.39774756441743275*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 0.3827327723098714*coeff[0]*fSkin[18]+0.3827327723098714*coeff[0]*fEdge[18]-0.39774756441743275*coeff[0]*fSkin[11]+0.39774756441743275*coeff[0]*fEdge[11]; 
  edgeSurf_incr[12] = 0.3827327723098714*coeff[0]*fSkin[20]+0.3827327723098714*coeff[0]*fEdge[20]-0.39774756441743275*coeff[0]*fSkin[12]+0.39774756441743275*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = 0.3827327723098714*coeff[0]*fSkin[24]+0.3827327723098714*coeff[0]*fEdge[24]-0.39774756441743275*coeff[0]*fSkin[13]+0.39774756441743275*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = 0.6123724356957944*coeff[1]*fSkin[14]-1.0164659979556616*coeff[0]*fSkin[14]-0.6123724356957944*coeff[1]*fEdge[14]-0.3093592167691142*coeff[0]*fEdge[14]-0.5303300858899105*coeff[1]*fSkin[7]+0.9951052080056654*coeff[0]*fSkin[7]-0.5303300858899105*coeff[1]*fEdge[7]-0.3827327723098711*coeff[0]*fEdge[7]; 
  edgeSurf_incr[15] = 0.6123724356957944*coeff[1]*fSkin[15]-1.0164659979556616*coeff[0]*fSkin[15]-0.6123724356957944*coeff[1]*fEdge[15]-0.3093592167691142*coeff[0]*fEdge[15]-0.5303300858899105*coeff[1]*fSkin[9]+0.9951052080056654*coeff[0]*fSkin[9]-0.5303300858899105*coeff[1]*fEdge[9]-0.3827327723098711*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = 0.6123724356957944*coeff[1]*fSkin[16]-1.0164659979556616*coeff[0]*fSkin[16]-0.6123724356957944*coeff[1]*fEdge[16]-0.3093592167691142*coeff[0]*fEdge[16]-0.5303300858899105*coeff[1]*fSkin[10]+0.9951052080056654*coeff[0]*fSkin[10]-0.5303300858899105*coeff[1]*fEdge[10]-0.3827327723098711*coeff[0]*fEdge[10]; 
  edgeSurf_incr[17] = 0.38273277230987135*coeff[0]*fSkin[27]+0.38273277230987135*coeff[0]*fEdge[27]-0.39774756441743275*coeff[0]*fSkin[17]+0.39774756441743275*coeff[0]*fEdge[17]; 
  edgeSurf_incr[18] = 0.6123724356957944*coeff[1]*fSkin[18]-1.0164659979556616*coeff[0]*fSkin[18]-0.6123724356957944*coeff[1]*fEdge[18]-0.3093592167691142*coeff[0]*fEdge[18]-0.5303300858899104*coeff[1]*fSkin[11]+0.9951052080056656*coeff[0]*fSkin[11]-0.5303300858899104*coeff[1]*fEdge[11]-0.38273277230987135*coeff[0]*fEdge[11]; 
  edgeSurf_incr[19] = 0.3827327723098714*coeff[0]*fSkin[28]+0.3827327723098714*coeff[0]*fEdge[28]-0.39774756441743275*coeff[0]*fSkin[19]+0.39774756441743275*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = 0.6123724356957944*coeff[1]*fSkin[20]-1.0164659979556616*coeff[0]*fSkin[20]-0.6123724356957944*coeff[1]*fEdge[20]-0.3093592167691142*coeff[0]*fEdge[20]-0.5303300858899104*coeff[1]*fSkin[12]+0.9951052080056656*coeff[0]*fSkin[12]-0.5303300858899104*coeff[1]*fEdge[12]-0.38273277230987135*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = 0.3827327723098714*coeff[0]*fSkin[29]+0.3827327723098714*coeff[0]*fEdge[29]-0.39774756441743275*coeff[0]*fSkin[21]+0.39774756441743275*coeff[0]*fEdge[21]; 
  edgeSurf_incr[22] = 0.3827327723098714*coeff[0]*fSkin[30]+0.3827327723098714*coeff[0]*fEdge[30]-0.39774756441743275*coeff[0]*fSkin[22]+0.39774756441743275*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 0.3827327723098714*coeff[0]*fSkin[32]+0.3827327723098714*coeff[0]*fEdge[32]-0.39774756441743275*coeff[0]*fSkin[23]+0.39774756441743275*coeff[0]*fEdge[23]; 
  edgeSurf_incr[24] = 0.6123724356957944*coeff[1]*fSkin[24]-1.0164659979556616*coeff[0]*fSkin[24]-0.6123724356957944*coeff[1]*fEdge[24]-0.3093592167691142*coeff[0]*fEdge[24]-0.5303300858899104*coeff[1]*fSkin[13]+0.9951052080056656*coeff[0]*fSkin[13]-0.5303300858899104*coeff[1]*fEdge[13]-0.38273277230987135*coeff[0]*fEdge[13]; 
  edgeSurf_incr[25] = 0.3827327723098714*coeff[0]*fSkin[34]+0.3827327723098714*coeff[0]*fEdge[34]-0.39774756441743275*coeff[0]*fSkin[25]+0.39774756441743275*coeff[0]*fEdge[25]; 
  edgeSurf_incr[26] = 0.3827327723098714*coeff[0]*fSkin[35]+0.3827327723098714*coeff[0]*fEdge[35]-0.39774756441743275*coeff[0]*fSkin[26]+0.39774756441743275*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = 0.6123724356957944*coeff[1]*fSkin[27]-1.0164659979556616*coeff[0]*fSkin[27]-0.6123724356957944*coeff[1]*fEdge[27]-0.3093592167691142*coeff[0]*fEdge[27]-0.5303300858899105*coeff[1]*fSkin[17]+0.9951052080056654*coeff[0]*fSkin[17]-0.5303300858899105*coeff[1]*fEdge[17]-0.3827327723098711*coeff[0]*fEdge[17]; 
  edgeSurf_incr[28] = 0.6123724356957944*coeff[1]*fSkin[28]-1.0164659979556616*coeff[0]*fSkin[28]-0.6123724356957944*coeff[1]*fEdge[28]-0.3093592167691142*coeff[0]*fEdge[28]-0.5303300858899104*coeff[1]*fSkin[19]+0.9951052080056656*coeff[0]*fSkin[19]-0.5303300858899104*coeff[1]*fEdge[19]-0.38273277230987135*coeff[0]*fEdge[19]; 
  edgeSurf_incr[29] = 0.6123724356957944*coeff[1]*fSkin[29]-1.0164659979556616*coeff[0]*fSkin[29]-0.6123724356957944*coeff[1]*fEdge[29]-0.3093592167691142*coeff[0]*fEdge[29]-0.5303300858899104*coeff[1]*fSkin[21]+0.9951052080056656*coeff[0]*fSkin[21]-0.5303300858899104*coeff[1]*fEdge[21]-0.38273277230987135*coeff[0]*fEdge[21]; 
  edgeSurf_incr[30] = 0.6123724356957944*coeff[1]*fSkin[30]-1.0164659979556616*coeff[0]*fSkin[30]-0.6123724356957944*coeff[1]*fEdge[30]-0.3093592167691142*coeff[0]*fEdge[30]-0.5303300858899104*coeff[1]*fSkin[22]+0.9951052080056656*coeff[0]*fSkin[22]-0.5303300858899104*coeff[1]*fEdge[22]-0.38273277230987135*coeff[0]*fEdge[22]; 
  edgeSurf_incr[31] = 0.3827327723098714*coeff[0]*fSkin[40]+0.3827327723098714*coeff[0]*fEdge[40]-0.39774756441743275*coeff[0]*fSkin[31]+0.39774756441743275*coeff[0]*fEdge[31]; 
  edgeSurf_incr[32] = 0.6123724356957944*coeff[1]*fSkin[32]-1.0164659979556616*coeff[0]*fSkin[32]-0.6123724356957944*coeff[1]*fEdge[32]-0.3093592167691142*coeff[0]*fEdge[32]-0.5303300858899104*coeff[1]*fSkin[23]+0.9951052080056656*coeff[0]*fSkin[23]-0.5303300858899104*coeff[1]*fEdge[23]-0.38273277230987135*coeff[0]*fEdge[23]; 
  edgeSurf_incr[33] = 0.3827327723098714*coeff[0]*fSkin[41]+0.3827327723098714*coeff[0]*fEdge[41]-0.39774756441743275*coeff[0]*fSkin[33]+0.39774756441743275*coeff[0]*fEdge[33]; 
  edgeSurf_incr[34] = 0.6123724356957944*coeff[1]*fSkin[34]-1.0164659979556616*coeff[0]*fSkin[34]-0.6123724356957944*coeff[1]*fEdge[34]-0.3093592167691142*coeff[0]*fEdge[34]-0.5303300858899104*coeff[1]*fSkin[25]+0.9951052080056656*coeff[0]*fSkin[25]-0.5303300858899104*coeff[1]*fEdge[25]-0.38273277230987135*coeff[0]*fEdge[25]; 
  edgeSurf_incr[35] = 0.6123724356957944*coeff[1]*fSkin[35]-1.0164659979556616*coeff[0]*fSkin[35]-0.6123724356957944*coeff[1]*fEdge[35]-0.3093592167691142*coeff[0]*fEdge[35]-0.5303300858899104*coeff[1]*fSkin[26]+0.9951052080056656*coeff[0]*fSkin[26]-0.5303300858899104*coeff[1]*fEdge[26]-0.38273277230987135*coeff[0]*fEdge[26]; 
  edgeSurf_incr[36] = 0.3827327723098714*coeff[0]*fSkin[42]+0.3827327723098714*coeff[0]*fEdge[42]-0.39774756441743275*coeff[0]*fSkin[36]+0.39774756441743275*coeff[0]*fEdge[36]; 
  edgeSurf_incr[37] = 0.38273277230987135*coeff[0]*fSkin[43]+0.38273277230987135*coeff[0]*fEdge[43]-0.39774756441743275*coeff[0]*fSkin[37]+0.39774756441743275*coeff[0]*fEdge[37]; 
  edgeSurf_incr[38] = 0.38273277230987135*coeff[0]*fSkin[45]+0.38273277230987135*coeff[0]*fEdge[45]-0.39774756441743275*coeff[0]*fSkin[38]+0.39774756441743275*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = 0.38273277230987135*coeff[0]*fSkin[47]+0.38273277230987135*coeff[0]*fEdge[47]-0.39774756441743275*coeff[0]*fSkin[39]+0.39774756441743275*coeff[0]*fEdge[39]; 
  edgeSurf_incr[40] = 0.6123724356957944*coeff[1]*fSkin[40]-1.0164659979556616*coeff[0]*fSkin[40]-0.6123724356957944*coeff[1]*fEdge[40]-0.3093592167691142*coeff[0]*fEdge[40]-0.5303300858899104*coeff[1]*fSkin[31]+0.9951052080056656*coeff[0]*fSkin[31]-0.5303300858899104*coeff[1]*fEdge[31]-0.38273277230987135*coeff[0]*fEdge[31]; 
  edgeSurf_incr[41] = 0.6123724356957944*coeff[1]*fSkin[41]-1.0164659979556616*coeff[0]*fSkin[41]-0.6123724356957944*coeff[1]*fEdge[41]-0.3093592167691142*coeff[0]*fEdge[41]-0.5303300858899104*coeff[1]*fSkin[33]+0.9951052080056656*coeff[0]*fSkin[33]-0.5303300858899104*coeff[1]*fEdge[33]-0.38273277230987135*coeff[0]*fEdge[33]; 
  edgeSurf_incr[42] = 0.6123724356957944*coeff[1]*fSkin[42]-1.0164659979556616*coeff[0]*fSkin[42]-0.6123724356957944*coeff[1]*fEdge[42]-0.3093592167691142*coeff[0]*fEdge[42]-0.5303300858899104*coeff[1]*fSkin[36]+0.9951052080056656*coeff[0]*fSkin[36]-0.5303300858899104*coeff[1]*fEdge[36]-0.38273277230987135*coeff[0]*fEdge[36]; 
  edgeSurf_incr[43] = 0.6123724356957944*coeff[1]*fSkin[43]-1.0164659979556616*coeff[0]*fSkin[43]-0.6123724356957944*coeff[1]*fEdge[43]-0.3093592167691142*coeff[0]*fEdge[43]-0.5303300858899105*coeff[1]*fSkin[37]+0.9951052080056654*coeff[0]*fSkin[37]-0.5303300858899105*coeff[1]*fEdge[37]-0.3827327723098711*coeff[0]*fEdge[37]; 
  edgeSurf_incr[44] = 0.38273277230987135*coeff[0]*fSkin[49]+0.38273277230987135*coeff[0]*fEdge[49]-0.39774756441743275*coeff[0]*fSkin[44]+0.39774756441743275*coeff[0]*fEdge[44]; 
  edgeSurf_incr[45] = 0.6123724356957944*coeff[1]*fSkin[45]-1.0164659979556616*coeff[0]*fSkin[45]-0.6123724356957944*coeff[1]*fEdge[45]-0.3093592167691142*coeff[0]*fEdge[45]-0.5303300858899105*coeff[1]*fSkin[38]+0.9951052080056654*coeff[0]*fSkin[38]-0.5303300858899105*coeff[1]*fEdge[38]-0.3827327723098711*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = 0.38273277230987135*coeff[0]*fSkin[50]+0.38273277230987135*coeff[0]*fEdge[50]-0.39774756441743275*coeff[0]*fSkin[46]+0.39774756441743275*coeff[0]*fEdge[46]; 
  edgeSurf_incr[47] = 0.6123724356957944*coeff[1]*fSkin[47]-1.0164659979556616*coeff[0]*fSkin[47]-0.6123724356957944*coeff[1]*fEdge[47]-0.3093592167691142*coeff[0]*fEdge[47]-0.5303300858899105*coeff[1]*fSkin[39]+0.9951052080056654*coeff[0]*fSkin[39]-0.5303300858899105*coeff[1]*fEdge[39]-0.3827327723098711*coeff[0]*fEdge[39]; 
  edgeSurf_incr[48] = 0.38273277230987135*coeff[0]*fSkin[51]+0.38273277230987135*coeff[0]*fEdge[51]-0.39774756441743275*coeff[0]*fSkin[48]+0.39774756441743275*coeff[0]*fEdge[48]; 
  edgeSurf_incr[49] = 0.6123724356957944*coeff[1]*fSkin[49]-1.0164659979556616*coeff[0]*fSkin[49]-0.6123724356957944*coeff[1]*fEdge[49]-0.3093592167691142*coeff[0]*fEdge[49]-0.5303300858899105*coeff[1]*fSkin[44]+0.9951052080056654*coeff[0]*fSkin[44]-0.5303300858899105*coeff[1]*fEdge[44]-0.3827327723098711*coeff[0]*fEdge[44]; 
  edgeSurf_incr[50] = 0.6123724356957944*coeff[1]*fSkin[50]-1.0164659979556616*coeff[0]*fSkin[50]-0.6123724356957944*coeff[1]*fEdge[50]-0.3093592167691142*coeff[0]*fEdge[50]-0.5303300858899105*coeff[1]*fSkin[46]+0.9951052080056654*coeff[0]*fSkin[46]-0.5303300858899105*coeff[1]*fEdge[46]-0.3827327723098711*coeff[0]*fEdge[46]; 
  edgeSurf_incr[51] = 0.6123724356957944*coeff[1]*fSkin[51]-1.0164659979556616*coeff[0]*fSkin[51]-0.6123724356957944*coeff[1]*fEdge[51]-0.3093592167691142*coeff[0]*fEdge[51]-0.5303300858899105*coeff[1]*fSkin[48]+0.9951052080056654*coeff[0]*fSkin[48]-0.5303300858899105*coeff[1]*fEdge[48]-0.3827327723098711*coeff[0]*fEdge[48]; 
  edgeSurf_incr[52] = 0.3827327723098714*coeff[0]*fSkin[53]+0.3827327723098714*coeff[0]*fEdge[53]-0.39774756441743275*coeff[0]*fSkin[52]+0.39774756441743275*coeff[0]*fEdge[52]; 
  edgeSurf_incr[53] = 0.6123724356957944*coeff[1]*fSkin[53]-1.0164659979556616*coeff[0]*fSkin[53]-0.6123724356957944*coeff[1]*fEdge[53]-0.3093592167691142*coeff[0]*fEdge[53]-0.5303300858899104*coeff[1]*fSkin[52]+0.9951052080056656*coeff[0]*fSkin[52]-0.5303300858899104*coeff[1]*fEdge[52]-0.38273277230987135*coeff[0]*fEdge[52]; 

  boundSurf_incr[1] = -(1.224744871391589*coeff[1]*fSkin[1])-0.7071067811865475*coeff[0]*fSkin[1]-1.060660171779821*fSkin[0]*coeff[1]-0.6123724356957944*coeff[0]*fSkin[0]; 
  boundSurf_incr[5] = -(1.224744871391589*coeff[1]*fSkin[5])-0.7071067811865475*coeff[0]*fSkin[5]-1.060660171779821*coeff[1]*fSkin[2]-0.6123724356957944*coeff[0]*fSkin[2]; 
  boundSurf_incr[6] = -(1.224744871391589*coeff[1]*fSkin[6])-0.7071067811865475*coeff[0]*fSkin[6]-1.060660171779821*coeff[1]*fSkin[3]-0.6123724356957944*coeff[0]*fSkin[3]; 
  boundSurf_incr[8] = -(1.224744871391589*coeff[1]*fSkin[8])-0.7071067811865475*coeff[0]*fSkin[8]-1.060660171779821*coeff[1]*fSkin[4]-0.6123724356957944*coeff[0]*fSkin[4]; 
  boundSurf_incr[14] = -(1.224744871391589*coeff[1]*fSkin[14])-0.7071067811865475*coeff[0]*fSkin[14]-1.060660171779821*coeff[1]*fSkin[7]-0.6123724356957944*coeff[0]*fSkin[7]; 
  boundSurf_incr[15] = -(1.224744871391589*coeff[1]*fSkin[15])-0.7071067811865475*coeff[0]*fSkin[15]-1.060660171779821*coeff[1]*fSkin[9]-0.6123724356957944*coeff[0]*fSkin[9]; 
  boundSurf_incr[16] = -(1.224744871391589*coeff[1]*fSkin[16])-0.7071067811865475*coeff[0]*fSkin[16]-1.060660171779821*coeff[1]*fSkin[10]-0.6123724356957944*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = -(1.224744871391589*coeff[1]*fSkin[18])-0.7071067811865475*coeff[0]*fSkin[18]-1.060660171779821*coeff[1]*fSkin[11]-0.6123724356957944*coeff[0]*fSkin[11]; 
  boundSurf_incr[20] = -(1.224744871391589*coeff[1]*fSkin[20])-0.7071067811865475*coeff[0]*fSkin[20]-1.060660171779821*coeff[1]*fSkin[12]-0.6123724356957944*coeff[0]*fSkin[12]; 
  boundSurf_incr[24] = -(1.224744871391589*coeff[1]*fSkin[24])-0.7071067811865475*coeff[0]*fSkin[24]-1.060660171779821*coeff[1]*fSkin[13]-0.6123724356957944*coeff[0]*fSkin[13]; 
  boundSurf_incr[27] = -(1.224744871391589*coeff[1]*fSkin[27])-0.7071067811865475*coeff[0]*fSkin[27]-1.060660171779821*coeff[1]*fSkin[17]-0.6123724356957944*coeff[0]*fSkin[17]; 
  boundSurf_incr[28] = -(1.224744871391589*coeff[1]*fSkin[28])-0.7071067811865475*coeff[0]*fSkin[28]-1.060660171779821*coeff[1]*fSkin[19]-0.6123724356957944*coeff[0]*fSkin[19]; 
  boundSurf_incr[29] = -(1.224744871391589*coeff[1]*fSkin[29])-0.7071067811865475*coeff[0]*fSkin[29]-1.060660171779821*coeff[1]*fSkin[21]-0.6123724356957944*coeff[0]*fSkin[21]; 
  boundSurf_incr[30] = -(1.224744871391589*coeff[1]*fSkin[30])-0.7071067811865475*coeff[0]*fSkin[30]-1.060660171779821*coeff[1]*fSkin[22]-0.6123724356957944*coeff[0]*fSkin[22]; 
  boundSurf_incr[32] = -(1.224744871391589*coeff[1]*fSkin[32])-0.7071067811865475*coeff[0]*fSkin[32]-1.060660171779821*coeff[1]*fSkin[23]-0.6123724356957944*coeff[0]*fSkin[23]; 
  boundSurf_incr[34] = -(1.224744871391589*coeff[1]*fSkin[34])-0.7071067811865475*coeff[0]*fSkin[34]-1.060660171779821*coeff[1]*fSkin[25]-0.6123724356957944*coeff[0]*fSkin[25]; 
  boundSurf_incr[35] = -(1.224744871391589*coeff[1]*fSkin[35])-0.7071067811865475*coeff[0]*fSkin[35]-1.060660171779821*coeff[1]*fSkin[26]-0.6123724356957944*coeff[0]*fSkin[26]; 
  boundSurf_incr[40] = -(1.224744871391589*coeff[1]*fSkin[40])-0.7071067811865475*coeff[0]*fSkin[40]-1.060660171779821*coeff[1]*fSkin[31]-0.6123724356957944*coeff[0]*fSkin[31]; 
  boundSurf_incr[41] = -(1.224744871391589*coeff[1]*fSkin[41])-0.7071067811865475*coeff[0]*fSkin[41]-1.060660171779821*coeff[1]*fSkin[33]-0.6123724356957944*coeff[0]*fSkin[33]; 
  boundSurf_incr[42] = -(1.224744871391589*coeff[1]*fSkin[42])-0.7071067811865475*coeff[0]*fSkin[42]-1.060660171779821*coeff[1]*fSkin[36]-0.6123724356957944*coeff[0]*fSkin[36]; 
  boundSurf_incr[43] = -(1.224744871391589*coeff[1]*fSkin[43])-0.7071067811865475*coeff[0]*fSkin[43]-1.060660171779821*coeff[1]*fSkin[37]-0.6123724356957944*coeff[0]*fSkin[37]; 
  boundSurf_incr[45] = -(1.224744871391589*coeff[1]*fSkin[45])-0.7071067811865475*coeff[0]*fSkin[45]-1.060660171779821*coeff[1]*fSkin[38]-0.6123724356957944*coeff[0]*fSkin[38]; 
  boundSurf_incr[47] = -(1.224744871391589*coeff[1]*fSkin[47])-0.7071067811865475*coeff[0]*fSkin[47]-1.060660171779821*coeff[1]*fSkin[39]-0.6123724356957944*coeff[0]*fSkin[39]; 
  boundSurf_incr[49] = -(1.224744871391589*coeff[1]*fSkin[49])-0.7071067811865475*coeff[0]*fSkin[49]-1.060660171779821*coeff[1]*fSkin[44]-0.6123724356957944*coeff[0]*fSkin[44]; 
  boundSurf_incr[50] = -(1.224744871391589*coeff[1]*fSkin[50])-0.7071067811865475*coeff[0]*fSkin[50]-1.060660171779821*coeff[1]*fSkin[46]-0.6123724356957944*coeff[0]*fSkin[46]; 
  boundSurf_incr[51] = -(1.224744871391589*coeff[1]*fSkin[51])-0.7071067811865475*coeff[0]*fSkin[51]-1.060660171779821*coeff[1]*fSkin[48]-0.6123724356957944*coeff[0]*fSkin[48]; 
  boundSurf_incr[53] = -(1.224744871391589*coeff[1]*fSkin[53])-0.7071067811865475*coeff[0]*fSkin[53]-1.060660171779821*coeff[1]*fSkin[52]-0.6123724356957944*coeff[0]*fSkin[52]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac; 
  out[24] += (vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac; 
  out[25] += (vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac; 
  out[26] += (vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac; 
  out[27] += (vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*Jfac; 
  out[28] += (vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*Jfac; 
  out[29] += (vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*Jfac; 
  out[30] += (vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*Jfac; 
  out[31] += (vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*Jfac; 
  out[32] += (vol_incr[32]+edgeSurf_incr[32]+boundSurf_incr[32])*Jfac; 
  out[33] += (vol_incr[33]+edgeSurf_incr[33]+boundSurf_incr[33])*Jfac; 
  out[34] += (vol_incr[34]+edgeSurf_incr[34]+boundSurf_incr[34])*Jfac; 
  out[35] += (vol_incr[35]+edgeSurf_incr[35]+boundSurf_incr[35])*Jfac; 
  out[36] += (vol_incr[36]+edgeSurf_incr[36]+boundSurf_incr[36])*Jfac; 
  out[37] += (vol_incr[37]+edgeSurf_incr[37]+boundSurf_incr[37])*Jfac; 
  out[38] += (vol_incr[38]+edgeSurf_incr[38]+boundSurf_incr[38])*Jfac; 
  out[39] += (vol_incr[39]+edgeSurf_incr[39]+boundSurf_incr[39])*Jfac; 
  out[40] += (vol_incr[40]+edgeSurf_incr[40]+boundSurf_incr[40])*Jfac; 
  out[41] += (vol_incr[41]+edgeSurf_incr[41]+boundSurf_incr[41])*Jfac; 
  out[42] += (vol_incr[42]+edgeSurf_incr[42]+boundSurf_incr[42])*Jfac; 
  out[43] += (vol_incr[43]+edgeSurf_incr[43]+boundSurf_incr[43])*Jfac; 
  out[44] += (vol_incr[44]+edgeSurf_incr[44]+boundSurf_incr[44])*Jfac; 
  out[45] += (vol_incr[45]+edgeSurf_incr[45]+boundSurf_incr[45])*Jfac; 
  out[46] += (vol_incr[46]+edgeSurf_incr[46]+boundSurf_incr[46])*Jfac; 
  out[47] += (vol_incr[47]+edgeSurf_incr[47]+boundSurf_incr[47])*Jfac; 
  out[48] += (vol_incr[48]+edgeSurf_incr[48]+boundSurf_incr[48])*Jfac; 
  out[49] += (vol_incr[49]+edgeSurf_incr[49]+boundSurf_incr[49])*Jfac; 
  out[50] += (vol_incr[50]+edgeSurf_incr[50]+boundSurf_incr[50])*Jfac; 
  out[51] += (vol_incr[51]+edgeSurf_incr[51]+boundSurf_incr[51])*Jfac; 
  out[52] += (vol_incr[52]+edgeSurf_incr[52]+boundSurf_incr[52])*Jfac; 
  out[53] += (vol_incr[53]+edgeSurf_incr[53]+boundSurf_incr[53])*Jfac; 

  return 0.;
}

