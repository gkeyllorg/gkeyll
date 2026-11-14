#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order2_boundary_surfy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  double vol_incr[36] = {0.0}; 

  double edgeSurf_incr[36] = {0.0}; 
  double boundSurf_incr[36] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.5412658773652741*coeff[1]*fSkin[2])-0.5412658773652741*coeff[1]*fEdge[2]-0.5625*fSkin[0]*coeff[1]+0.5625*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(0.5412658773652741*coeff[1]*fSkin[5])-0.5412658773652741*coeff[1]*fEdge[5]-0.5625*coeff[1]*fSkin[1]+0.5625*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(1.4375*coeff[1]*fSkin[2])-0.4375*coeff[1]*fEdge[2]-1.4072912811497125*fSkin[0]*coeff[1]+0.5412658773652739*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(0.5412658773652741*coeff[1]*fSkin[7])-0.5412658773652741*coeff[1]*fEdge[7]-0.5625*coeff[1]*fSkin[3]+0.5625*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = -(0.5412658773652741*coeff[1]*fSkin[9])-0.5412658773652741*coeff[1]*fEdge[9]-0.5625*coeff[1]*fSkin[4]+0.5625*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = -(1.4375*coeff[1]*fSkin[5])-0.4375*coeff[1]*fEdge[5]-1.4072912811497125*coeff[1]*fSkin[1]+0.5412658773652739*coeff[1]*fEdge[1]; 
  edgeSurf_incr[6] = -(0.5412658773652741*coeff[1]*fSkin[13])-0.5412658773652741*coeff[1]*fEdge[13]-0.5625*coeff[1]*fSkin[6]+0.5625*coeff[1]*fEdge[6]; 
  edgeSurf_incr[7] = -(1.4375*coeff[1]*fSkin[7])-0.4375*coeff[1]*fEdge[7]-1.4072912811497125*coeff[1]*fSkin[3]+0.5412658773652739*coeff[1]*fEdge[3]; 
  edgeSurf_incr[8] = -(0.5412658773652741*coeff[1]*fSkin[14])-0.5412658773652741*coeff[1]*fEdge[14]-0.5625*coeff[1]*fSkin[8]+0.5625*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = -(1.4375*coeff[1]*fSkin[9])-0.4375*coeff[1]*fEdge[9]-1.4072912811497125*coeff[1]*fSkin[4]+0.5412658773652739*coeff[1]*fEdge[4]; 
  edgeSurf_incr[10] = -(0.5412658773652741*coeff[1]*fSkin[16])-0.5412658773652741*coeff[1]*fEdge[16]-0.5625*coeff[1]*fSkin[10]+0.5625*coeff[1]*fEdge[10]; 
  edgeSurf_incr[11] = -(0.5412658773652742*coeff[1]*fSkin[18])-0.5412658773652742*coeff[1]*fEdge[18]-0.5625*coeff[1]*fSkin[11]+0.5625*coeff[1]*fEdge[11]; 
  edgeSurf_incr[12] = -(0.5412658773652742*coeff[1]*fSkin[21])-0.5412658773652742*coeff[1]*fEdge[21]-0.5625*coeff[1]*fSkin[12]+0.5625*coeff[1]*fEdge[12]; 
  edgeSurf_incr[13] = -(1.4375*coeff[1]*fSkin[13])-0.4375*coeff[1]*fEdge[13]-1.4072912811497125*coeff[1]*fSkin[6]+0.5412658773652739*coeff[1]*fEdge[6]; 
  edgeSurf_incr[14] = -(1.4375*coeff[1]*fSkin[14])-0.4375*coeff[1]*fEdge[14]-1.4072912811497125*coeff[1]*fSkin[8]+0.5412658773652739*coeff[1]*fEdge[8]; 
  edgeSurf_incr[15] = -(0.5412658773652741*coeff[1]*fSkin[23])-0.5412658773652741*coeff[1]*fEdge[23]-0.5625*coeff[1]*fSkin[15]+0.5625*coeff[1]*fEdge[15]; 
  edgeSurf_incr[16] = -(1.4375*coeff[1]*fSkin[16])-0.4375*coeff[1]*fEdge[16]-1.4072912811497125*coeff[1]*fSkin[10]+0.5412658773652739*coeff[1]*fEdge[10]; 
  edgeSurf_incr[17] = -(0.5412658773652742*coeff[1]*fSkin[24])-0.5412658773652742*coeff[1]*fEdge[24]-0.5625*coeff[1]*fSkin[17]+0.5625*coeff[1]*fEdge[17]; 
  edgeSurf_incr[18] = -(1.4375*coeff[1]*fSkin[18])-0.4375*coeff[1]*fEdge[18]-1.4072912811497127*coeff[1]*fSkin[11]+0.5412658773652742*coeff[1]*fEdge[11]; 
  edgeSurf_incr[19] = -(0.5412658773652742*coeff[1]*fSkin[26])-0.5412658773652742*coeff[1]*fEdge[26]-0.5625*coeff[1]*fSkin[19]+0.5625*coeff[1]*fEdge[19]; 
  edgeSurf_incr[20] = -(0.5412658773652742*coeff[1]*fSkin[27])-0.5412658773652742*coeff[1]*fEdge[27]-0.5625*coeff[1]*fSkin[20]+0.5625*coeff[1]*fEdge[20]; 
  edgeSurf_incr[21] = -(1.4375*coeff[1]*fSkin[21])-0.4375*coeff[1]*fEdge[21]-1.4072912811497127*coeff[1]*fSkin[12]+0.5412658773652742*coeff[1]*fEdge[12]; 
  edgeSurf_incr[22] = -(0.5412658773652742*coeff[1]*fSkin[29])-0.5412658773652742*coeff[1]*fEdge[29]-0.5625*coeff[1]*fSkin[22]+0.5625*coeff[1]*fEdge[22]; 
  edgeSurf_incr[23] = -(1.4375*coeff[1]*fSkin[23])-0.4375*coeff[1]*fEdge[23]-1.4072912811497125*coeff[1]*fSkin[15]+0.5412658773652739*coeff[1]*fEdge[15]; 
  edgeSurf_incr[24] = -(1.4375*coeff[1]*fSkin[24])-0.4375*coeff[1]*fEdge[24]-1.4072912811497127*coeff[1]*fSkin[17]+0.5412658773652742*coeff[1]*fEdge[17]; 
  edgeSurf_incr[25] = -(0.5412658773652742*coeff[1]*fSkin[31])-0.5412658773652742*coeff[1]*fEdge[31]-0.5625*coeff[1]*fSkin[25]+0.5625*coeff[1]*fEdge[25]; 
  edgeSurf_incr[26] = -(1.4375*coeff[1]*fSkin[26])-0.4375*coeff[1]*fEdge[26]-1.4072912811497127*coeff[1]*fSkin[19]+0.5412658773652742*coeff[1]*fEdge[19]; 
  edgeSurf_incr[27] = -(1.4375*coeff[1]*fSkin[27])-0.4375*coeff[1]*fEdge[27]-1.4072912811497127*coeff[1]*fSkin[20]+0.5412658773652742*coeff[1]*fEdge[20]; 
  edgeSurf_incr[28] = -(0.5412658773652742*coeff[1]*fSkin[32])-0.5412658773652742*coeff[1]*fEdge[32]-0.5625*coeff[1]*fSkin[28]+0.5625*coeff[1]*fEdge[28]; 
  edgeSurf_incr[29] = -(1.4375*coeff[1]*fSkin[29])-0.4375*coeff[1]*fEdge[29]-1.4072912811497127*coeff[1]*fSkin[22]+0.5412658773652742*coeff[1]*fEdge[22]; 
  edgeSurf_incr[30] = -(0.5412658773652741*coeff[1]*fSkin[34])-0.5412658773652741*coeff[1]*fEdge[34]-0.5625*coeff[1]*fSkin[30]+0.5625*coeff[1]*fEdge[30]; 
  edgeSurf_incr[31] = -(1.4375*coeff[1]*fSkin[31])-0.4375*coeff[1]*fEdge[31]-1.4072912811497127*coeff[1]*fSkin[25]+0.5412658773652742*coeff[1]*fEdge[25]; 
  edgeSurf_incr[32] = -(1.4375*coeff[1]*fSkin[32])-0.4375*coeff[1]*fEdge[32]-1.4072912811497127*coeff[1]*fSkin[28]+0.5412658773652742*coeff[1]*fEdge[28]; 
  edgeSurf_incr[33] = -(0.5412658773652741*coeff[1]*fSkin[35])-0.5412658773652741*coeff[1]*fEdge[35]-0.5625*coeff[1]*fSkin[33]+0.5625*coeff[1]*fEdge[33]; 
  edgeSurf_incr[34] = -(1.4375*coeff[1]*fSkin[34])-0.4375*coeff[1]*fEdge[34]-1.4072912811497125*coeff[1]*fSkin[30]+0.5412658773652739*coeff[1]*fEdge[30]; 
  edgeSurf_incr[35] = -(1.4375*coeff[1]*fSkin[35])-0.4375*coeff[1]*fEdge[35]-1.4072912811497125*coeff[1]*fSkin[33]+0.5412658773652739*coeff[1]*fEdge[33]; 

  boundSurf_incr[2] = 0.8660254037844386*fSkin[0]*coeff[1]-1.0*coeff[1]*fSkin[2]; 
  boundSurf_incr[5] = 0.8660254037844386*coeff[1]*fSkin[1]-1.0*coeff[1]*fSkin[5]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[1]*fSkin[3]-1.0*coeff[1]*fSkin[7]; 
  boundSurf_incr[9] = 0.8660254037844386*coeff[1]*fSkin[4]-1.0*coeff[1]*fSkin[9]; 
  boundSurf_incr[13] = 0.8660254037844386*coeff[1]*fSkin[6]-1.0*coeff[1]*fSkin[13]; 
  boundSurf_incr[14] = 0.8660254037844386*coeff[1]*fSkin[8]-1.0*coeff[1]*fSkin[14]; 
  boundSurf_incr[16] = 0.8660254037844386*coeff[1]*fSkin[10]-1.0*coeff[1]*fSkin[16]; 
  boundSurf_incr[18] = 0.8660254037844387*coeff[1]*fSkin[11]-1.0*coeff[1]*fSkin[18]; 
  boundSurf_incr[21] = 0.8660254037844387*coeff[1]*fSkin[12]-1.0*coeff[1]*fSkin[21]; 
  boundSurf_incr[23] = 0.8660254037844386*coeff[1]*fSkin[15]-1.0*coeff[1]*fSkin[23]; 
  boundSurf_incr[24] = 0.8660254037844387*coeff[1]*fSkin[17]-1.0*coeff[1]*fSkin[24]; 
  boundSurf_incr[26] = 0.8660254037844387*coeff[1]*fSkin[19]-1.0*coeff[1]*fSkin[26]; 
  boundSurf_incr[27] = 0.8660254037844387*coeff[1]*fSkin[20]-1.0*coeff[1]*fSkin[27]; 
  boundSurf_incr[29] = 0.8660254037844387*coeff[1]*fSkin[22]-1.0*coeff[1]*fSkin[29]; 
  boundSurf_incr[31] = 0.8660254037844387*coeff[1]*fSkin[25]-1.0*coeff[1]*fSkin[31]; 
  boundSurf_incr[32] = 0.8660254037844387*coeff[1]*fSkin[28]-1.0*coeff[1]*fSkin[32]; 
  boundSurf_incr[34] = 0.8660254037844386*coeff[1]*fSkin[30]-1.0*coeff[1]*fSkin[34]; 
  boundSurf_incr[35] = 0.8660254037844386*coeff[1]*fSkin[33]-1.0*coeff[1]*fSkin[35]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[1]*fSkin[2]+0.5412658773652741*coeff[1]*fEdge[2]-0.5625*fSkin[0]*coeff[1]+0.5625*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 0.5412658773652741*coeff[1]*fSkin[5]+0.5412658773652741*coeff[1]*fEdge[5]-0.5625*coeff[1]*fSkin[1]+0.5625*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(1.4375*coeff[1]*fSkin[2])-0.4375*coeff[1]*fEdge[2]+1.4072912811497125*fSkin[0]*coeff[1]-0.5412658773652739*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 0.5412658773652741*coeff[1]*fSkin[7]+0.5412658773652741*coeff[1]*fEdge[7]-0.5625*coeff[1]*fSkin[3]+0.5625*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = 0.5412658773652741*coeff[1]*fSkin[9]+0.5412658773652741*coeff[1]*fEdge[9]-0.5625*coeff[1]*fSkin[4]+0.5625*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = -(1.4375*coeff[1]*fSkin[5])-0.4375*coeff[1]*fEdge[5]+1.4072912811497125*coeff[1]*fSkin[1]-0.5412658773652739*coeff[1]*fEdge[1]; 
  edgeSurf_incr[6] = 0.5412658773652741*coeff[1]*fSkin[13]+0.5412658773652741*coeff[1]*fEdge[13]-0.5625*coeff[1]*fSkin[6]+0.5625*coeff[1]*fEdge[6]; 
  edgeSurf_incr[7] = -(1.4375*coeff[1]*fSkin[7])-0.4375*coeff[1]*fEdge[7]+1.4072912811497125*coeff[1]*fSkin[3]-0.5412658773652739*coeff[1]*fEdge[3]; 
  edgeSurf_incr[8] = 0.5412658773652741*coeff[1]*fSkin[14]+0.5412658773652741*coeff[1]*fEdge[14]-0.5625*coeff[1]*fSkin[8]+0.5625*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = -(1.4375*coeff[1]*fSkin[9])-0.4375*coeff[1]*fEdge[9]+1.4072912811497125*coeff[1]*fSkin[4]-0.5412658773652739*coeff[1]*fEdge[4]; 
  edgeSurf_incr[10] = 0.5412658773652741*coeff[1]*fSkin[16]+0.5412658773652741*coeff[1]*fEdge[16]-0.5625*coeff[1]*fSkin[10]+0.5625*coeff[1]*fEdge[10]; 
  edgeSurf_incr[11] = 0.5412658773652742*coeff[1]*fSkin[18]+0.5412658773652742*coeff[1]*fEdge[18]-0.5625*coeff[1]*fSkin[11]+0.5625*coeff[1]*fEdge[11]; 
  edgeSurf_incr[12] = 0.5412658773652742*coeff[1]*fSkin[21]+0.5412658773652742*coeff[1]*fEdge[21]-0.5625*coeff[1]*fSkin[12]+0.5625*coeff[1]*fEdge[12]; 
  edgeSurf_incr[13] = -(1.4375*coeff[1]*fSkin[13])-0.4375*coeff[1]*fEdge[13]+1.4072912811497125*coeff[1]*fSkin[6]-0.5412658773652739*coeff[1]*fEdge[6]; 
  edgeSurf_incr[14] = -(1.4375*coeff[1]*fSkin[14])-0.4375*coeff[1]*fEdge[14]+1.4072912811497125*coeff[1]*fSkin[8]-0.5412658773652739*coeff[1]*fEdge[8]; 
  edgeSurf_incr[15] = 0.5412658773652741*coeff[1]*fSkin[23]+0.5412658773652741*coeff[1]*fEdge[23]-0.5625*coeff[1]*fSkin[15]+0.5625*coeff[1]*fEdge[15]; 
  edgeSurf_incr[16] = -(1.4375*coeff[1]*fSkin[16])-0.4375*coeff[1]*fEdge[16]+1.4072912811497125*coeff[1]*fSkin[10]-0.5412658773652739*coeff[1]*fEdge[10]; 
  edgeSurf_incr[17] = 0.5412658773652742*coeff[1]*fSkin[24]+0.5412658773652742*coeff[1]*fEdge[24]-0.5625*coeff[1]*fSkin[17]+0.5625*coeff[1]*fEdge[17]; 
  edgeSurf_incr[18] = -(1.4375*coeff[1]*fSkin[18])-0.4375*coeff[1]*fEdge[18]+1.4072912811497127*coeff[1]*fSkin[11]-0.5412658773652742*coeff[1]*fEdge[11]; 
  edgeSurf_incr[19] = 0.5412658773652742*coeff[1]*fSkin[26]+0.5412658773652742*coeff[1]*fEdge[26]-0.5625*coeff[1]*fSkin[19]+0.5625*coeff[1]*fEdge[19]; 
  edgeSurf_incr[20] = 0.5412658773652742*coeff[1]*fSkin[27]+0.5412658773652742*coeff[1]*fEdge[27]-0.5625*coeff[1]*fSkin[20]+0.5625*coeff[1]*fEdge[20]; 
  edgeSurf_incr[21] = -(1.4375*coeff[1]*fSkin[21])-0.4375*coeff[1]*fEdge[21]+1.4072912811497127*coeff[1]*fSkin[12]-0.5412658773652742*coeff[1]*fEdge[12]; 
  edgeSurf_incr[22] = 0.5412658773652742*coeff[1]*fSkin[29]+0.5412658773652742*coeff[1]*fEdge[29]-0.5625*coeff[1]*fSkin[22]+0.5625*coeff[1]*fEdge[22]; 
  edgeSurf_incr[23] = -(1.4375*coeff[1]*fSkin[23])-0.4375*coeff[1]*fEdge[23]+1.4072912811497125*coeff[1]*fSkin[15]-0.5412658773652739*coeff[1]*fEdge[15]; 
  edgeSurf_incr[24] = -(1.4375*coeff[1]*fSkin[24])-0.4375*coeff[1]*fEdge[24]+1.4072912811497127*coeff[1]*fSkin[17]-0.5412658773652742*coeff[1]*fEdge[17]; 
  edgeSurf_incr[25] = 0.5412658773652742*coeff[1]*fSkin[31]+0.5412658773652742*coeff[1]*fEdge[31]-0.5625*coeff[1]*fSkin[25]+0.5625*coeff[1]*fEdge[25]; 
  edgeSurf_incr[26] = -(1.4375*coeff[1]*fSkin[26])-0.4375*coeff[1]*fEdge[26]+1.4072912811497127*coeff[1]*fSkin[19]-0.5412658773652742*coeff[1]*fEdge[19]; 
  edgeSurf_incr[27] = -(1.4375*coeff[1]*fSkin[27])-0.4375*coeff[1]*fEdge[27]+1.4072912811497127*coeff[1]*fSkin[20]-0.5412658773652742*coeff[1]*fEdge[20]; 
  edgeSurf_incr[28] = 0.5412658773652742*coeff[1]*fSkin[32]+0.5412658773652742*coeff[1]*fEdge[32]-0.5625*coeff[1]*fSkin[28]+0.5625*coeff[1]*fEdge[28]; 
  edgeSurf_incr[29] = -(1.4375*coeff[1]*fSkin[29])-0.4375*coeff[1]*fEdge[29]+1.4072912811497127*coeff[1]*fSkin[22]-0.5412658773652742*coeff[1]*fEdge[22]; 
  edgeSurf_incr[30] = 0.5412658773652741*coeff[1]*fSkin[34]+0.5412658773652741*coeff[1]*fEdge[34]-0.5625*coeff[1]*fSkin[30]+0.5625*coeff[1]*fEdge[30]; 
  edgeSurf_incr[31] = -(1.4375*coeff[1]*fSkin[31])-0.4375*coeff[1]*fEdge[31]+1.4072912811497127*coeff[1]*fSkin[25]-0.5412658773652742*coeff[1]*fEdge[25]; 
  edgeSurf_incr[32] = -(1.4375*coeff[1]*fSkin[32])-0.4375*coeff[1]*fEdge[32]+1.4072912811497127*coeff[1]*fSkin[28]-0.5412658773652742*coeff[1]*fEdge[28]; 
  edgeSurf_incr[33] = 0.5412658773652741*coeff[1]*fSkin[35]+0.5412658773652741*coeff[1]*fEdge[35]-0.5625*coeff[1]*fSkin[33]+0.5625*coeff[1]*fEdge[33]; 
  edgeSurf_incr[34] = -(1.4375*coeff[1]*fSkin[34])-0.4375*coeff[1]*fEdge[34]+1.4072912811497125*coeff[1]*fSkin[30]-0.5412658773652739*coeff[1]*fEdge[30]; 
  edgeSurf_incr[35] = -(1.4375*coeff[1]*fSkin[35])-0.4375*coeff[1]*fEdge[35]+1.4072912811497125*coeff[1]*fSkin[33]-0.5412658773652739*coeff[1]*fEdge[33]; 

  boundSurf_incr[2] = -(1.0*coeff[1]*fSkin[2])-0.8660254037844386*fSkin[0]*coeff[1]; 
  boundSurf_incr[5] = -(1.0*coeff[1]*fSkin[5])-0.8660254037844386*coeff[1]*fSkin[1]; 
  boundSurf_incr[7] = -(1.0*coeff[1]*fSkin[7])-0.8660254037844386*coeff[1]*fSkin[3]; 
  boundSurf_incr[9] = -(1.0*coeff[1]*fSkin[9])-0.8660254037844386*coeff[1]*fSkin[4]; 
  boundSurf_incr[13] = -(1.0*coeff[1]*fSkin[13])-0.8660254037844386*coeff[1]*fSkin[6]; 
  boundSurf_incr[14] = -(1.0*coeff[1]*fSkin[14])-0.8660254037844386*coeff[1]*fSkin[8]; 
  boundSurf_incr[16] = -(1.0*coeff[1]*fSkin[16])-0.8660254037844386*coeff[1]*fSkin[10]; 
  boundSurf_incr[18] = -(1.0*coeff[1]*fSkin[18])-0.8660254037844387*coeff[1]*fSkin[11]; 
  boundSurf_incr[21] = -(1.0*coeff[1]*fSkin[21])-0.8660254037844387*coeff[1]*fSkin[12]; 
  boundSurf_incr[23] = -(1.0*coeff[1]*fSkin[23])-0.8660254037844386*coeff[1]*fSkin[15]; 
  boundSurf_incr[24] = -(1.0*coeff[1]*fSkin[24])-0.8660254037844387*coeff[1]*fSkin[17]; 
  boundSurf_incr[26] = -(1.0*coeff[1]*fSkin[26])-0.8660254037844387*coeff[1]*fSkin[19]; 
  boundSurf_incr[27] = -(1.0*coeff[1]*fSkin[27])-0.8660254037844387*coeff[1]*fSkin[20]; 
  boundSurf_incr[29] = -(1.0*coeff[1]*fSkin[29])-0.8660254037844387*coeff[1]*fSkin[22]; 
  boundSurf_incr[31] = -(1.0*coeff[1]*fSkin[31])-0.8660254037844387*coeff[1]*fSkin[25]; 
  boundSurf_incr[32] = -(1.0*coeff[1]*fSkin[32])-0.8660254037844387*coeff[1]*fSkin[28]; 
  boundSurf_incr[34] = -(1.0*coeff[1]*fSkin[34])-0.8660254037844386*coeff[1]*fSkin[30]; 
  boundSurf_incr[35] = -(1.0*coeff[1]*fSkin[35])-0.8660254037844386*coeff[1]*fSkin[33]; 

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

  return 0.;
}

GKYL_CU_DH double dg_diffusion_vlasov_order2_boundary_surfy_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  double vol_incr[36] = {0.0}; 
  vol_incr[2] = 1.5*fSkin[1]*coeff[7]+1.5*fSkin[0]*coeff[6]; 
  vol_incr[5] = 1.5*fSkin[0]*coeff[7]+1.5*fSkin[1]*coeff[6]; 
  vol_incr[7] = 1.5*fSkin[6]*coeff[7]+1.5*fSkin[3]*coeff[6]; 
  vol_incr[9] = 1.5*coeff[7]*fSkin[8]+1.5*fSkin[4]*coeff[6]; 
  vol_incr[13] = 1.5*fSkin[3]*coeff[7]+1.5*coeff[6]*fSkin[6]; 
  vol_incr[14] = 1.5*coeff[6]*fSkin[8]+1.5*fSkin[4]*coeff[7]; 
  vol_incr[16] = 1.5*coeff[7]*fSkin[15]+1.5*coeff[6]*fSkin[10]; 
  vol_incr[18] = 1.5*coeff[7]*fSkin[17]+1.5*coeff[6]*fSkin[11]; 
  vol_incr[21] = 1.5*coeff[7]*fSkin[20]+1.5*coeff[6]*fSkin[12]; 
  vol_incr[23] = 1.5*coeff[6]*fSkin[15]+1.5*coeff[7]*fSkin[10]; 
  vol_incr[24] = 1.5*coeff[6]*fSkin[17]+1.5*coeff[7]*fSkin[11]; 
  vol_incr[26] = 1.5*coeff[7]*fSkin[25]+1.5*coeff[6]*fSkin[19]; 
  vol_incr[27] = 1.5*coeff[6]*fSkin[20]+1.5*coeff[7]*fSkin[12]; 
  vol_incr[29] = 1.5*coeff[7]*fSkin[28]+1.5*coeff[6]*fSkin[22]; 
  vol_incr[31] = 1.5*coeff[6]*fSkin[25]+1.5*coeff[7]*fSkin[19]; 
  vol_incr[32] = 1.5*coeff[6]*fSkin[28]+1.5*coeff[7]*fSkin[22]; 
  vol_incr[34] = 1.5*coeff[7]*fSkin[33]+1.5*coeff[6]*fSkin[30]; 
  vol_incr[35] = 1.5*coeff[6]*fSkin[33]+1.5*coeff[7]*fSkin[30]; 

  double edgeSurf_incr[36] = {0.0}; 
  double boundSurf_incr[36] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.27063293868263705*coeff[5]*fSkin[5])-0.27063293868263705*coeff[5]*fEdge[5]-0.28125*fSkin[1]*coeff[5]+0.28125*fEdge[1]*coeff[5]-0.27063293868263705*fSkin[2]*coeff[4]-0.27063293868263705*fEdge[2]*coeff[4]-0.28125*fSkin[0]*coeff[4]+0.28125*fEdge[0]*coeff[4]; 
  edgeSurf_incr[1] = -(0.27063293868263705*coeff[4]*fSkin[5])-0.27063293868263705*coeff[4]*fEdge[5]-0.27063293868263705*fSkin[2]*coeff[5]-0.27063293868263705*fEdge[2]*coeff[5]-0.28125*fSkin[0]*coeff[5]+0.28125*fEdge[0]*coeff[5]-0.28125*fSkin[1]*coeff[4]+0.28125*fEdge[1]*coeff[4]; 
  edgeSurf_incr[2] = -(0.4330127018922193*fSkin[5]*coeff[7])+0.4330127018922193*fEdge[5]*coeff[7]-0.375*fSkin[1]*coeff[7]-0.375*fEdge[1]*coeff[7]-0.4330127018922193*fSkin[2]*coeff[6]+0.4330127018922193*fEdge[2]*coeff[6]-0.375*fSkin[0]*coeff[6]-0.375*fEdge[0]*coeff[6]-0.71875*coeff[5]*fSkin[5]-0.21875*coeff[5]*fEdge[5]-0.7036456405748562*fSkin[1]*coeff[5]+0.27063293868263694*fEdge[1]*coeff[5]-0.71875*fSkin[2]*coeff[4]-0.21875*fEdge[2]*coeff[4]-0.7036456405748562*fSkin[0]*coeff[4]+0.27063293868263694*fEdge[0]*coeff[4]; 
  edgeSurf_incr[3] = -(0.27063293868263705*coeff[5]*fSkin[13])-0.27063293868263705*coeff[5]*fEdge[13]-0.27063293868263705*coeff[4]*fSkin[7]-0.27063293868263705*coeff[4]*fEdge[7]-0.28125*coeff[5]*fSkin[6]+0.28125*coeff[5]*fEdge[6]-0.28125*fSkin[3]*coeff[4]+0.28125*fEdge[3]*coeff[4]; 
  edgeSurf_incr[4] = -(0.27063293868263705*coeff[5]*fSkin[14])-0.27063293868263705*coeff[5]*fEdge[14]-0.27063293868263705*coeff[4]*fSkin[9]-0.27063293868263705*coeff[4]*fEdge[9]-0.28125*coeff[5]*fSkin[8]+0.28125*coeff[5]*fEdge[8]-0.28125*coeff[4]*fSkin[4]+0.28125*coeff[4]*fEdge[4]; 
  edgeSurf_incr[5] = -(0.4330127018922193*fSkin[2]*coeff[7])+0.4330127018922193*fEdge[2]*coeff[7]-0.375*fSkin[0]*coeff[7]-0.375*fEdge[0]*coeff[7]-0.4330127018922193*fSkin[5]*coeff[6]+0.4330127018922193*fEdge[5]*coeff[6]-0.375*fSkin[1]*coeff[6]-0.375*fEdge[1]*coeff[6]-0.71875*coeff[4]*fSkin[5]-0.21875*coeff[4]*fEdge[5]-0.71875*fSkin[2]*coeff[5]-0.21875*fEdge[2]*coeff[5]-0.7036456405748562*fSkin[0]*coeff[5]+0.27063293868263694*fEdge[0]*coeff[5]-0.7036456405748562*fSkin[1]*coeff[4]+0.27063293868263694*fEdge[1]*coeff[4]; 
  edgeSurf_incr[6] = -(0.27063293868263705*coeff[4]*fSkin[13])-0.27063293868263705*coeff[4]*fEdge[13]-0.27063293868263705*coeff[5]*fSkin[7]-0.27063293868263705*coeff[5]*fEdge[7]-0.28125*coeff[4]*fSkin[6]+0.28125*coeff[4]*fEdge[6]-0.28125*fSkin[3]*coeff[5]+0.28125*fEdge[3]*coeff[5]; 
  edgeSurf_incr[7] = -(0.4330127018922193*coeff[7]*fSkin[13])-0.71875*coeff[5]*fSkin[13]+0.4330127018922193*coeff[7]*fEdge[13]-0.21875*coeff[5]*fEdge[13]-0.4330127018922193*coeff[6]*fSkin[7]-0.71875*coeff[4]*fSkin[7]+0.4330127018922193*coeff[6]*fEdge[7]-0.21875*coeff[4]*fEdge[7]-0.375*fSkin[6]*coeff[7]-0.375*fEdge[6]*coeff[7]-0.7036456405748562*coeff[5]*fSkin[6]+0.27063293868263694*coeff[5]*fEdge[6]-0.375*fSkin[3]*coeff[6]-0.375*fEdge[3]*coeff[6]-0.7036456405748562*fSkin[3]*coeff[4]+0.27063293868263694*fEdge[3]*coeff[4]; 
  edgeSurf_incr[8] = -(0.27063293868263705*coeff[4]*fSkin[14])-0.27063293868263705*coeff[4]*fEdge[14]-0.27063293868263705*coeff[5]*fSkin[9]-0.27063293868263705*coeff[5]*fEdge[9]-0.28125*coeff[4]*fSkin[8]+0.28125*coeff[4]*fEdge[8]-0.28125*fSkin[4]*coeff[5]+0.28125*fEdge[4]*coeff[5]; 
  edgeSurf_incr[9] = -(0.4330127018922193*coeff[7]*fSkin[14])-0.71875*coeff[5]*fSkin[14]+0.4330127018922193*coeff[7]*fEdge[14]-0.21875*coeff[5]*fEdge[14]-0.4330127018922193*coeff[6]*fSkin[9]-0.71875*coeff[4]*fSkin[9]+0.4330127018922193*coeff[6]*fEdge[9]-0.21875*coeff[4]*fEdge[9]-0.375*coeff[7]*fSkin[8]-0.7036456405748562*coeff[5]*fSkin[8]-0.375*coeff[7]*fEdge[8]+0.27063293868263694*coeff[5]*fEdge[8]-0.375*fSkin[4]*coeff[6]-0.375*fEdge[4]*coeff[6]-0.7036456405748562*coeff[4]*fSkin[4]+0.27063293868263694*coeff[4]*fEdge[4]; 
  edgeSurf_incr[10] = -(0.27063293868263705*coeff[5]*fSkin[23])-0.27063293868263705*coeff[5]*fEdge[23]-0.27063293868263705*coeff[4]*fSkin[16]-0.27063293868263705*coeff[4]*fEdge[16]-0.28125*coeff[5]*fSkin[15]+0.28125*coeff[5]*fEdge[15]-0.28125*coeff[4]*fSkin[10]+0.28125*coeff[4]*fEdge[10]; 
  edgeSurf_incr[11] = -(0.27063293868263705*coeff[5]*fSkin[24])-0.27063293868263705*coeff[5]*fEdge[24]-0.2706329386826371*coeff[4]*fSkin[18]-0.2706329386826371*coeff[4]*fEdge[18]-0.28124999999999994*coeff[5]*fSkin[17]+0.28124999999999994*coeff[5]*fEdge[17]-0.28125*coeff[4]*fSkin[11]+0.28125*coeff[4]*fEdge[11]; 
  edgeSurf_incr[12] = -(0.27063293868263705*coeff[5]*fSkin[27])-0.27063293868263705*coeff[5]*fEdge[27]-0.2706329386826371*coeff[4]*fSkin[21]-0.2706329386826371*coeff[4]*fEdge[21]-0.28124999999999994*coeff[5]*fSkin[20]+0.28124999999999994*coeff[5]*fEdge[20]-0.28125*coeff[4]*fSkin[12]+0.28125*coeff[4]*fEdge[12]; 
  edgeSurf_incr[13] = -(0.4330127018922193*coeff[6]*fSkin[13])-0.71875*coeff[4]*fSkin[13]+0.4330127018922193*coeff[6]*fEdge[13]-0.21875*coeff[4]*fEdge[13]-0.4330127018922193*coeff[7]*fSkin[7]-0.71875*coeff[5]*fSkin[7]+0.4330127018922193*coeff[7]*fEdge[7]-0.21875*coeff[5]*fEdge[7]-0.375*fSkin[3]*coeff[7]-0.375*fEdge[3]*coeff[7]-0.375*coeff[6]*fSkin[6]-0.7036456405748562*coeff[4]*fSkin[6]-0.375*coeff[6]*fEdge[6]+0.27063293868263694*coeff[4]*fEdge[6]-0.7036456405748562*fSkin[3]*coeff[5]+0.27063293868263694*fEdge[3]*coeff[5]; 
  edgeSurf_incr[14] = -(0.4330127018922193*coeff[6]*fSkin[14])-0.71875*coeff[4]*fSkin[14]+0.4330127018922193*coeff[6]*fEdge[14]-0.21875*coeff[4]*fEdge[14]-0.4330127018922193*coeff[7]*fSkin[9]-0.71875*coeff[5]*fSkin[9]+0.4330127018922193*coeff[7]*fEdge[9]-0.21875*coeff[5]*fEdge[9]-0.375*coeff[6]*fSkin[8]-0.7036456405748562*coeff[4]*fSkin[8]-0.375*coeff[6]*fEdge[8]+0.27063293868263694*coeff[4]*fEdge[8]-0.375*fSkin[4]*coeff[7]-0.375*fEdge[4]*coeff[7]-0.7036456405748562*fSkin[4]*coeff[5]+0.27063293868263694*fEdge[4]*coeff[5]; 
  edgeSurf_incr[15] = -(0.27063293868263705*coeff[4]*fSkin[23])-0.27063293868263705*coeff[4]*fEdge[23]-0.27063293868263705*coeff[5]*fSkin[16]-0.27063293868263705*coeff[5]*fEdge[16]-0.28125*coeff[4]*fSkin[15]+0.28125*coeff[4]*fEdge[15]-0.28125*coeff[5]*fSkin[10]+0.28125*coeff[5]*fEdge[10]; 
  edgeSurf_incr[16] = -(0.4330127018922193*coeff[7]*fSkin[23])-0.71875*coeff[5]*fSkin[23]+0.4330127018922193*coeff[7]*fEdge[23]-0.21875*coeff[5]*fEdge[23]-0.4330127018922193*coeff[6]*fSkin[16]-0.71875*coeff[4]*fSkin[16]+0.4330127018922193*coeff[6]*fEdge[16]-0.21875*coeff[4]*fEdge[16]-0.375*coeff[7]*fSkin[15]-0.7036456405748562*coeff[5]*fSkin[15]-0.375*coeff[7]*fEdge[15]+0.27063293868263694*coeff[5]*fEdge[15]-0.375*coeff[6]*fSkin[10]-0.7036456405748562*coeff[4]*fSkin[10]-0.375*coeff[6]*fEdge[10]+0.27063293868263694*coeff[4]*fEdge[10]; 
  edgeSurf_incr[17] = -(0.2706329386826371*coeff[4]*fSkin[24])-0.2706329386826371*coeff[4]*fEdge[24]-0.27063293868263705*coeff[5]*fSkin[18]-0.27063293868263705*coeff[5]*fEdge[18]-0.28125*coeff[4]*fSkin[17]+0.28125*coeff[4]*fEdge[17]-0.28124999999999994*coeff[5]*fSkin[11]+0.28124999999999994*coeff[5]*fEdge[11]; 
  edgeSurf_incr[18] = -(0.43301270189221935*coeff[7]*fSkin[24])-0.7187500000000001*coeff[5]*fSkin[24]+0.43301270189221935*coeff[7]*fEdge[24]-0.21875*coeff[5]*fEdge[24]-0.4330127018922193*coeff[6]*fSkin[18]-0.71875*coeff[4]*fSkin[18]+0.4330127018922193*coeff[6]*fEdge[18]-0.21875*coeff[4]*fEdge[18]-0.375*coeff[7]*fSkin[17]-0.7036456405748562*coeff[5]*fSkin[17]-0.375*coeff[7]*fEdge[17]+0.27063293868263694*coeff[5]*fEdge[17]-0.375*coeff[6]*fSkin[11]-0.7036456405748563*coeff[4]*fSkin[11]-0.375*coeff[6]*fEdge[11]+0.2706329386826371*coeff[4]*fEdge[11]; 
  edgeSurf_incr[19] = -(0.27063293868263705*coeff[5]*fSkin[31])-0.27063293868263705*coeff[5]*fEdge[31]-0.2706329386826371*coeff[4]*fSkin[26]-0.2706329386826371*coeff[4]*fEdge[26]-0.28124999999999994*coeff[5]*fSkin[25]+0.28124999999999994*coeff[5]*fEdge[25]-0.28125*coeff[4]*fSkin[19]+0.28125*coeff[4]*fEdge[19]; 
  edgeSurf_incr[20] = -(0.2706329386826371*coeff[4]*fSkin[27])-0.2706329386826371*coeff[4]*fEdge[27]-0.27063293868263705*coeff[5]*fSkin[21]-0.27063293868263705*coeff[5]*fEdge[21]-0.28125*coeff[4]*fSkin[20]+0.28125*coeff[4]*fEdge[20]-0.28124999999999994*coeff[5]*fSkin[12]+0.28124999999999994*coeff[5]*fEdge[12]; 
  edgeSurf_incr[21] = -(0.43301270189221935*coeff[7]*fSkin[27])-0.7187500000000001*coeff[5]*fSkin[27]+0.43301270189221935*coeff[7]*fEdge[27]-0.21875*coeff[5]*fEdge[27]-0.4330127018922193*coeff[6]*fSkin[21]-0.71875*coeff[4]*fSkin[21]+0.4330127018922193*coeff[6]*fEdge[21]-0.21875*coeff[4]*fEdge[21]-0.375*coeff[7]*fSkin[20]-0.7036456405748562*coeff[5]*fSkin[20]-0.375*coeff[7]*fEdge[20]+0.27063293868263694*coeff[5]*fEdge[20]-0.375*coeff[6]*fSkin[12]-0.7036456405748563*coeff[4]*fSkin[12]-0.375*coeff[6]*fEdge[12]+0.2706329386826371*coeff[4]*fEdge[12]; 
  edgeSurf_incr[22] = -(0.27063293868263705*coeff[5]*fSkin[32])-0.27063293868263705*coeff[5]*fEdge[32]-0.2706329386826371*coeff[4]*fSkin[29]-0.2706329386826371*coeff[4]*fEdge[29]-0.28124999999999994*coeff[5]*fSkin[28]+0.28124999999999994*coeff[5]*fEdge[28]-0.28125*coeff[4]*fSkin[22]+0.28125*coeff[4]*fEdge[22]; 
  edgeSurf_incr[23] = -(0.4330127018922193*coeff[6]*fSkin[23])-0.71875*coeff[4]*fSkin[23]+0.4330127018922193*coeff[6]*fEdge[23]-0.21875*coeff[4]*fEdge[23]-0.4330127018922193*coeff[7]*fSkin[16]-0.71875*coeff[5]*fSkin[16]+0.4330127018922193*coeff[7]*fEdge[16]-0.21875*coeff[5]*fEdge[16]-0.375*coeff[6]*fSkin[15]-0.7036456405748562*coeff[4]*fSkin[15]-0.375*coeff[6]*fEdge[15]+0.27063293868263694*coeff[4]*fEdge[15]-0.375*coeff[7]*fSkin[10]-0.7036456405748562*coeff[5]*fSkin[10]-0.375*coeff[7]*fEdge[10]+0.27063293868263694*coeff[5]*fEdge[10]; 
  edgeSurf_incr[24] = -(0.4330127018922193*coeff[6]*fSkin[24])-0.71875*coeff[4]*fSkin[24]+0.4330127018922193*coeff[6]*fEdge[24]-0.21875*coeff[4]*fEdge[24]-0.43301270189221935*coeff[7]*fSkin[18]-0.7187500000000001*coeff[5]*fSkin[18]+0.43301270189221935*coeff[7]*fEdge[18]-0.21875*coeff[5]*fEdge[18]-0.375*coeff[6]*fSkin[17]-0.7036456405748563*coeff[4]*fSkin[17]-0.375*coeff[6]*fEdge[17]+0.2706329386826371*coeff[4]*fEdge[17]-0.375*coeff[7]*fSkin[11]-0.7036456405748562*coeff[5]*fSkin[11]-0.375*coeff[7]*fEdge[11]+0.27063293868263694*coeff[5]*fEdge[11]; 
  edgeSurf_incr[25] = -(0.2706329386826371*coeff[4]*fSkin[31])-0.2706329386826371*coeff[4]*fEdge[31]-0.27063293868263705*coeff[5]*fSkin[26]-0.27063293868263705*coeff[5]*fEdge[26]-0.28125*coeff[4]*fSkin[25]+0.28125*coeff[4]*fEdge[25]-0.28124999999999994*coeff[5]*fSkin[19]+0.28124999999999994*coeff[5]*fEdge[19]; 
  edgeSurf_incr[26] = -(0.43301270189221935*coeff[7]*fSkin[31])-0.7187500000000001*coeff[5]*fSkin[31]+0.43301270189221935*coeff[7]*fEdge[31]-0.21875*coeff[5]*fEdge[31]-0.4330127018922193*coeff[6]*fSkin[26]-0.71875*coeff[4]*fSkin[26]+0.4330127018922193*coeff[6]*fEdge[26]-0.21875*coeff[4]*fEdge[26]-0.375*coeff[7]*fSkin[25]-0.7036456405748562*coeff[5]*fSkin[25]-0.375*coeff[7]*fEdge[25]+0.27063293868263694*coeff[5]*fEdge[25]-0.375*coeff[6]*fSkin[19]-0.7036456405748563*coeff[4]*fSkin[19]-0.375*coeff[6]*fEdge[19]+0.2706329386826371*coeff[4]*fEdge[19]; 
  edgeSurf_incr[27] = -(0.4330127018922193*coeff[6]*fSkin[27])-0.71875*coeff[4]*fSkin[27]+0.4330127018922193*coeff[6]*fEdge[27]-0.21875*coeff[4]*fEdge[27]-0.43301270189221935*coeff[7]*fSkin[21]-0.7187500000000001*coeff[5]*fSkin[21]+0.43301270189221935*coeff[7]*fEdge[21]-0.21875*coeff[5]*fEdge[21]-0.375*coeff[6]*fSkin[20]-0.7036456405748563*coeff[4]*fSkin[20]-0.375*coeff[6]*fEdge[20]+0.2706329386826371*coeff[4]*fEdge[20]-0.375*coeff[7]*fSkin[12]-0.7036456405748562*coeff[5]*fSkin[12]-0.375*coeff[7]*fEdge[12]+0.27063293868263694*coeff[5]*fEdge[12]; 
  edgeSurf_incr[28] = -(0.2706329386826371*coeff[4]*fSkin[32])-0.2706329386826371*coeff[4]*fEdge[32]-0.27063293868263705*coeff[5]*fSkin[29]-0.27063293868263705*coeff[5]*fEdge[29]-0.28125*coeff[4]*fSkin[28]+0.28125*coeff[4]*fEdge[28]-0.28124999999999994*coeff[5]*fSkin[22]+0.28124999999999994*coeff[5]*fEdge[22]; 
  edgeSurf_incr[29] = -(0.43301270189221935*coeff[7]*fSkin[32])-0.7187500000000001*coeff[5]*fSkin[32]+0.43301270189221935*coeff[7]*fEdge[32]-0.21875*coeff[5]*fEdge[32]-0.4330127018922193*coeff[6]*fSkin[29]-0.71875*coeff[4]*fSkin[29]+0.4330127018922193*coeff[6]*fEdge[29]-0.21875*coeff[4]*fEdge[29]-0.375*coeff[7]*fSkin[28]-0.7036456405748562*coeff[5]*fSkin[28]-0.375*coeff[7]*fEdge[28]+0.27063293868263694*coeff[5]*fEdge[28]-0.375*coeff[6]*fSkin[22]-0.7036456405748563*coeff[4]*fSkin[22]-0.375*coeff[6]*fEdge[22]+0.2706329386826371*coeff[4]*fEdge[22]; 
  edgeSurf_incr[30] = -(0.27063293868263705*coeff[5]*fSkin[35])-0.27063293868263705*coeff[5]*fEdge[35]-0.27063293868263705*coeff[4]*fSkin[34]-0.27063293868263705*coeff[4]*fEdge[34]-0.28125*coeff[5]*fSkin[33]+0.28125*coeff[5]*fEdge[33]-0.28125*coeff[4]*fSkin[30]+0.28125*coeff[4]*fEdge[30]; 
  edgeSurf_incr[31] = -(0.4330127018922193*coeff[6]*fSkin[31])-0.71875*coeff[4]*fSkin[31]+0.4330127018922193*coeff[6]*fEdge[31]-0.21875*coeff[4]*fEdge[31]-0.43301270189221935*coeff[7]*fSkin[26]-0.7187500000000001*coeff[5]*fSkin[26]+0.43301270189221935*coeff[7]*fEdge[26]-0.21875*coeff[5]*fEdge[26]-0.375*coeff[6]*fSkin[25]-0.7036456405748563*coeff[4]*fSkin[25]-0.375*coeff[6]*fEdge[25]+0.2706329386826371*coeff[4]*fEdge[25]-0.375*coeff[7]*fSkin[19]-0.7036456405748562*coeff[5]*fSkin[19]-0.375*coeff[7]*fEdge[19]+0.27063293868263694*coeff[5]*fEdge[19]; 
  edgeSurf_incr[32] = -(0.4330127018922193*coeff[6]*fSkin[32])-0.71875*coeff[4]*fSkin[32]+0.4330127018922193*coeff[6]*fEdge[32]-0.21875*coeff[4]*fEdge[32]-0.43301270189221935*coeff[7]*fSkin[29]-0.7187500000000001*coeff[5]*fSkin[29]+0.43301270189221935*coeff[7]*fEdge[29]-0.21875*coeff[5]*fEdge[29]-0.375*coeff[6]*fSkin[28]-0.7036456405748563*coeff[4]*fSkin[28]-0.375*coeff[6]*fEdge[28]+0.2706329386826371*coeff[4]*fEdge[28]-0.375*coeff[7]*fSkin[22]-0.7036456405748562*coeff[5]*fSkin[22]-0.375*coeff[7]*fEdge[22]+0.27063293868263694*coeff[5]*fEdge[22]; 
  edgeSurf_incr[33] = -(0.27063293868263705*coeff[4]*fSkin[35])-0.27063293868263705*coeff[4]*fEdge[35]-0.27063293868263705*coeff[5]*fSkin[34]-0.27063293868263705*coeff[5]*fEdge[34]-0.28125*coeff[4]*fSkin[33]+0.28125*coeff[4]*fEdge[33]-0.28125*coeff[5]*fSkin[30]+0.28125*coeff[5]*fEdge[30]; 
  edgeSurf_incr[34] = -(0.4330127018922193*coeff[7]*fSkin[35])-0.71875*coeff[5]*fSkin[35]+0.4330127018922193*coeff[7]*fEdge[35]-0.21875*coeff[5]*fEdge[35]-0.4330127018922193*coeff[6]*fSkin[34]-0.71875*coeff[4]*fSkin[34]+0.4330127018922193*coeff[6]*fEdge[34]-0.21875*coeff[4]*fEdge[34]-0.375*coeff[7]*fSkin[33]-0.7036456405748562*coeff[5]*fSkin[33]-0.375*coeff[7]*fEdge[33]+0.27063293868263694*coeff[5]*fEdge[33]-0.375*coeff[6]*fSkin[30]-0.7036456405748562*coeff[4]*fSkin[30]-0.375*coeff[6]*fEdge[30]+0.27063293868263694*coeff[4]*fEdge[30]; 
  edgeSurf_incr[35] = -(0.4330127018922193*coeff[6]*fSkin[35])-0.71875*coeff[4]*fSkin[35]+0.4330127018922193*coeff[6]*fEdge[35]-0.21875*coeff[4]*fEdge[35]-0.4330127018922193*coeff[7]*fSkin[34]-0.71875*coeff[5]*fSkin[34]+0.4330127018922193*coeff[7]*fEdge[34]-0.21875*coeff[5]*fEdge[34]-0.375*coeff[6]*fSkin[33]-0.7036456405748562*coeff[4]*fSkin[33]-0.375*coeff[6]*fEdge[33]+0.27063293868263694*coeff[4]*fEdge[33]-0.375*coeff[7]*fSkin[30]-0.7036456405748562*coeff[5]*fSkin[30]-0.375*coeff[7]*fEdge[30]+0.27063293868263694*coeff[5]*fEdge[30]; 

  boundSurf_incr[2] = 0.8660254037844386*fSkin[5]*coeff[7]-0.75*fSkin[1]*coeff[7]+0.8660254037844386*fSkin[2]*coeff[6]-0.75*fSkin[0]*coeff[6]-0.5*coeff[5]*fSkin[5]+0.4330127018922193*fSkin[1]*coeff[5]-0.5*fSkin[2]*coeff[4]+0.4330127018922193*fSkin[0]*coeff[4]; 
  boundSurf_incr[5] = 0.8660254037844386*fSkin[2]*coeff[7]-0.75*fSkin[0]*coeff[7]+0.8660254037844386*fSkin[5]*coeff[6]-0.75*fSkin[1]*coeff[6]-0.5*coeff[4]*fSkin[5]-0.5*fSkin[2]*coeff[5]+0.4330127018922193*fSkin[0]*coeff[5]+0.4330127018922193*fSkin[1]*coeff[4]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[7]*fSkin[13]-0.5*coeff[5]*fSkin[13]+0.8660254037844386*coeff[6]*fSkin[7]-0.5*coeff[4]*fSkin[7]-0.75*fSkin[6]*coeff[7]+0.4330127018922193*coeff[5]*fSkin[6]-0.75*fSkin[3]*coeff[6]+0.4330127018922193*fSkin[3]*coeff[4]; 
  boundSurf_incr[9] = 0.8660254037844386*coeff[7]*fSkin[14]-0.5*coeff[5]*fSkin[14]+0.8660254037844386*coeff[6]*fSkin[9]-0.5*coeff[4]*fSkin[9]-0.75*coeff[7]*fSkin[8]+0.4330127018922193*coeff[5]*fSkin[8]-0.75*fSkin[4]*coeff[6]+0.4330127018922193*coeff[4]*fSkin[4]; 
  boundSurf_incr[13] = 0.8660254037844386*coeff[6]*fSkin[13]-0.5*coeff[4]*fSkin[13]+0.8660254037844386*coeff[7]*fSkin[7]-0.5*coeff[5]*fSkin[7]-0.75*fSkin[3]*coeff[7]-0.75*coeff[6]*fSkin[6]+0.4330127018922193*coeff[4]*fSkin[6]+0.4330127018922193*fSkin[3]*coeff[5]; 
  boundSurf_incr[14] = 0.8660254037844386*coeff[6]*fSkin[14]-0.5*coeff[4]*fSkin[14]+0.8660254037844386*coeff[7]*fSkin[9]-0.5*coeff[5]*fSkin[9]-0.75*coeff[6]*fSkin[8]+0.4330127018922193*coeff[4]*fSkin[8]-0.75*fSkin[4]*coeff[7]+0.4330127018922193*fSkin[4]*coeff[5]; 
  boundSurf_incr[16] = 0.8660254037844386*coeff[7]*fSkin[23]-0.5*coeff[5]*fSkin[23]+0.8660254037844386*coeff[6]*fSkin[16]-0.5*coeff[4]*fSkin[16]-0.75*coeff[7]*fSkin[15]+0.4330127018922193*coeff[5]*fSkin[15]-0.75*coeff[6]*fSkin[10]+0.4330127018922193*coeff[4]*fSkin[10]; 
  boundSurf_incr[18] = 0.8660254037844387*coeff[7]*fSkin[24]-0.5000000000000001*coeff[5]*fSkin[24]+0.8660254037844386*coeff[6]*fSkin[18]-0.5*coeff[4]*fSkin[18]-0.75*coeff[7]*fSkin[17]+0.4330127018922193*coeff[5]*fSkin[17]-0.75*coeff[6]*fSkin[11]+0.43301270189221935*coeff[4]*fSkin[11]; 
  boundSurf_incr[21] = 0.8660254037844387*coeff[7]*fSkin[27]-0.5000000000000001*coeff[5]*fSkin[27]+0.8660254037844386*coeff[6]*fSkin[21]-0.5*coeff[4]*fSkin[21]-0.75*coeff[7]*fSkin[20]+0.4330127018922193*coeff[5]*fSkin[20]-0.75*coeff[6]*fSkin[12]+0.43301270189221935*coeff[4]*fSkin[12]; 
  boundSurf_incr[23] = 0.8660254037844386*coeff[6]*fSkin[23]-0.5*coeff[4]*fSkin[23]+0.8660254037844386*coeff[7]*fSkin[16]-0.5*coeff[5]*fSkin[16]-0.75*coeff[6]*fSkin[15]+0.4330127018922193*coeff[4]*fSkin[15]-0.75*coeff[7]*fSkin[10]+0.4330127018922193*coeff[5]*fSkin[10]; 
  boundSurf_incr[24] = 0.8660254037844386*coeff[6]*fSkin[24]-0.5*coeff[4]*fSkin[24]+0.8660254037844387*coeff[7]*fSkin[18]-0.5000000000000001*coeff[5]*fSkin[18]-0.75*coeff[6]*fSkin[17]+0.43301270189221935*coeff[4]*fSkin[17]-0.75*coeff[7]*fSkin[11]+0.4330127018922193*coeff[5]*fSkin[11]; 
  boundSurf_incr[26] = 0.8660254037844387*coeff[7]*fSkin[31]-0.5000000000000001*coeff[5]*fSkin[31]+0.8660254037844386*coeff[6]*fSkin[26]-0.5*coeff[4]*fSkin[26]-0.75*coeff[7]*fSkin[25]+0.4330127018922193*coeff[5]*fSkin[25]-0.75*coeff[6]*fSkin[19]+0.43301270189221935*coeff[4]*fSkin[19]; 
  boundSurf_incr[27] = 0.8660254037844386*coeff[6]*fSkin[27]-0.5*coeff[4]*fSkin[27]+0.8660254037844387*coeff[7]*fSkin[21]-0.5000000000000001*coeff[5]*fSkin[21]-0.75*coeff[6]*fSkin[20]+0.43301270189221935*coeff[4]*fSkin[20]-0.75*coeff[7]*fSkin[12]+0.4330127018922193*coeff[5]*fSkin[12]; 
  boundSurf_incr[29] = 0.8660254037844387*coeff[7]*fSkin[32]-0.5000000000000001*coeff[5]*fSkin[32]+0.8660254037844386*coeff[6]*fSkin[29]-0.5*coeff[4]*fSkin[29]-0.75*coeff[7]*fSkin[28]+0.4330127018922193*coeff[5]*fSkin[28]-0.75*coeff[6]*fSkin[22]+0.43301270189221935*coeff[4]*fSkin[22]; 
  boundSurf_incr[31] = 0.8660254037844386*coeff[6]*fSkin[31]-0.5*coeff[4]*fSkin[31]+0.8660254037844387*coeff[7]*fSkin[26]-0.5000000000000001*coeff[5]*fSkin[26]-0.75*coeff[6]*fSkin[25]+0.43301270189221935*coeff[4]*fSkin[25]-0.75*coeff[7]*fSkin[19]+0.4330127018922193*coeff[5]*fSkin[19]; 
  boundSurf_incr[32] = 0.8660254037844386*coeff[6]*fSkin[32]-0.5*coeff[4]*fSkin[32]+0.8660254037844387*coeff[7]*fSkin[29]-0.5000000000000001*coeff[5]*fSkin[29]-0.75*coeff[6]*fSkin[28]+0.43301270189221935*coeff[4]*fSkin[28]-0.75*coeff[7]*fSkin[22]+0.4330127018922193*coeff[5]*fSkin[22]; 
  boundSurf_incr[34] = 0.8660254037844386*coeff[7]*fSkin[35]-0.5*coeff[5]*fSkin[35]+0.8660254037844386*coeff[6]*fSkin[34]-0.5*coeff[4]*fSkin[34]-0.75*coeff[7]*fSkin[33]+0.4330127018922193*coeff[5]*fSkin[33]-0.75*coeff[6]*fSkin[30]+0.4330127018922193*coeff[4]*fSkin[30]; 
  boundSurf_incr[35] = 0.8660254037844386*coeff[6]*fSkin[35]-0.5*coeff[4]*fSkin[35]+0.8660254037844386*coeff[7]*fSkin[34]-0.5*coeff[5]*fSkin[34]-0.75*coeff[6]*fSkin[33]+0.4330127018922193*coeff[4]*fSkin[33]-0.75*coeff[7]*fSkin[30]+0.4330127018922193*coeff[5]*fSkin[30]; 

  } else { 

  edgeSurf_incr[0] = 0.27063293868263705*coeff[5]*fSkin[5]+0.27063293868263705*coeff[5]*fEdge[5]-0.28125*fSkin[1]*coeff[5]+0.28125*fEdge[1]*coeff[5]+0.27063293868263705*fSkin[2]*coeff[4]+0.27063293868263705*fEdge[2]*coeff[4]-0.28125*fSkin[0]*coeff[4]+0.28125*fEdge[0]*coeff[4]; 
  edgeSurf_incr[1] = 0.27063293868263705*coeff[4]*fSkin[5]+0.27063293868263705*coeff[4]*fEdge[5]+0.27063293868263705*fSkin[2]*coeff[5]+0.27063293868263705*fEdge[2]*coeff[5]-0.28125*fSkin[0]*coeff[5]+0.28125*fEdge[0]*coeff[5]-0.28125*fSkin[1]*coeff[4]+0.28125*fEdge[1]*coeff[4]; 
  edgeSurf_incr[2] = 0.4330127018922193*fSkin[5]*coeff[7]-0.4330127018922193*fEdge[5]*coeff[7]-0.375*fSkin[1]*coeff[7]-0.375*fEdge[1]*coeff[7]+0.4330127018922193*fSkin[2]*coeff[6]-0.4330127018922193*fEdge[2]*coeff[6]-0.375*fSkin[0]*coeff[6]-0.375*fEdge[0]*coeff[6]-0.71875*coeff[5]*fSkin[5]-0.21875*coeff[5]*fEdge[5]+0.7036456405748562*fSkin[1]*coeff[5]-0.27063293868263694*fEdge[1]*coeff[5]-0.71875*fSkin[2]*coeff[4]-0.21875*fEdge[2]*coeff[4]+0.7036456405748562*fSkin[0]*coeff[4]-0.27063293868263694*fEdge[0]*coeff[4]; 
  edgeSurf_incr[3] = 0.27063293868263705*coeff[5]*fSkin[13]+0.27063293868263705*coeff[5]*fEdge[13]+0.27063293868263705*coeff[4]*fSkin[7]+0.27063293868263705*coeff[4]*fEdge[7]-0.28125*coeff[5]*fSkin[6]+0.28125*coeff[5]*fEdge[6]-0.28125*fSkin[3]*coeff[4]+0.28125*fEdge[3]*coeff[4]; 
  edgeSurf_incr[4] = 0.27063293868263705*coeff[5]*fSkin[14]+0.27063293868263705*coeff[5]*fEdge[14]+0.27063293868263705*coeff[4]*fSkin[9]+0.27063293868263705*coeff[4]*fEdge[9]-0.28125*coeff[5]*fSkin[8]+0.28125*coeff[5]*fEdge[8]-0.28125*coeff[4]*fSkin[4]+0.28125*coeff[4]*fEdge[4]; 
  edgeSurf_incr[5] = 0.4330127018922193*fSkin[2]*coeff[7]-0.4330127018922193*fEdge[2]*coeff[7]-0.375*fSkin[0]*coeff[7]-0.375*fEdge[0]*coeff[7]+0.4330127018922193*fSkin[5]*coeff[6]-0.4330127018922193*fEdge[5]*coeff[6]-0.375*fSkin[1]*coeff[6]-0.375*fEdge[1]*coeff[6]-0.71875*coeff[4]*fSkin[5]-0.21875*coeff[4]*fEdge[5]-0.71875*fSkin[2]*coeff[5]-0.21875*fEdge[2]*coeff[5]+0.7036456405748562*fSkin[0]*coeff[5]-0.27063293868263694*fEdge[0]*coeff[5]+0.7036456405748562*fSkin[1]*coeff[4]-0.27063293868263694*fEdge[1]*coeff[4]; 
  edgeSurf_incr[6] = 0.27063293868263705*coeff[4]*fSkin[13]+0.27063293868263705*coeff[4]*fEdge[13]+0.27063293868263705*coeff[5]*fSkin[7]+0.27063293868263705*coeff[5]*fEdge[7]-0.28125*coeff[4]*fSkin[6]+0.28125*coeff[4]*fEdge[6]-0.28125*fSkin[3]*coeff[5]+0.28125*fEdge[3]*coeff[5]; 
  edgeSurf_incr[7] = 0.4330127018922193*coeff[7]*fSkin[13]-0.71875*coeff[5]*fSkin[13]-0.4330127018922193*coeff[7]*fEdge[13]-0.21875*coeff[5]*fEdge[13]+0.4330127018922193*coeff[6]*fSkin[7]-0.71875*coeff[4]*fSkin[7]-0.4330127018922193*coeff[6]*fEdge[7]-0.21875*coeff[4]*fEdge[7]-0.375*fSkin[6]*coeff[7]-0.375*fEdge[6]*coeff[7]+0.7036456405748562*coeff[5]*fSkin[6]-0.27063293868263694*coeff[5]*fEdge[6]-0.375*fSkin[3]*coeff[6]-0.375*fEdge[3]*coeff[6]+0.7036456405748562*fSkin[3]*coeff[4]-0.27063293868263694*fEdge[3]*coeff[4]; 
  edgeSurf_incr[8] = 0.27063293868263705*coeff[4]*fSkin[14]+0.27063293868263705*coeff[4]*fEdge[14]+0.27063293868263705*coeff[5]*fSkin[9]+0.27063293868263705*coeff[5]*fEdge[9]-0.28125*coeff[4]*fSkin[8]+0.28125*coeff[4]*fEdge[8]-0.28125*fSkin[4]*coeff[5]+0.28125*fEdge[4]*coeff[5]; 
  edgeSurf_incr[9] = 0.4330127018922193*coeff[7]*fSkin[14]-0.71875*coeff[5]*fSkin[14]-0.4330127018922193*coeff[7]*fEdge[14]-0.21875*coeff[5]*fEdge[14]+0.4330127018922193*coeff[6]*fSkin[9]-0.71875*coeff[4]*fSkin[9]-0.4330127018922193*coeff[6]*fEdge[9]-0.21875*coeff[4]*fEdge[9]-0.375*coeff[7]*fSkin[8]+0.7036456405748562*coeff[5]*fSkin[8]-0.375*coeff[7]*fEdge[8]-0.27063293868263694*coeff[5]*fEdge[8]-0.375*fSkin[4]*coeff[6]-0.375*fEdge[4]*coeff[6]+0.7036456405748562*coeff[4]*fSkin[4]-0.27063293868263694*coeff[4]*fEdge[4]; 
  edgeSurf_incr[10] = 0.27063293868263705*coeff[5]*fSkin[23]+0.27063293868263705*coeff[5]*fEdge[23]+0.27063293868263705*coeff[4]*fSkin[16]+0.27063293868263705*coeff[4]*fEdge[16]-0.28125*coeff[5]*fSkin[15]+0.28125*coeff[5]*fEdge[15]-0.28125*coeff[4]*fSkin[10]+0.28125*coeff[4]*fEdge[10]; 
  edgeSurf_incr[11] = 0.27063293868263705*coeff[5]*fSkin[24]+0.27063293868263705*coeff[5]*fEdge[24]+0.2706329386826371*coeff[4]*fSkin[18]+0.2706329386826371*coeff[4]*fEdge[18]-0.28124999999999994*coeff[5]*fSkin[17]+0.28124999999999994*coeff[5]*fEdge[17]-0.28125*coeff[4]*fSkin[11]+0.28125*coeff[4]*fEdge[11]; 
  edgeSurf_incr[12] = 0.27063293868263705*coeff[5]*fSkin[27]+0.27063293868263705*coeff[5]*fEdge[27]+0.2706329386826371*coeff[4]*fSkin[21]+0.2706329386826371*coeff[4]*fEdge[21]-0.28124999999999994*coeff[5]*fSkin[20]+0.28124999999999994*coeff[5]*fEdge[20]-0.28125*coeff[4]*fSkin[12]+0.28125*coeff[4]*fEdge[12]; 
  edgeSurf_incr[13] = 0.4330127018922193*coeff[6]*fSkin[13]-0.71875*coeff[4]*fSkin[13]-0.4330127018922193*coeff[6]*fEdge[13]-0.21875*coeff[4]*fEdge[13]+0.4330127018922193*coeff[7]*fSkin[7]-0.71875*coeff[5]*fSkin[7]-0.4330127018922193*coeff[7]*fEdge[7]-0.21875*coeff[5]*fEdge[7]-0.375*fSkin[3]*coeff[7]-0.375*fEdge[3]*coeff[7]-0.375*coeff[6]*fSkin[6]+0.7036456405748562*coeff[4]*fSkin[6]-0.375*coeff[6]*fEdge[6]-0.27063293868263694*coeff[4]*fEdge[6]+0.7036456405748562*fSkin[3]*coeff[5]-0.27063293868263694*fEdge[3]*coeff[5]; 
  edgeSurf_incr[14] = 0.4330127018922193*coeff[6]*fSkin[14]-0.71875*coeff[4]*fSkin[14]-0.4330127018922193*coeff[6]*fEdge[14]-0.21875*coeff[4]*fEdge[14]+0.4330127018922193*coeff[7]*fSkin[9]-0.71875*coeff[5]*fSkin[9]-0.4330127018922193*coeff[7]*fEdge[9]-0.21875*coeff[5]*fEdge[9]-0.375*coeff[6]*fSkin[8]+0.7036456405748562*coeff[4]*fSkin[8]-0.375*coeff[6]*fEdge[8]-0.27063293868263694*coeff[4]*fEdge[8]-0.375*fSkin[4]*coeff[7]-0.375*fEdge[4]*coeff[7]+0.7036456405748562*fSkin[4]*coeff[5]-0.27063293868263694*fEdge[4]*coeff[5]; 
  edgeSurf_incr[15] = 0.27063293868263705*coeff[4]*fSkin[23]+0.27063293868263705*coeff[4]*fEdge[23]+0.27063293868263705*coeff[5]*fSkin[16]+0.27063293868263705*coeff[5]*fEdge[16]-0.28125*coeff[4]*fSkin[15]+0.28125*coeff[4]*fEdge[15]-0.28125*coeff[5]*fSkin[10]+0.28125*coeff[5]*fEdge[10]; 
  edgeSurf_incr[16] = 0.4330127018922193*coeff[7]*fSkin[23]-0.71875*coeff[5]*fSkin[23]-0.4330127018922193*coeff[7]*fEdge[23]-0.21875*coeff[5]*fEdge[23]+0.4330127018922193*coeff[6]*fSkin[16]-0.71875*coeff[4]*fSkin[16]-0.4330127018922193*coeff[6]*fEdge[16]-0.21875*coeff[4]*fEdge[16]-0.375*coeff[7]*fSkin[15]+0.7036456405748562*coeff[5]*fSkin[15]-0.375*coeff[7]*fEdge[15]-0.27063293868263694*coeff[5]*fEdge[15]-0.375*coeff[6]*fSkin[10]+0.7036456405748562*coeff[4]*fSkin[10]-0.375*coeff[6]*fEdge[10]-0.27063293868263694*coeff[4]*fEdge[10]; 
  edgeSurf_incr[17] = 0.2706329386826371*coeff[4]*fSkin[24]+0.2706329386826371*coeff[4]*fEdge[24]+0.27063293868263705*coeff[5]*fSkin[18]+0.27063293868263705*coeff[5]*fEdge[18]-0.28125*coeff[4]*fSkin[17]+0.28125*coeff[4]*fEdge[17]-0.28124999999999994*coeff[5]*fSkin[11]+0.28124999999999994*coeff[5]*fEdge[11]; 
  edgeSurf_incr[18] = 0.43301270189221935*coeff[7]*fSkin[24]-0.7187500000000001*coeff[5]*fSkin[24]-0.43301270189221935*coeff[7]*fEdge[24]-0.21875*coeff[5]*fEdge[24]+0.4330127018922193*coeff[6]*fSkin[18]-0.71875*coeff[4]*fSkin[18]-0.4330127018922193*coeff[6]*fEdge[18]-0.21875*coeff[4]*fEdge[18]-0.375*coeff[7]*fSkin[17]+0.7036456405748562*coeff[5]*fSkin[17]-0.375*coeff[7]*fEdge[17]-0.27063293868263694*coeff[5]*fEdge[17]-0.375*coeff[6]*fSkin[11]+0.7036456405748563*coeff[4]*fSkin[11]-0.375*coeff[6]*fEdge[11]-0.2706329386826371*coeff[4]*fEdge[11]; 
  edgeSurf_incr[19] = 0.27063293868263705*coeff[5]*fSkin[31]+0.27063293868263705*coeff[5]*fEdge[31]+0.2706329386826371*coeff[4]*fSkin[26]+0.2706329386826371*coeff[4]*fEdge[26]-0.28124999999999994*coeff[5]*fSkin[25]+0.28124999999999994*coeff[5]*fEdge[25]-0.28125*coeff[4]*fSkin[19]+0.28125*coeff[4]*fEdge[19]; 
  edgeSurf_incr[20] = 0.2706329386826371*coeff[4]*fSkin[27]+0.2706329386826371*coeff[4]*fEdge[27]+0.27063293868263705*coeff[5]*fSkin[21]+0.27063293868263705*coeff[5]*fEdge[21]-0.28125*coeff[4]*fSkin[20]+0.28125*coeff[4]*fEdge[20]-0.28124999999999994*coeff[5]*fSkin[12]+0.28124999999999994*coeff[5]*fEdge[12]; 
  edgeSurf_incr[21] = 0.43301270189221935*coeff[7]*fSkin[27]-0.7187500000000001*coeff[5]*fSkin[27]-0.43301270189221935*coeff[7]*fEdge[27]-0.21875*coeff[5]*fEdge[27]+0.4330127018922193*coeff[6]*fSkin[21]-0.71875*coeff[4]*fSkin[21]-0.4330127018922193*coeff[6]*fEdge[21]-0.21875*coeff[4]*fEdge[21]-0.375*coeff[7]*fSkin[20]+0.7036456405748562*coeff[5]*fSkin[20]-0.375*coeff[7]*fEdge[20]-0.27063293868263694*coeff[5]*fEdge[20]-0.375*coeff[6]*fSkin[12]+0.7036456405748563*coeff[4]*fSkin[12]-0.375*coeff[6]*fEdge[12]-0.2706329386826371*coeff[4]*fEdge[12]; 
  edgeSurf_incr[22] = 0.27063293868263705*coeff[5]*fSkin[32]+0.27063293868263705*coeff[5]*fEdge[32]+0.2706329386826371*coeff[4]*fSkin[29]+0.2706329386826371*coeff[4]*fEdge[29]-0.28124999999999994*coeff[5]*fSkin[28]+0.28124999999999994*coeff[5]*fEdge[28]-0.28125*coeff[4]*fSkin[22]+0.28125*coeff[4]*fEdge[22]; 
  edgeSurf_incr[23] = 0.4330127018922193*coeff[6]*fSkin[23]-0.71875*coeff[4]*fSkin[23]-0.4330127018922193*coeff[6]*fEdge[23]-0.21875*coeff[4]*fEdge[23]+0.4330127018922193*coeff[7]*fSkin[16]-0.71875*coeff[5]*fSkin[16]-0.4330127018922193*coeff[7]*fEdge[16]-0.21875*coeff[5]*fEdge[16]-0.375*coeff[6]*fSkin[15]+0.7036456405748562*coeff[4]*fSkin[15]-0.375*coeff[6]*fEdge[15]-0.27063293868263694*coeff[4]*fEdge[15]-0.375*coeff[7]*fSkin[10]+0.7036456405748562*coeff[5]*fSkin[10]-0.375*coeff[7]*fEdge[10]-0.27063293868263694*coeff[5]*fEdge[10]; 
  edgeSurf_incr[24] = 0.4330127018922193*coeff[6]*fSkin[24]-0.71875*coeff[4]*fSkin[24]-0.4330127018922193*coeff[6]*fEdge[24]-0.21875*coeff[4]*fEdge[24]+0.43301270189221935*coeff[7]*fSkin[18]-0.7187500000000001*coeff[5]*fSkin[18]-0.43301270189221935*coeff[7]*fEdge[18]-0.21875*coeff[5]*fEdge[18]-0.375*coeff[6]*fSkin[17]+0.7036456405748563*coeff[4]*fSkin[17]-0.375*coeff[6]*fEdge[17]-0.2706329386826371*coeff[4]*fEdge[17]-0.375*coeff[7]*fSkin[11]+0.7036456405748562*coeff[5]*fSkin[11]-0.375*coeff[7]*fEdge[11]-0.27063293868263694*coeff[5]*fEdge[11]; 
  edgeSurf_incr[25] = 0.2706329386826371*coeff[4]*fSkin[31]+0.2706329386826371*coeff[4]*fEdge[31]+0.27063293868263705*coeff[5]*fSkin[26]+0.27063293868263705*coeff[5]*fEdge[26]-0.28125*coeff[4]*fSkin[25]+0.28125*coeff[4]*fEdge[25]-0.28124999999999994*coeff[5]*fSkin[19]+0.28124999999999994*coeff[5]*fEdge[19]; 
  edgeSurf_incr[26] = 0.43301270189221935*coeff[7]*fSkin[31]-0.7187500000000001*coeff[5]*fSkin[31]-0.43301270189221935*coeff[7]*fEdge[31]-0.21875*coeff[5]*fEdge[31]+0.4330127018922193*coeff[6]*fSkin[26]-0.71875*coeff[4]*fSkin[26]-0.4330127018922193*coeff[6]*fEdge[26]-0.21875*coeff[4]*fEdge[26]-0.375*coeff[7]*fSkin[25]+0.7036456405748562*coeff[5]*fSkin[25]-0.375*coeff[7]*fEdge[25]-0.27063293868263694*coeff[5]*fEdge[25]-0.375*coeff[6]*fSkin[19]+0.7036456405748563*coeff[4]*fSkin[19]-0.375*coeff[6]*fEdge[19]-0.2706329386826371*coeff[4]*fEdge[19]; 
  edgeSurf_incr[27] = 0.4330127018922193*coeff[6]*fSkin[27]-0.71875*coeff[4]*fSkin[27]-0.4330127018922193*coeff[6]*fEdge[27]-0.21875*coeff[4]*fEdge[27]+0.43301270189221935*coeff[7]*fSkin[21]-0.7187500000000001*coeff[5]*fSkin[21]-0.43301270189221935*coeff[7]*fEdge[21]-0.21875*coeff[5]*fEdge[21]-0.375*coeff[6]*fSkin[20]+0.7036456405748563*coeff[4]*fSkin[20]-0.375*coeff[6]*fEdge[20]-0.2706329386826371*coeff[4]*fEdge[20]-0.375*coeff[7]*fSkin[12]+0.7036456405748562*coeff[5]*fSkin[12]-0.375*coeff[7]*fEdge[12]-0.27063293868263694*coeff[5]*fEdge[12]; 
  edgeSurf_incr[28] = 0.2706329386826371*coeff[4]*fSkin[32]+0.2706329386826371*coeff[4]*fEdge[32]+0.27063293868263705*coeff[5]*fSkin[29]+0.27063293868263705*coeff[5]*fEdge[29]-0.28125*coeff[4]*fSkin[28]+0.28125*coeff[4]*fEdge[28]-0.28124999999999994*coeff[5]*fSkin[22]+0.28124999999999994*coeff[5]*fEdge[22]; 
  edgeSurf_incr[29] = 0.43301270189221935*coeff[7]*fSkin[32]-0.7187500000000001*coeff[5]*fSkin[32]-0.43301270189221935*coeff[7]*fEdge[32]-0.21875*coeff[5]*fEdge[32]+0.4330127018922193*coeff[6]*fSkin[29]-0.71875*coeff[4]*fSkin[29]-0.4330127018922193*coeff[6]*fEdge[29]-0.21875*coeff[4]*fEdge[29]-0.375*coeff[7]*fSkin[28]+0.7036456405748562*coeff[5]*fSkin[28]-0.375*coeff[7]*fEdge[28]-0.27063293868263694*coeff[5]*fEdge[28]-0.375*coeff[6]*fSkin[22]+0.7036456405748563*coeff[4]*fSkin[22]-0.375*coeff[6]*fEdge[22]-0.2706329386826371*coeff[4]*fEdge[22]; 
  edgeSurf_incr[30] = 0.27063293868263705*coeff[5]*fSkin[35]+0.27063293868263705*coeff[5]*fEdge[35]+0.27063293868263705*coeff[4]*fSkin[34]+0.27063293868263705*coeff[4]*fEdge[34]-0.28125*coeff[5]*fSkin[33]+0.28125*coeff[5]*fEdge[33]-0.28125*coeff[4]*fSkin[30]+0.28125*coeff[4]*fEdge[30]; 
  edgeSurf_incr[31] = 0.4330127018922193*coeff[6]*fSkin[31]-0.71875*coeff[4]*fSkin[31]-0.4330127018922193*coeff[6]*fEdge[31]-0.21875*coeff[4]*fEdge[31]+0.43301270189221935*coeff[7]*fSkin[26]-0.7187500000000001*coeff[5]*fSkin[26]-0.43301270189221935*coeff[7]*fEdge[26]-0.21875*coeff[5]*fEdge[26]-0.375*coeff[6]*fSkin[25]+0.7036456405748563*coeff[4]*fSkin[25]-0.375*coeff[6]*fEdge[25]-0.2706329386826371*coeff[4]*fEdge[25]-0.375*coeff[7]*fSkin[19]+0.7036456405748562*coeff[5]*fSkin[19]-0.375*coeff[7]*fEdge[19]-0.27063293868263694*coeff[5]*fEdge[19]; 
  edgeSurf_incr[32] = 0.4330127018922193*coeff[6]*fSkin[32]-0.71875*coeff[4]*fSkin[32]-0.4330127018922193*coeff[6]*fEdge[32]-0.21875*coeff[4]*fEdge[32]+0.43301270189221935*coeff[7]*fSkin[29]-0.7187500000000001*coeff[5]*fSkin[29]-0.43301270189221935*coeff[7]*fEdge[29]-0.21875*coeff[5]*fEdge[29]-0.375*coeff[6]*fSkin[28]+0.7036456405748563*coeff[4]*fSkin[28]-0.375*coeff[6]*fEdge[28]-0.2706329386826371*coeff[4]*fEdge[28]-0.375*coeff[7]*fSkin[22]+0.7036456405748562*coeff[5]*fSkin[22]-0.375*coeff[7]*fEdge[22]-0.27063293868263694*coeff[5]*fEdge[22]; 
  edgeSurf_incr[33] = 0.27063293868263705*coeff[4]*fSkin[35]+0.27063293868263705*coeff[4]*fEdge[35]+0.27063293868263705*coeff[5]*fSkin[34]+0.27063293868263705*coeff[5]*fEdge[34]-0.28125*coeff[4]*fSkin[33]+0.28125*coeff[4]*fEdge[33]-0.28125*coeff[5]*fSkin[30]+0.28125*coeff[5]*fEdge[30]; 
  edgeSurf_incr[34] = 0.4330127018922193*coeff[7]*fSkin[35]-0.71875*coeff[5]*fSkin[35]-0.4330127018922193*coeff[7]*fEdge[35]-0.21875*coeff[5]*fEdge[35]+0.4330127018922193*coeff[6]*fSkin[34]-0.71875*coeff[4]*fSkin[34]-0.4330127018922193*coeff[6]*fEdge[34]-0.21875*coeff[4]*fEdge[34]-0.375*coeff[7]*fSkin[33]+0.7036456405748562*coeff[5]*fSkin[33]-0.375*coeff[7]*fEdge[33]-0.27063293868263694*coeff[5]*fEdge[33]-0.375*coeff[6]*fSkin[30]+0.7036456405748562*coeff[4]*fSkin[30]-0.375*coeff[6]*fEdge[30]-0.27063293868263694*coeff[4]*fEdge[30]; 
  edgeSurf_incr[35] = 0.4330127018922193*coeff[6]*fSkin[35]-0.71875*coeff[4]*fSkin[35]-0.4330127018922193*coeff[6]*fEdge[35]-0.21875*coeff[4]*fEdge[35]+0.4330127018922193*coeff[7]*fSkin[34]-0.71875*coeff[5]*fSkin[34]-0.4330127018922193*coeff[7]*fEdge[34]-0.21875*coeff[5]*fEdge[34]-0.375*coeff[6]*fSkin[33]+0.7036456405748562*coeff[4]*fSkin[33]-0.375*coeff[6]*fEdge[33]-0.27063293868263694*coeff[4]*fEdge[33]-0.375*coeff[7]*fSkin[30]+0.7036456405748562*coeff[5]*fSkin[30]-0.375*coeff[7]*fEdge[30]-0.27063293868263694*coeff[5]*fEdge[30]; 

  boundSurf_incr[2] = -(0.8660254037844386*fSkin[5]*coeff[7])-0.75*fSkin[1]*coeff[7]-0.8660254037844386*fSkin[2]*coeff[6]-0.75*fSkin[0]*coeff[6]-0.5*coeff[5]*fSkin[5]-0.4330127018922193*fSkin[1]*coeff[5]-0.5*fSkin[2]*coeff[4]-0.4330127018922193*fSkin[0]*coeff[4]; 
  boundSurf_incr[5] = -(0.8660254037844386*fSkin[2]*coeff[7])-0.75*fSkin[0]*coeff[7]-0.8660254037844386*fSkin[5]*coeff[6]-0.75*fSkin[1]*coeff[6]-0.5*coeff[4]*fSkin[5]-0.5*fSkin[2]*coeff[5]-0.4330127018922193*fSkin[0]*coeff[5]-0.4330127018922193*fSkin[1]*coeff[4]; 
  boundSurf_incr[7] = -(0.8660254037844386*coeff[7]*fSkin[13])-0.5*coeff[5]*fSkin[13]-0.8660254037844386*coeff[6]*fSkin[7]-0.5*coeff[4]*fSkin[7]-0.75*fSkin[6]*coeff[7]-0.4330127018922193*coeff[5]*fSkin[6]-0.75*fSkin[3]*coeff[6]-0.4330127018922193*fSkin[3]*coeff[4]; 
  boundSurf_incr[9] = -(0.8660254037844386*coeff[7]*fSkin[14])-0.5*coeff[5]*fSkin[14]-0.8660254037844386*coeff[6]*fSkin[9]-0.5*coeff[4]*fSkin[9]-0.75*coeff[7]*fSkin[8]-0.4330127018922193*coeff[5]*fSkin[8]-0.75*fSkin[4]*coeff[6]-0.4330127018922193*coeff[4]*fSkin[4]; 
  boundSurf_incr[13] = -(0.8660254037844386*coeff[6]*fSkin[13])-0.5*coeff[4]*fSkin[13]-0.8660254037844386*coeff[7]*fSkin[7]-0.5*coeff[5]*fSkin[7]-0.75*fSkin[3]*coeff[7]-0.75*coeff[6]*fSkin[6]-0.4330127018922193*coeff[4]*fSkin[6]-0.4330127018922193*fSkin[3]*coeff[5]; 
  boundSurf_incr[14] = -(0.8660254037844386*coeff[6]*fSkin[14])-0.5*coeff[4]*fSkin[14]-0.8660254037844386*coeff[7]*fSkin[9]-0.5*coeff[5]*fSkin[9]-0.75*coeff[6]*fSkin[8]-0.4330127018922193*coeff[4]*fSkin[8]-0.75*fSkin[4]*coeff[7]-0.4330127018922193*fSkin[4]*coeff[5]; 
  boundSurf_incr[16] = -(0.8660254037844386*coeff[7]*fSkin[23])-0.5*coeff[5]*fSkin[23]-0.8660254037844386*coeff[6]*fSkin[16]-0.5*coeff[4]*fSkin[16]-0.75*coeff[7]*fSkin[15]-0.4330127018922193*coeff[5]*fSkin[15]-0.75*coeff[6]*fSkin[10]-0.4330127018922193*coeff[4]*fSkin[10]; 
  boundSurf_incr[18] = -(0.8660254037844387*coeff[7]*fSkin[24])-0.5000000000000001*coeff[5]*fSkin[24]-0.8660254037844386*coeff[6]*fSkin[18]-0.5*coeff[4]*fSkin[18]-0.75*coeff[7]*fSkin[17]-0.4330127018922193*coeff[5]*fSkin[17]-0.75*coeff[6]*fSkin[11]-0.43301270189221935*coeff[4]*fSkin[11]; 
  boundSurf_incr[21] = -(0.8660254037844387*coeff[7]*fSkin[27])-0.5000000000000001*coeff[5]*fSkin[27]-0.8660254037844386*coeff[6]*fSkin[21]-0.5*coeff[4]*fSkin[21]-0.75*coeff[7]*fSkin[20]-0.4330127018922193*coeff[5]*fSkin[20]-0.75*coeff[6]*fSkin[12]-0.43301270189221935*coeff[4]*fSkin[12]; 
  boundSurf_incr[23] = -(0.8660254037844386*coeff[6]*fSkin[23])-0.5*coeff[4]*fSkin[23]-0.8660254037844386*coeff[7]*fSkin[16]-0.5*coeff[5]*fSkin[16]-0.75*coeff[6]*fSkin[15]-0.4330127018922193*coeff[4]*fSkin[15]-0.75*coeff[7]*fSkin[10]-0.4330127018922193*coeff[5]*fSkin[10]; 
  boundSurf_incr[24] = -(0.8660254037844386*coeff[6]*fSkin[24])-0.5*coeff[4]*fSkin[24]-0.8660254037844387*coeff[7]*fSkin[18]-0.5000000000000001*coeff[5]*fSkin[18]-0.75*coeff[6]*fSkin[17]-0.43301270189221935*coeff[4]*fSkin[17]-0.75*coeff[7]*fSkin[11]-0.4330127018922193*coeff[5]*fSkin[11]; 
  boundSurf_incr[26] = -(0.8660254037844387*coeff[7]*fSkin[31])-0.5000000000000001*coeff[5]*fSkin[31]-0.8660254037844386*coeff[6]*fSkin[26]-0.5*coeff[4]*fSkin[26]-0.75*coeff[7]*fSkin[25]-0.4330127018922193*coeff[5]*fSkin[25]-0.75*coeff[6]*fSkin[19]-0.43301270189221935*coeff[4]*fSkin[19]; 
  boundSurf_incr[27] = -(0.8660254037844386*coeff[6]*fSkin[27])-0.5*coeff[4]*fSkin[27]-0.8660254037844387*coeff[7]*fSkin[21]-0.5000000000000001*coeff[5]*fSkin[21]-0.75*coeff[6]*fSkin[20]-0.43301270189221935*coeff[4]*fSkin[20]-0.75*coeff[7]*fSkin[12]-0.4330127018922193*coeff[5]*fSkin[12]; 
  boundSurf_incr[29] = -(0.8660254037844387*coeff[7]*fSkin[32])-0.5000000000000001*coeff[5]*fSkin[32]-0.8660254037844386*coeff[6]*fSkin[29]-0.5*coeff[4]*fSkin[29]-0.75*coeff[7]*fSkin[28]-0.4330127018922193*coeff[5]*fSkin[28]-0.75*coeff[6]*fSkin[22]-0.43301270189221935*coeff[4]*fSkin[22]; 
  boundSurf_incr[31] = -(0.8660254037844386*coeff[6]*fSkin[31])-0.5*coeff[4]*fSkin[31]-0.8660254037844387*coeff[7]*fSkin[26]-0.5000000000000001*coeff[5]*fSkin[26]-0.75*coeff[6]*fSkin[25]-0.43301270189221935*coeff[4]*fSkin[25]-0.75*coeff[7]*fSkin[19]-0.4330127018922193*coeff[5]*fSkin[19]; 
  boundSurf_incr[32] = -(0.8660254037844386*coeff[6]*fSkin[32])-0.5*coeff[4]*fSkin[32]-0.8660254037844387*coeff[7]*fSkin[29]-0.5000000000000001*coeff[5]*fSkin[29]-0.75*coeff[6]*fSkin[28]-0.43301270189221935*coeff[4]*fSkin[28]-0.75*coeff[7]*fSkin[22]-0.4330127018922193*coeff[5]*fSkin[22]; 
  boundSurf_incr[34] = -(0.8660254037844386*coeff[7]*fSkin[35])-0.5*coeff[5]*fSkin[35]-0.8660254037844386*coeff[6]*fSkin[34]-0.5*coeff[4]*fSkin[34]-0.75*coeff[7]*fSkin[33]-0.4330127018922193*coeff[5]*fSkin[33]-0.75*coeff[6]*fSkin[30]-0.4330127018922193*coeff[4]*fSkin[30]; 
  boundSurf_incr[35] = -(0.8660254037844386*coeff[6]*fSkin[35])-0.5*coeff[4]*fSkin[35]-0.8660254037844386*coeff[7]*fSkin[34]-0.5*coeff[5]*fSkin[34]-0.75*coeff[6]*fSkin[33]-0.4330127018922193*coeff[4]*fSkin[33]-0.75*coeff[7]*fSkin[30]-0.4330127018922193*coeff[5]*fSkin[30]; 

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

  return 0.;
}

