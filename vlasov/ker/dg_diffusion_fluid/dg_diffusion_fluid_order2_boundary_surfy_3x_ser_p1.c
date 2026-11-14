#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.5412658773652741*coeff[1]*fSkin[2])-0.5412658773652741*coeff[1]*fEdge[2]-0.5625*fSkin[0]*coeff[1]+0.5625*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(0.5412658773652741*coeff[1]*fSkin[4])-0.5412658773652741*coeff[1]*fEdge[4]-0.5625*coeff[1]*fSkin[1]+0.5625*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(1.4375*coeff[1]*fSkin[2])-0.4375*coeff[1]*fEdge[2]-1.4072912811497125*fSkin[0]*coeff[1]+0.5412658773652739*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(0.5412658773652741*coeff[1]*fSkin[6])-0.5412658773652741*coeff[1]*fEdge[6]-0.5625*coeff[1]*fSkin[3]+0.5625*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = -(1.4375*coeff[1]*fSkin[4])-0.4375*coeff[1]*fEdge[4]-1.4072912811497125*coeff[1]*fSkin[1]+0.5412658773652739*coeff[1]*fEdge[1]; 
  edgeSurf_incr[5] = -(0.5412658773652741*coeff[1]*fSkin[7])-0.5412658773652741*coeff[1]*fEdge[7]-0.5625*coeff[1]*fSkin[5]+0.5625*coeff[1]*fEdge[5]; 
  edgeSurf_incr[6] = -(1.4375*coeff[1]*fSkin[6])-0.4375*coeff[1]*fEdge[6]-1.4072912811497125*coeff[1]*fSkin[3]+0.5412658773652739*coeff[1]*fEdge[3]; 
  edgeSurf_incr[7] = -(1.4375*coeff[1]*fSkin[7])-0.4375*coeff[1]*fEdge[7]-1.4072912811497125*coeff[1]*fSkin[5]+0.5412658773652739*coeff[1]*fEdge[5]; 

  boundSurf_incr[2] = 0.8660254037844386*fSkin[0]*coeff[1]-1.0*coeff[1]*fSkin[2]; 
  boundSurf_incr[4] = 0.8660254037844386*coeff[1]*fSkin[1]-1.0*coeff[1]*fSkin[4]; 
  boundSurf_incr[6] = 0.8660254037844386*coeff[1]*fSkin[3]-1.0*coeff[1]*fSkin[6]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[1]*fSkin[5]-1.0*coeff[1]*fSkin[7]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[1]*fSkin[2]+0.5412658773652741*coeff[1]*fEdge[2]-0.5625*fSkin[0]*coeff[1]+0.5625*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 0.5412658773652741*coeff[1]*fSkin[4]+0.5412658773652741*coeff[1]*fEdge[4]-0.5625*coeff[1]*fSkin[1]+0.5625*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(1.4375*coeff[1]*fSkin[2])-0.4375*coeff[1]*fEdge[2]+1.4072912811497125*fSkin[0]*coeff[1]-0.5412658773652739*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 0.5412658773652741*coeff[1]*fSkin[6]+0.5412658773652741*coeff[1]*fEdge[6]-0.5625*coeff[1]*fSkin[3]+0.5625*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = -(1.4375*coeff[1]*fSkin[4])-0.4375*coeff[1]*fEdge[4]+1.4072912811497125*coeff[1]*fSkin[1]-0.5412658773652739*coeff[1]*fEdge[1]; 
  edgeSurf_incr[5] = 0.5412658773652741*coeff[1]*fSkin[7]+0.5412658773652741*coeff[1]*fEdge[7]-0.5625*coeff[1]*fSkin[5]+0.5625*coeff[1]*fEdge[5]; 
  edgeSurf_incr[6] = -(1.4375*coeff[1]*fSkin[6])-0.4375*coeff[1]*fEdge[6]+1.4072912811497125*coeff[1]*fSkin[3]-0.5412658773652739*coeff[1]*fEdge[3]; 
  edgeSurf_incr[7] = -(1.4375*coeff[1]*fSkin[7])-0.4375*coeff[1]*fEdge[7]+1.4072912811497125*coeff[1]*fSkin[5]-0.5412658773652739*coeff[1]*fEdge[5]; 

  boundSurf_incr[2] = -(1.0*coeff[1]*fSkin[2])-0.8660254037844386*fSkin[0]*coeff[1]; 
  boundSurf_incr[4] = -(1.0*coeff[1]*fSkin[4])-0.8660254037844386*coeff[1]*fSkin[1]; 
  boundSurf_incr[6] = -(1.0*coeff[1]*fSkin[6])-0.8660254037844386*coeff[1]*fSkin[3]; 
  boundSurf_incr[7] = -(1.0*coeff[1]*fSkin[7])-0.8660254037844386*coeff[1]*fSkin[5]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 

  return 0.;
}

