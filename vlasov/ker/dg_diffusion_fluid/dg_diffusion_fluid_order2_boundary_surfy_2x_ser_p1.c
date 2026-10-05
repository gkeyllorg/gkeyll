#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_2x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  double vol_incr[4] = {0.0}; 

  double edgeSurf_incr[4] = {0.0}; 
  double boundSurf_incr[4] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.5412658773652741*coeff[1]*fSkin[2])-0.5412658773652741*coeff[1]*fEdge[2]-0.5625*fSkin[0]*coeff[1]+0.5625*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(0.5412658773652741*coeff[1]*fSkin[3])-0.5412658773652741*coeff[1]*fEdge[3]-0.5625*coeff[1]*fSkin[1]+0.5625*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(1.4375*coeff[1]*fSkin[2])-0.4375*coeff[1]*fEdge[2]-1.4072912811497125*fSkin[0]*coeff[1]+0.5412658773652739*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(1.4375*coeff[1]*fSkin[3])-0.4375*coeff[1]*fEdge[3]-1.4072912811497125*coeff[1]*fSkin[1]+0.5412658773652739*coeff[1]*fEdge[1]; 

  boundSurf_incr[2] = 0.8660254037844386*fSkin[0]*coeff[1]-1.0*coeff[1]*fSkin[2]; 
  boundSurf_incr[3] = 0.8660254037844386*coeff[1]*fSkin[1]-1.0*coeff[1]*fSkin[3]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[1]*fSkin[2]+0.5412658773652741*coeff[1]*fEdge[2]-0.5625*fSkin[0]*coeff[1]+0.5625*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 0.5412658773652741*coeff[1]*fSkin[3]+0.5412658773652741*coeff[1]*fEdge[3]-0.5625*coeff[1]*fSkin[1]+0.5625*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(1.4375*coeff[1]*fSkin[2])-0.4375*coeff[1]*fEdge[2]+1.4072912811497125*fSkin[0]*coeff[1]-0.5412658773652739*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(1.4375*coeff[1]*fSkin[3])-0.4375*coeff[1]*fEdge[3]+1.4072912811497125*coeff[1]*fSkin[1]-0.5412658773652739*coeff[1]*fEdge[1]; 

  boundSurf_incr[2] = -(1.0*coeff[1]*fSkin[2])-0.8660254037844386*fSkin[0]*coeff[1]; 
  boundSurf_incr[3] = -(1.0*coeff[1]*fSkin[3])-0.8660254037844386*coeff[1]*fSkin[1]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 

  return 0.;
}

