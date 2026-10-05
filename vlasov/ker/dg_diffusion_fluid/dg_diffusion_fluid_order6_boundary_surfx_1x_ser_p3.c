#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  double vol_incr[4] = {0.0}; 

  double edgeSurf_incr[4] = {0.0}; 
  double boundSurf_incr[4] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(214.86394436340956*coeff[0]*fSkin[3])-214.86394436340956*coeff[0]*fEdge[3]-313.6609416875682*coeff[0]*fSkin[2]+313.6609416875682*coeff[0]*fEdge[2]-268.53553340784646*coeff[0]*fSkin[1]-268.53553340784646*coeff[0]*fEdge[1]-155.0390625*coeff[0]*fSkin[0]+155.0390625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(423.7092449203309*coeff[0]*fSkin[3])-320.60129178382454*coeff[0]*fEdge[3]-581.4013671669961*coeff[0]*fSkin[2]+505.15200753853765*coeff[0]*fEdge[2]-474.9609375*coeff[0]*fSkin[1]-455.2734375*coeff[0]*fEdge[1]-268.53553340784646*coeff[0]*fSkin[0]+268.53553340784646*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(577.1644241520037*coeff[0]*fSkin[3])-177.82903879277956*coeff[0]*fEdge[3]-670.60546875*coeff[0]*fSkin[2]+375.29296875*coeff[0]*fEdge[2]-471.79291270108683*coeff[0]*fSkin[1]-395.54355307262836*coeff[0]*fEdge[1]-250.37847099621678*coeff[0]*fSkin[0]+250.37847099621678*coeff[0]*fEdge[0]; 
  edgeSurf_incr[3] = -(319.39453125*coeff[0]*fSkin[3])+400.60546875*coeff[0]*fEdge[3]-3.119807698118848*coeff[0]*fSkin[2]-454.45198802599225*coeff[0]*fEdge[2]+234.40948720877657*coeff[0]*fSkin[1]+318.18469913218803*coeff[0]*fEdge[1]+159.52020111828895*coeff[0]*fSkin[0]-159.52020111828895*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = -(103.10795313650637*coeff[0]*fSkin[3])+76.2493596284585*coeff[0]*fSkin[2]-19.6875*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 399.3353853592242*coeff[0]*fSkin[3]-295.3125*coeff[0]*fSkin[2]+76.2493596284585*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = -(720.0*coeff[0]*fSkin[3])+457.5717957241111*coeff[0]*fSkin[2]-83.77521192341143*coeff[0]*fSkin[1]; 

  } else { 

  edgeSurf_incr[0] = 214.86394436340956*coeff[0]*fSkin[3]+214.86394436340956*coeff[0]*fEdge[3]-313.6609416875682*coeff[0]*fSkin[2]+313.6609416875682*coeff[0]*fEdge[2]+268.53553340784646*coeff[0]*fSkin[1]+268.53553340784646*coeff[0]*fEdge[1]-155.0390625*coeff[0]*fSkin[0]+155.0390625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(423.7092449203309*coeff[0]*fSkin[3])-320.60129178382454*coeff[0]*fEdge[3]+581.4013671669961*coeff[0]*fSkin[2]-505.15200753853765*coeff[0]*fEdge[2]-474.9609375*coeff[0]*fSkin[1]-455.2734375*coeff[0]*fEdge[1]+268.53553340784646*coeff[0]*fSkin[0]-268.53553340784646*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 577.1644241520037*coeff[0]*fSkin[3]+177.82903879277956*coeff[0]*fEdge[3]-670.60546875*coeff[0]*fSkin[2]+375.29296875*coeff[0]*fEdge[2]+471.79291270108683*coeff[0]*fSkin[1]+395.54355307262836*coeff[0]*fEdge[1]-250.37847099621678*coeff[0]*fSkin[0]+250.37847099621678*coeff[0]*fEdge[0]; 
  edgeSurf_incr[3] = -(319.39453125*coeff[0]*fSkin[3])+400.60546875*coeff[0]*fEdge[3]+3.119807698118848*coeff[0]*fSkin[2]+454.45198802599225*coeff[0]*fEdge[2]+234.40948720877657*coeff[0]*fSkin[1]+318.18469913218803*coeff[0]*fEdge[1]-159.52020111828895*coeff[0]*fSkin[0]+159.52020111828895*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = -(103.10795313650637*coeff[0]*fSkin[3])-76.2493596284585*coeff[0]*fSkin[2]-19.6875*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = -(399.3353853592242*coeff[0]*fSkin[3])-295.3125*coeff[0]*fSkin[2]-76.2493596284585*coeff[0]*fSkin[1]; 
  boundSurf_incr[3] = -(720.0*coeff[0]*fSkin[3])-457.5717957241111*coeff[0]*fSkin[2]-83.77521192341143*coeff[0]*fSkin[1]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 

  return 0.;
}

