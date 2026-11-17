#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_1x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[4] = {0.0}; 
  vol_incr[2] = 6.708203932499369*coeff[0]*fSkin[0]; 
  vol_incr[3] = 22.9128784747792*coeff[0]*fSkin[1]; 

  double edgeSurf_incr[4] = {0.0}; 
  double boundSurf_incr[4] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.6821077598838398*coeff[0]*fSkin[3])-0.6821077598838398*coeff[0]*fEdge[3]-1.5547660156053227*coeff[0]*fSkin[2]+1.5547660156053227*coeff[0]*fEdge[2]-1.935025511580855*coeff[0]*fSkin[1]-1.935025511580855*coeff[0]*fEdge[1]-1.3671875*coeff[0]*fSkin[0]+1.3671875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.508772131709791*coeff[0]*fSkin[3])-0.8541184610018139*coeff[0]*fEdge[3]-3.2980873807547533*coeff[0]*fSkin[2]+2.0877800850649355*coeff[0]*fEdge[2]-4.0078125*coeff[0]*fSkin[1]-2.6953125*coeff[0]*fEdge[1]-2.801050915365293*coeff[0]*fSkin[0]+1.935025511580855*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(2.792970701173145*coeff[0]*fSkin[3])-0.257507936987595*coeff[0]*fEdge[3]-5.8203125*coeff[0]*fSkin[2]+1.1328125*coeff[0]*fEdge[2]-6.868493903039716*coeff[0]*fSkin[1]-1.785203261142482*coeff[0]*fEdge[1]-4.734175171112836*coeff[0]*fSkin[0]+1.380073204863152*coeff[0]*fEdge[0]; 
  edgeSurf_incr[3] = -(4.8046875*coeff[0]*fSkin[3])+1.1953125*coeff[0]*fEdge[3]-9.659849020842344*coeff[0]*fSkin[2]-1.4328005724694384*coeff[0]*fEdge[2]-11.134226883838018*coeff[0]*fSkin[1]+0.8950343154210625*coeff[0]*fEdge[1]-7.585865087193007*coeff[0]*fSkin[0]-0.3513888460007659*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = -(0.9165151389911681*coeff[0]*fSkin[3])+1.3555441711725957*coeff[0]*fSkin[2]-1.35*coeff[0]*fSkin[1]+0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[2] = 3.5496478698597698*coeff[0]*fSkin[3]-5.25*coeff[0]*fSkin[2]+5.2285275173800105*coeff[0]*fSkin[1]-3.3541019662496847*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = -(8.4*coeff[0]*fSkin[3])+12.423767544509195*coeff[0]*fSkin[2]-12.372954376380765*coeff[0]*fSkin[1]+7.937253933193772*coeff[0]*fSkin[0]; 

  } else { 

  edgeSurf_incr[0] = 0.6821077598838398*coeff[0]*fSkin[3]+0.6821077598838398*coeff[0]*fEdge[3]-1.5547660156053227*coeff[0]*fSkin[2]+1.5547660156053227*coeff[0]*fEdge[2]+1.935025511580855*coeff[0]*fSkin[1]+1.935025511580855*coeff[0]*fEdge[1]-1.3671875*coeff[0]*fSkin[0]+1.3671875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.508772131709791*coeff[0]*fSkin[3])-0.8541184610018139*coeff[0]*fEdge[3]+3.2980873807547533*coeff[0]*fSkin[2]-2.0877800850649355*coeff[0]*fEdge[2]-4.0078125*coeff[0]*fSkin[1]-2.6953125*coeff[0]*fEdge[1]+2.801050915365293*coeff[0]*fSkin[0]-1.935025511580855*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 2.792970701173145*coeff[0]*fSkin[3]+0.257507936987595*coeff[0]*fEdge[3]-5.8203125*coeff[0]*fSkin[2]+1.1328125*coeff[0]*fEdge[2]+6.868493903039716*coeff[0]*fSkin[1]+1.785203261142482*coeff[0]*fEdge[1]-4.734175171112836*coeff[0]*fSkin[0]+1.380073204863152*coeff[0]*fEdge[0]; 
  edgeSurf_incr[3] = -(4.8046875*coeff[0]*fSkin[3])+1.1953125*coeff[0]*fEdge[3]+9.659849020842344*coeff[0]*fSkin[2]+1.4328005724694384*coeff[0]*fEdge[2]-11.134226883838018*coeff[0]*fSkin[1]+0.8950343154210625*coeff[0]*fEdge[1]+7.585865087193007*coeff[0]*fSkin[0]+0.3513888460007659*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = -(0.9165151389911681*coeff[0]*fSkin[3])-1.3555441711725957*coeff[0]*fSkin[2]-1.35*coeff[0]*fSkin[1]-0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[2] = -(3.5496478698597698*coeff[0]*fSkin[3])-5.25*coeff[0]*fSkin[2]-5.2285275173800105*coeff[0]*fSkin[1]-3.3541019662496847*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = -(8.4*coeff[0]*fSkin[3])-12.423767544509195*coeff[0]*fSkin[2]-12.372954376380765*coeff[0]*fSkin[1]-7.937253933193772*coeff[0]*fSkin[0]; 

  }

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 

  return 0.;
}

