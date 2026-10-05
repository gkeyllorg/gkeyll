#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfx_2x_ser_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[12] = {0.0}; 
  vol_incr[4] = 6.708203932499369*coeff[0]*fSkin[0]; 
  vol_incr[6] = 6.7082039324993685*coeff[0]*fSkin[2]; 
  vol_incr[8] = 22.9128784747792*coeff[0]*fSkin[1]; 
  vol_incr[10] = 22.9128784747792*coeff[0]*fSkin[3]; 

  double edgeSurf_incr[12] = {0.0}; 
  double boundSurf_incr[12] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.6821077598838398*coeff[0]*fSkin[8])-0.6821077598838398*coeff[0]*fEdge[8]-1.5547660156053227*coeff[0]*fSkin[4]+1.5547660156053227*coeff[0]*fEdge[4]-1.935025511580855*coeff[0]*fSkin[1]-1.935025511580855*coeff[0]*fEdge[1]-1.3671875*coeff[0]*fSkin[0]+1.3671875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.508772131709791*coeff[0]*fSkin[8])-0.8541184610018139*coeff[0]*fEdge[8]-3.2980873807547533*coeff[0]*fSkin[4]+2.0877800850649355*coeff[0]*fEdge[4]-4.0078125*coeff[0]*fSkin[1]-2.6953125*coeff[0]*fEdge[1]-2.801050915365293*coeff[0]*fSkin[0]+1.935025511580855*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = -(0.6821077598838398*coeff[0]*fSkin[10])-0.6821077598838398*coeff[0]*fEdge[10]-1.5547660156053227*coeff[0]*fSkin[6]+1.5547660156053227*coeff[0]*fEdge[6]-1.935025511580855*coeff[0]*fSkin[3]-1.935025511580855*coeff[0]*fEdge[3]-1.3671875*coeff[0]*fSkin[2]+1.3671875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(1.508772131709791*coeff[0]*fSkin[10])-0.8541184610018139*coeff[0]*fEdge[10]-3.2980873807547537*coeff[0]*fSkin[6]+2.087780085064936*coeff[0]*fEdge[6]-4.0078125*coeff[0]*fSkin[3]-2.6953125*coeff[0]*fEdge[3]-2.801050915365293*coeff[0]*fSkin[2]+1.935025511580855*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = -(2.792970701173145*coeff[0]*fSkin[8])-0.257507936987595*coeff[0]*fEdge[8]-5.8203125*coeff[0]*fSkin[4]+1.1328125*coeff[0]*fEdge[4]-6.868493903039716*coeff[0]*fSkin[1]-1.785203261142482*coeff[0]*fEdge[1]-4.734175171112836*coeff[0]*fSkin[0]+1.380073204863152*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = -(1.935025511580855*coeff[0]*fSkin[7])-1.935025511580855*coeff[0]*fEdge[7]-1.3671875*coeff[0]*fSkin[5]+1.3671875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = -(2.792970701173145*coeff[0]*fSkin[10])-0.25750793698759494*coeff[0]*fEdge[10]-5.8203125*coeff[0]*fSkin[6]+1.1328125*coeff[0]*fEdge[6]-6.868493903039716*coeff[0]*fSkin[3]-1.7852032611424813*coeff[0]*fEdge[3]-4.734175171112836*coeff[0]*fSkin[2]+1.3800732048631523*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = -(4.0078125*coeff[0]*fSkin[7])-2.6953125*coeff[0]*fEdge[7]-2.8010509153652943*coeff[0]*fSkin[5]+1.9350255115808557*coeff[0]*fEdge[5]; 
  edgeSurf_incr[8] = -(4.8046875*coeff[0]*fSkin[8])+1.1953125*coeff[0]*fEdge[8]-9.659849020842344*coeff[0]*fSkin[4]-1.4328005724694384*coeff[0]*fEdge[4]-11.134226883838018*coeff[0]*fSkin[1]+0.8950343154210625*coeff[0]*fEdge[1]-7.585865087193007*coeff[0]*fSkin[0]-0.3513888460007659*coeff[0]*fEdge[0]; 
  edgeSurf_incr[9] = -(1.9350255115808548*coeff[0]*fSkin[11])-1.9350255115808548*coeff[0]*fEdge[11]-1.3671875*coeff[0]*fSkin[9]+1.3671875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = -(4.8046875*coeff[0]*fSkin[10])+1.1953125*coeff[0]*fEdge[10]-9.659849020842342*coeff[0]*fSkin[6]-1.4328005724694384*coeff[0]*fEdge[6]-11.134226883838018*coeff[0]*fSkin[3]+0.8950343154210616*coeff[0]*fEdge[3]-7.585865087193007*coeff[0]*fSkin[2]-0.35138884600076503*coeff[0]*fEdge[2]; 
  edgeSurf_incr[11] = -(4.0078125*coeff[0]*fSkin[11])-2.6953125*coeff[0]*fEdge[11]-2.801050915365294*coeff[0]*fSkin[9]+1.9350255115808555*coeff[0]*fEdge[9]; 

  boundSurf_incr[1] = -(0.9165151389911681*coeff[0]*fSkin[8])+1.3555441711725957*coeff[0]*fSkin[4]-1.35*coeff[0]*fSkin[1]+0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = -(0.916515138991168*coeff[0]*fSkin[10])+1.355544171172596*coeff[0]*fSkin[6]-1.35*coeff[0]*fSkin[3]+0.8660254037844386*coeff[0]*fSkin[2]; 
  boundSurf_incr[4] = 3.5496478698597698*coeff[0]*fSkin[8]-5.25*coeff[0]*fSkin[4]+5.2285275173800105*coeff[0]*fSkin[1]-3.3541019662496847*coeff[0]*fSkin[0]; 
  boundSurf_incr[6] = 3.5496478698597698*coeff[0]*fSkin[10]-5.25*coeff[0]*fSkin[6]+5.228527517380013*coeff[0]*fSkin[3]-3.3541019662496843*coeff[0]*fSkin[2]; 
  boundSurf_incr[7] = 0.8660254037844387*coeff[0]*fSkin[5]-1.35*coeff[0]*fSkin[7]; 
  boundSurf_incr[8] = -(8.4*coeff[0]*fSkin[8])+12.423767544509195*coeff[0]*fSkin[4]-12.372954376380765*coeff[0]*fSkin[1]+7.937253933193772*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = -(8.4*coeff[0]*fSkin[10])+12.423767544509193*coeff[0]*fSkin[6]-12.372954376380768*coeff[0]*fSkin[3]+7.937253933193771*coeff[0]*fSkin[2]; 
  boundSurf_incr[11] = 0.8660254037844386*coeff[0]*fSkin[9]-1.35*coeff[0]*fSkin[11]; 

  } else { 

  edgeSurf_incr[0] = 0.6821077598838398*coeff[0]*fSkin[8]+0.6821077598838398*coeff[0]*fEdge[8]-1.5547660156053227*coeff[0]*fSkin[4]+1.5547660156053227*coeff[0]*fEdge[4]+1.935025511580855*coeff[0]*fSkin[1]+1.935025511580855*coeff[0]*fEdge[1]-1.3671875*coeff[0]*fSkin[0]+1.3671875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = -(1.508772131709791*coeff[0]*fSkin[8])-0.8541184610018139*coeff[0]*fEdge[8]+3.2980873807547533*coeff[0]*fSkin[4]-2.0877800850649355*coeff[0]*fEdge[4]-4.0078125*coeff[0]*fSkin[1]-2.6953125*coeff[0]*fEdge[1]+2.801050915365293*coeff[0]*fSkin[0]-1.935025511580855*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.6821077598838398*coeff[0]*fSkin[10]+0.6821077598838398*coeff[0]*fEdge[10]-1.5547660156053227*coeff[0]*fSkin[6]+1.5547660156053227*coeff[0]*fEdge[6]+1.935025511580855*coeff[0]*fSkin[3]+1.935025511580855*coeff[0]*fEdge[3]-1.3671875*coeff[0]*fSkin[2]+1.3671875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = -(1.508772131709791*coeff[0]*fSkin[10])-0.8541184610018139*coeff[0]*fEdge[10]+3.2980873807547537*coeff[0]*fSkin[6]-2.087780085064936*coeff[0]*fEdge[6]-4.0078125*coeff[0]*fSkin[3]-2.6953125*coeff[0]*fEdge[3]+2.801050915365293*coeff[0]*fSkin[2]-1.935025511580855*coeff[0]*fEdge[2]; 
  edgeSurf_incr[4] = 2.792970701173145*coeff[0]*fSkin[8]+0.257507936987595*coeff[0]*fEdge[8]-5.8203125*coeff[0]*fSkin[4]+1.1328125*coeff[0]*fEdge[4]+6.868493903039716*coeff[0]*fSkin[1]+1.785203261142482*coeff[0]*fEdge[1]-4.734175171112836*coeff[0]*fSkin[0]+1.380073204863152*coeff[0]*fEdge[0]; 
  edgeSurf_incr[5] = 1.935025511580855*coeff[0]*fSkin[7]+1.935025511580855*coeff[0]*fEdge[7]-1.3671875*coeff[0]*fSkin[5]+1.3671875*coeff[0]*fEdge[5]; 
  edgeSurf_incr[6] = 2.792970701173145*coeff[0]*fSkin[10]+0.25750793698759494*coeff[0]*fEdge[10]-5.8203125*coeff[0]*fSkin[6]+1.1328125*coeff[0]*fEdge[6]+6.868493903039716*coeff[0]*fSkin[3]+1.7852032611424813*coeff[0]*fEdge[3]-4.734175171112836*coeff[0]*fSkin[2]+1.3800732048631523*coeff[0]*fEdge[2]; 
  edgeSurf_incr[7] = -(4.0078125*coeff[0]*fSkin[7])-2.6953125*coeff[0]*fEdge[7]+2.8010509153652943*coeff[0]*fSkin[5]-1.9350255115808557*coeff[0]*fEdge[5]; 
  edgeSurf_incr[8] = -(4.8046875*coeff[0]*fSkin[8])+1.1953125*coeff[0]*fEdge[8]+9.659849020842344*coeff[0]*fSkin[4]+1.4328005724694384*coeff[0]*fEdge[4]-11.134226883838018*coeff[0]*fSkin[1]+0.8950343154210625*coeff[0]*fEdge[1]+7.585865087193007*coeff[0]*fSkin[0]+0.3513888460007659*coeff[0]*fEdge[0]; 
  edgeSurf_incr[9] = 1.9350255115808548*coeff[0]*fSkin[11]+1.9350255115808548*coeff[0]*fEdge[11]-1.3671875*coeff[0]*fSkin[9]+1.3671875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = -(4.8046875*coeff[0]*fSkin[10])+1.1953125*coeff[0]*fEdge[10]+9.659849020842342*coeff[0]*fSkin[6]+1.4328005724694384*coeff[0]*fEdge[6]-11.134226883838018*coeff[0]*fSkin[3]+0.8950343154210616*coeff[0]*fEdge[3]+7.585865087193007*coeff[0]*fSkin[2]+0.35138884600076503*coeff[0]*fEdge[2]; 
  edgeSurf_incr[11] = -(4.0078125*coeff[0]*fSkin[11])-2.6953125*coeff[0]*fEdge[11]+2.801050915365294*coeff[0]*fSkin[9]-1.9350255115808555*coeff[0]*fEdge[9]; 

  boundSurf_incr[1] = -(0.9165151389911681*coeff[0]*fSkin[8])-1.3555441711725957*coeff[0]*fSkin[4]-1.35*coeff[0]*fSkin[1]-0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[3] = -(0.916515138991168*coeff[0]*fSkin[10])-1.355544171172596*coeff[0]*fSkin[6]-1.35*coeff[0]*fSkin[3]-0.8660254037844386*coeff[0]*fSkin[2]; 
  boundSurf_incr[4] = -(3.5496478698597698*coeff[0]*fSkin[8])-5.25*coeff[0]*fSkin[4]-5.2285275173800105*coeff[0]*fSkin[1]-3.3541019662496847*coeff[0]*fSkin[0]; 
  boundSurf_incr[6] = -(3.5496478698597698*coeff[0]*fSkin[10])-5.25*coeff[0]*fSkin[6]-5.228527517380013*coeff[0]*fSkin[3]-3.3541019662496843*coeff[0]*fSkin[2]; 
  boundSurf_incr[7] = -(1.35*coeff[0]*fSkin[7])-0.8660254037844387*coeff[0]*fSkin[5]; 
  boundSurf_incr[8] = -(8.4*coeff[0]*fSkin[8])-12.423767544509195*coeff[0]*fSkin[4]-12.372954376380765*coeff[0]*fSkin[1]-7.937253933193772*coeff[0]*fSkin[0]; 
  boundSurf_incr[10] = -(8.4*coeff[0]*fSkin[10])-12.423767544509193*coeff[0]*fSkin[6]-12.372954376380768*coeff[0]*fSkin[3]-7.937253933193771*coeff[0]*fSkin[2]; 
  boundSurf_incr[11] = -(1.35*coeff[0]*fSkin[11])-0.8660254037844386*coeff[0]*fSkin[9]; 

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

  return 0.;
}

