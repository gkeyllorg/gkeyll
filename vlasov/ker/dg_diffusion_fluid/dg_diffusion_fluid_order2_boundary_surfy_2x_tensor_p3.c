#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order2_boundary_surfy_2x_tensor_p3_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  double vol_incr[16] = {0.0}; 
  vol_incr[5] = 6.708203932499369*fSkin[0]*coeff[1]; 
  vol_incr[7] = 6.7082039324993685*coeff[1]*fSkin[1]; 
  vol_incr[9] = 22.9128784747792*coeff[1]*fSkin[2]; 
  vol_incr[10] = 6.708203932499369*coeff[1]*fSkin[4]; 
  vol_incr[12] = 22.9128784747792*coeff[1]*fSkin[3]; 
  vol_incr[13] = 6.708203932499368*coeff[1]*fSkin[8]; 
  vol_incr[14] = 22.912878474779202*coeff[1]*fSkin[6]; 
  vol_incr[15] = 22.9128784747792*coeff[1]*fSkin[11]; 

  double edgeSurf_incr[16] = {0.0}; 
  double boundSurf_incr[16] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = -(0.6821077598838398*coeff[1]*fSkin[9])-0.6821077598838398*coeff[1]*fEdge[9]-1.5547660156053227*coeff[1]*fSkin[5]+1.5547660156053227*coeff[1]*fEdge[5]-1.935025511580855*coeff[1]*fSkin[2]-1.935025511580855*coeff[1]*fEdge[2]-1.3671875*fSkin[0]*coeff[1]+1.3671875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = -(0.6821077598838398*coeff[1]*fSkin[12])-0.6821077598838398*coeff[1]*fEdge[12]-1.5547660156053227*coeff[1]*fSkin[7]+1.5547660156053227*coeff[1]*fEdge[7]-1.935025511580855*coeff[1]*fSkin[3]-1.935025511580855*coeff[1]*fEdge[3]-1.3671875*coeff[1]*fSkin[1]+1.3671875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(1.508772131709791*coeff[1]*fSkin[9])-0.8541184610018139*coeff[1]*fEdge[9]-3.2980873807547533*coeff[1]*fSkin[5]+2.0877800850649355*coeff[1]*fEdge[5]-4.0078125*coeff[1]*fSkin[2]-2.6953125*coeff[1]*fEdge[2]-2.801050915365293*fSkin[0]*coeff[1]+1.935025511580855*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(1.508772131709791*coeff[1]*fSkin[12])-0.8541184610018139*coeff[1]*fEdge[12]-3.2980873807547537*coeff[1]*fSkin[7]+2.087780085064936*coeff[1]*fEdge[7]-4.0078125*coeff[1]*fSkin[3]-2.6953125*coeff[1]*fEdge[3]-2.801050915365293*coeff[1]*fSkin[1]+1.935025511580855*coeff[1]*fEdge[1]; 
  edgeSurf_incr[4] = -(0.6821077598838398*coeff[1]*fSkin[14])-0.6821077598838398*coeff[1]*fEdge[14]-1.5547660156053227*coeff[1]*fSkin[10]+1.5547660156053227*coeff[1]*fEdge[10]-1.935025511580855*coeff[1]*fSkin[6]-1.935025511580855*coeff[1]*fEdge[6]-1.3671875*coeff[1]*fSkin[4]+1.3671875*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = -(2.792970701173145*coeff[1]*fSkin[9])-0.257507936987595*coeff[1]*fEdge[9]-5.8203125*coeff[1]*fSkin[5]+1.1328125*coeff[1]*fEdge[5]-6.868493903039716*coeff[1]*fSkin[2]-1.785203261142482*coeff[1]*fEdge[2]-4.734175171112836*fSkin[0]*coeff[1]+1.380073204863152*fEdge[0]*coeff[1]; 
  edgeSurf_incr[6] = -(1.508772131709791*coeff[1]*fSkin[14])-0.8541184610018139*coeff[1]*fEdge[14]-3.2980873807547537*coeff[1]*fSkin[10]+2.087780085064936*coeff[1]*fEdge[10]-4.0078125*coeff[1]*fSkin[6]-2.6953125*coeff[1]*fEdge[6]-2.8010509153652943*coeff[1]*fSkin[4]+1.9350255115808557*coeff[1]*fEdge[4]; 
  edgeSurf_incr[7] = -(2.792970701173145*coeff[1]*fSkin[12])-0.25750793698759494*coeff[1]*fEdge[12]-5.8203125*coeff[1]*fSkin[7]+1.1328125*coeff[1]*fEdge[7]-6.868493903039716*coeff[1]*fSkin[3]-1.7852032611424813*coeff[1]*fEdge[3]-4.734175171112836*coeff[1]*fSkin[1]+1.3800732048631523*coeff[1]*fEdge[1]; 
  edgeSurf_incr[8] = -(0.6821077598838398*coeff[1]*fSkin[15])-0.6821077598838398*coeff[1]*fEdge[15]-1.5547660156053225*coeff[1]*fSkin[13]+1.5547660156053225*coeff[1]*fEdge[13]-1.9350255115808548*coeff[1]*fSkin[11]-1.9350255115808548*coeff[1]*fEdge[11]-1.3671875*coeff[1]*fSkin[8]+1.3671875*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = -(4.8046875*coeff[1]*fSkin[9])+1.1953125*coeff[1]*fEdge[9]-9.659849020842344*coeff[1]*fSkin[5]-1.4328005724694384*coeff[1]*fEdge[5]-11.134226883838018*coeff[1]*fSkin[2]+0.8950343154210625*coeff[1]*fEdge[2]-7.585865087193007*fSkin[0]*coeff[1]-0.3513888460007659*fEdge[0]*coeff[1]; 
  edgeSurf_incr[10] = -(2.792970701173145*coeff[1]*fSkin[14])-0.2575079369875949*coeff[1]*fEdge[14]-5.8203125*coeff[1]*fSkin[10]+1.1328125*coeff[1]*fEdge[10]-6.868493903039716*coeff[1]*fSkin[6]-1.7852032611424813*coeff[1]*fEdge[6]-4.734175171112836*coeff[1]*fSkin[4]+1.380073204863152*coeff[1]*fEdge[4]; 
  edgeSurf_incr[11] = -(1.508772131709791*coeff[1]*fSkin[15])-0.8541184610018139*coeff[1]*fEdge[15]-3.2980873807547533*coeff[1]*fSkin[13]+2.087780085064936*coeff[1]*fEdge[13]-4.0078125*coeff[1]*fSkin[11]-2.6953125*coeff[1]*fEdge[11]-2.801050915365294*coeff[1]*fSkin[8]+1.9350255115808555*coeff[1]*fEdge[8]; 
  edgeSurf_incr[12] = -(4.8046875*coeff[1]*fSkin[12])+1.1953125*coeff[1]*fEdge[12]-9.659849020842342*coeff[1]*fSkin[7]-1.4328005724694384*coeff[1]*fEdge[7]-11.134226883838018*coeff[1]*fSkin[3]+0.8950343154210616*coeff[1]*fEdge[3]-7.585865087193007*coeff[1]*fSkin[1]-0.35138884600076503*coeff[1]*fEdge[1]; 
  edgeSurf_incr[13] = -(2.792970701173145*coeff[1]*fSkin[15])-0.2575079369875949*coeff[1]*fEdge[15]-5.8203125*coeff[1]*fSkin[13]+1.1328125*coeff[1]*fEdge[13]-6.868493903039715*coeff[1]*fSkin[11]-1.785203261142481*coeff[1]*fEdge[11]-4.734175171112836*coeff[1]*fSkin[8]+1.380073204863152*coeff[1]*fEdge[8]; 
  edgeSurf_incr[14] = -(4.8046875*coeff[1]*fSkin[14])+1.1953125*coeff[1]*fEdge[14]-9.659849020842342*coeff[1]*fSkin[10]-1.4328005724694384*coeff[1]*fEdge[10]-11.13422688383802*coeff[1]*fSkin[6]+0.8950343154210625*coeff[1]*fEdge[6]-7.585865087193007*coeff[1]*fSkin[4]-0.35138884600076503*coeff[1]*fEdge[4]; 
  edgeSurf_incr[15] = -(4.8046875*coeff[1]*fSkin[15])+1.1953125*coeff[1]*fEdge[15]-9.659849020842342*coeff[1]*fSkin[13]-1.4328005724694384*coeff[1]*fEdge[13]-11.134226883838018*coeff[1]*fSkin[11]+0.8950343154210616*coeff[1]*fEdge[11]-7.585865087193007*coeff[1]*fSkin[8]-0.3513888460007659*coeff[1]*fEdge[8]; 

  boundSurf_incr[2] = -(0.9165151389911681*coeff[1]*fSkin[9])+1.3555441711725957*coeff[1]*fSkin[5]-1.35*coeff[1]*fSkin[2]+0.8660254037844386*fSkin[0]*coeff[1]; 
  boundSurf_incr[3] = -(0.916515138991168*coeff[1]*fSkin[12])+1.355544171172596*coeff[1]*fSkin[7]-1.35*coeff[1]*fSkin[3]+0.8660254037844386*coeff[1]*fSkin[1]; 
  boundSurf_incr[5] = 3.5496478698597698*coeff[1]*fSkin[9]-5.25*coeff[1]*fSkin[5]+5.2285275173800105*coeff[1]*fSkin[2]-3.3541019662496847*fSkin[0]*coeff[1]; 
  boundSurf_incr[6] = -(0.916515138991168*coeff[1]*fSkin[14])+1.355544171172596*coeff[1]*fSkin[10]-1.35*coeff[1]*fSkin[6]+0.8660254037844387*coeff[1]*fSkin[4]; 
  boundSurf_incr[7] = 3.5496478698597698*coeff[1]*fSkin[12]-5.25*coeff[1]*fSkin[7]+5.228527517380013*coeff[1]*fSkin[3]-3.3541019662496843*coeff[1]*fSkin[1]; 
  boundSurf_incr[9] = -(8.4*coeff[1]*fSkin[9])+12.423767544509195*coeff[1]*fSkin[5]-12.372954376380765*coeff[1]*fSkin[2]+7.937253933193772*fSkin[0]*coeff[1]; 
  boundSurf_incr[10] = 3.5496478698597698*coeff[1]*fSkin[14]-5.25*coeff[1]*fSkin[10]+5.228527517380013*coeff[1]*fSkin[6]-3.3541019662496847*coeff[1]*fSkin[4]; 
  boundSurf_incr[11] = -(0.916515138991168*coeff[1]*fSkin[15])+1.355544171172596*coeff[1]*fSkin[13]-1.35*coeff[1]*fSkin[11]+0.8660254037844386*coeff[1]*fSkin[8]; 
  boundSurf_incr[12] = -(8.4*coeff[1]*fSkin[12])+12.423767544509193*coeff[1]*fSkin[7]-12.372954376380768*coeff[1]*fSkin[3]+7.937253933193771*coeff[1]*fSkin[1]; 
  boundSurf_incr[13] = 3.5496478698597698*coeff[1]*fSkin[15]-5.25*coeff[1]*fSkin[13]+5.228527517380012*coeff[1]*fSkin[11]-3.354101966249684*coeff[1]*fSkin[8]; 
  boundSurf_incr[14] = -(8.4*coeff[1]*fSkin[14])+12.423767544509193*coeff[1]*fSkin[10]-12.37295437638077*coeff[1]*fSkin[6]+7.937253933193772*coeff[1]*fSkin[4]; 
  boundSurf_incr[15] = -(8.4*coeff[1]*fSkin[15])+12.423767544509195*coeff[1]*fSkin[13]-12.372954376380768*coeff[1]*fSkin[11]+7.937253933193772*coeff[1]*fSkin[8]; 

  } else { 

  edgeSurf_incr[0] = 0.6821077598838398*coeff[1]*fSkin[9]+0.6821077598838398*coeff[1]*fEdge[9]-1.5547660156053227*coeff[1]*fSkin[5]+1.5547660156053227*coeff[1]*fEdge[5]+1.935025511580855*coeff[1]*fSkin[2]+1.935025511580855*coeff[1]*fEdge[2]-1.3671875*fSkin[0]*coeff[1]+1.3671875*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 0.6821077598838398*coeff[1]*fSkin[12]+0.6821077598838398*coeff[1]*fEdge[12]-1.5547660156053227*coeff[1]*fSkin[7]+1.5547660156053227*coeff[1]*fEdge[7]+1.935025511580855*coeff[1]*fSkin[3]+1.935025511580855*coeff[1]*fEdge[3]-1.3671875*coeff[1]*fSkin[1]+1.3671875*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = -(1.508772131709791*coeff[1]*fSkin[9])-0.8541184610018139*coeff[1]*fEdge[9]+3.2980873807547533*coeff[1]*fSkin[5]-2.0877800850649355*coeff[1]*fEdge[5]-4.0078125*coeff[1]*fSkin[2]-2.6953125*coeff[1]*fEdge[2]+2.801050915365293*fSkin[0]*coeff[1]-1.935025511580855*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = -(1.508772131709791*coeff[1]*fSkin[12])-0.8541184610018139*coeff[1]*fEdge[12]+3.2980873807547537*coeff[1]*fSkin[7]-2.087780085064936*coeff[1]*fEdge[7]-4.0078125*coeff[1]*fSkin[3]-2.6953125*coeff[1]*fEdge[3]+2.801050915365293*coeff[1]*fSkin[1]-1.935025511580855*coeff[1]*fEdge[1]; 
  edgeSurf_incr[4] = 0.6821077598838398*coeff[1]*fSkin[14]+0.6821077598838398*coeff[1]*fEdge[14]-1.5547660156053227*coeff[1]*fSkin[10]+1.5547660156053227*coeff[1]*fEdge[10]+1.935025511580855*coeff[1]*fSkin[6]+1.935025511580855*coeff[1]*fEdge[6]-1.3671875*coeff[1]*fSkin[4]+1.3671875*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = 2.792970701173145*coeff[1]*fSkin[9]+0.257507936987595*coeff[1]*fEdge[9]-5.8203125*coeff[1]*fSkin[5]+1.1328125*coeff[1]*fEdge[5]+6.868493903039716*coeff[1]*fSkin[2]+1.785203261142482*coeff[1]*fEdge[2]-4.734175171112836*fSkin[0]*coeff[1]+1.380073204863152*fEdge[0]*coeff[1]; 
  edgeSurf_incr[6] = -(1.508772131709791*coeff[1]*fSkin[14])-0.8541184610018139*coeff[1]*fEdge[14]+3.2980873807547537*coeff[1]*fSkin[10]-2.087780085064936*coeff[1]*fEdge[10]-4.0078125*coeff[1]*fSkin[6]-2.6953125*coeff[1]*fEdge[6]+2.8010509153652943*coeff[1]*fSkin[4]-1.9350255115808557*coeff[1]*fEdge[4]; 
  edgeSurf_incr[7] = 2.792970701173145*coeff[1]*fSkin[12]+0.25750793698759494*coeff[1]*fEdge[12]-5.8203125*coeff[1]*fSkin[7]+1.1328125*coeff[1]*fEdge[7]+6.868493903039716*coeff[1]*fSkin[3]+1.7852032611424813*coeff[1]*fEdge[3]-4.734175171112836*coeff[1]*fSkin[1]+1.3800732048631523*coeff[1]*fEdge[1]; 
  edgeSurf_incr[8] = 0.6821077598838398*coeff[1]*fSkin[15]+0.6821077598838398*coeff[1]*fEdge[15]-1.5547660156053225*coeff[1]*fSkin[13]+1.5547660156053225*coeff[1]*fEdge[13]+1.9350255115808548*coeff[1]*fSkin[11]+1.9350255115808548*coeff[1]*fEdge[11]-1.3671875*coeff[1]*fSkin[8]+1.3671875*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = -(4.8046875*coeff[1]*fSkin[9])+1.1953125*coeff[1]*fEdge[9]+9.659849020842344*coeff[1]*fSkin[5]+1.4328005724694384*coeff[1]*fEdge[5]-11.134226883838018*coeff[1]*fSkin[2]+0.8950343154210625*coeff[1]*fEdge[2]+7.585865087193007*fSkin[0]*coeff[1]+0.3513888460007659*fEdge[0]*coeff[1]; 
  edgeSurf_incr[10] = 2.792970701173145*coeff[1]*fSkin[14]+0.2575079369875949*coeff[1]*fEdge[14]-5.8203125*coeff[1]*fSkin[10]+1.1328125*coeff[1]*fEdge[10]+6.868493903039716*coeff[1]*fSkin[6]+1.7852032611424813*coeff[1]*fEdge[6]-4.734175171112836*coeff[1]*fSkin[4]+1.380073204863152*coeff[1]*fEdge[4]; 
  edgeSurf_incr[11] = -(1.508772131709791*coeff[1]*fSkin[15])-0.8541184610018139*coeff[1]*fEdge[15]+3.2980873807547533*coeff[1]*fSkin[13]-2.087780085064936*coeff[1]*fEdge[13]-4.0078125*coeff[1]*fSkin[11]-2.6953125*coeff[1]*fEdge[11]+2.801050915365294*coeff[1]*fSkin[8]-1.9350255115808555*coeff[1]*fEdge[8]; 
  edgeSurf_incr[12] = -(4.8046875*coeff[1]*fSkin[12])+1.1953125*coeff[1]*fEdge[12]+9.659849020842342*coeff[1]*fSkin[7]+1.4328005724694384*coeff[1]*fEdge[7]-11.134226883838018*coeff[1]*fSkin[3]+0.8950343154210616*coeff[1]*fEdge[3]+7.585865087193007*coeff[1]*fSkin[1]+0.35138884600076503*coeff[1]*fEdge[1]; 
  edgeSurf_incr[13] = 2.792970701173145*coeff[1]*fSkin[15]+0.2575079369875949*coeff[1]*fEdge[15]-5.8203125*coeff[1]*fSkin[13]+1.1328125*coeff[1]*fEdge[13]+6.868493903039715*coeff[1]*fSkin[11]+1.785203261142481*coeff[1]*fEdge[11]-4.734175171112836*coeff[1]*fSkin[8]+1.380073204863152*coeff[1]*fEdge[8]; 
  edgeSurf_incr[14] = -(4.8046875*coeff[1]*fSkin[14])+1.1953125*coeff[1]*fEdge[14]+9.659849020842342*coeff[1]*fSkin[10]+1.4328005724694384*coeff[1]*fEdge[10]-11.13422688383802*coeff[1]*fSkin[6]+0.8950343154210625*coeff[1]*fEdge[6]+7.585865087193007*coeff[1]*fSkin[4]+0.35138884600076503*coeff[1]*fEdge[4]; 
  edgeSurf_incr[15] = -(4.8046875*coeff[1]*fSkin[15])+1.1953125*coeff[1]*fEdge[15]+9.659849020842342*coeff[1]*fSkin[13]+1.4328005724694384*coeff[1]*fEdge[13]-11.134226883838018*coeff[1]*fSkin[11]+0.8950343154210616*coeff[1]*fEdge[11]+7.585865087193007*coeff[1]*fSkin[8]+0.3513888460007659*coeff[1]*fEdge[8]; 

  boundSurf_incr[2] = -(0.9165151389911681*coeff[1]*fSkin[9])-1.3555441711725957*coeff[1]*fSkin[5]-1.35*coeff[1]*fSkin[2]-0.8660254037844386*fSkin[0]*coeff[1]; 
  boundSurf_incr[3] = -(0.916515138991168*coeff[1]*fSkin[12])-1.355544171172596*coeff[1]*fSkin[7]-1.35*coeff[1]*fSkin[3]-0.8660254037844386*coeff[1]*fSkin[1]; 
  boundSurf_incr[5] = -(3.5496478698597698*coeff[1]*fSkin[9])-5.25*coeff[1]*fSkin[5]-5.2285275173800105*coeff[1]*fSkin[2]-3.3541019662496847*fSkin[0]*coeff[1]; 
  boundSurf_incr[6] = -(0.916515138991168*coeff[1]*fSkin[14])-1.355544171172596*coeff[1]*fSkin[10]-1.35*coeff[1]*fSkin[6]-0.8660254037844387*coeff[1]*fSkin[4]; 
  boundSurf_incr[7] = -(3.5496478698597698*coeff[1]*fSkin[12])-5.25*coeff[1]*fSkin[7]-5.228527517380013*coeff[1]*fSkin[3]-3.3541019662496843*coeff[1]*fSkin[1]; 
  boundSurf_incr[9] = -(8.4*coeff[1]*fSkin[9])-12.423767544509195*coeff[1]*fSkin[5]-12.372954376380765*coeff[1]*fSkin[2]-7.937253933193772*fSkin[0]*coeff[1]; 
  boundSurf_incr[10] = -(3.5496478698597698*coeff[1]*fSkin[14])-5.25*coeff[1]*fSkin[10]-5.228527517380013*coeff[1]*fSkin[6]-3.3541019662496847*coeff[1]*fSkin[4]; 
  boundSurf_incr[11] = -(0.916515138991168*coeff[1]*fSkin[15])-1.355544171172596*coeff[1]*fSkin[13]-1.35*coeff[1]*fSkin[11]-0.8660254037844386*coeff[1]*fSkin[8]; 
  boundSurf_incr[12] = -(8.4*coeff[1]*fSkin[12])-12.423767544509193*coeff[1]*fSkin[7]-12.372954376380768*coeff[1]*fSkin[3]-7.937253933193771*coeff[1]*fSkin[1]; 
  boundSurf_incr[13] = -(3.5496478698597698*coeff[1]*fSkin[15])-5.25*coeff[1]*fSkin[13]-5.228527517380012*coeff[1]*fSkin[11]-3.354101966249684*coeff[1]*fSkin[8]; 
  boundSurf_incr[14] = -(8.4*coeff[1]*fSkin[14])-12.423767544509193*coeff[1]*fSkin[10]-12.37295437638077*coeff[1]*fSkin[6]-7.937253933193772*coeff[1]*fSkin[4]; 
  boundSurf_incr[15] = -(8.4*coeff[1]*fSkin[15])-12.423767544509195*coeff[1]*fSkin[13]-12.372954376380768*coeff[1]*fSkin[11]-7.937253933193772*coeff[1]*fSkin[8]; 

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

  return 0.;
}

