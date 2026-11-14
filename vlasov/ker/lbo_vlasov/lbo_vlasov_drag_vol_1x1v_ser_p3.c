#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_1x1v_ser_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[1]; 

  double alphaDrag[12]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (1.4142135623730951*nuUSum[0]-1.4142135623730951*nuSum[0]*w[1])*rdvx2; 
  alphaDrag[1] = (1.4142135623730951*nuUSum[1]-1.4142135623730951*nuSum[1]*w[1])*rdvx2; 
  alphaDrag[2] = -(0.408248290463863*nuSum[0]*dxv[1]*rdvx2); 
  alphaDrag[3] = -(0.408248290463863*dxv[1]*nuSum[1]*rdvx2); 
  alphaDrag[4] = (1.4142135623730951*nuUSum[2]-1.4142135623730951*w[1]*nuSum[2])*rdvx2; 
  alphaDrag[6] = -(0.408248290463863*dxv[1]*nuSum[2]*rdvx2); 
  alphaDrag[8] = (1.4142135623730951*nuUSum[3]-1.4142135623730951*w[1]*nuSum[3])*rdvx2; 
  alphaDrag[10] = -(0.408248290463863*dxv[1]*nuSum[3]*rdvx2); 

  out[2] += 0.8660254037844386*(alphaDrag[10]*f[10]+alphaDrag[8]*f[8]+alphaDrag[6]*f[6]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.7606388292556648*(alphaDrag[6]*f[10]+f[6]*alphaDrag[10]+alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+0.7745966692414833*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 1.9364916731037085*(alphaDrag[8]*f[10]+f[8]*alphaDrag[10])+1.7320508075688772*alphaDrag[3]*f[7]+1.9364916731037085*(alphaDrag[4]*f[6]+f[4]*alphaDrag[6])+1.7320508075688772*alphaDrag[2]*f[5]+1.9364916731037085*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[6] += 0.5163977794943223*alphaDrag[10]*f[10]+0.7606388292556648*(alphaDrag[3]*f[10]+f[3]*alphaDrag[10])+0.5163977794943223*alphaDrag[8]*f[8]+0.7606388292556648*(alphaDrag[1]*f[8]+f[1]*alphaDrag[8])+0.5532833351724881*alphaDrag[6]*f[6]+0.8660254037844386*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6])+0.5532833351724881*alphaDrag[4]*f[4]+0.8660254037844386*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4])+0.7745966692414833*(alphaDrag[3]*f[3]+alphaDrag[1]*f[1]); 
  out[7] += 1.7008401285415224*(alphaDrag[4]*f[10]+f[4]*alphaDrag[10]+alphaDrag[6]*f[8]+f[6]*alphaDrag[8])+1.5491933384829666*alphaDrag[6]*f[7]+1.7320508075688772*(alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[3]*(f[5]+f[4])+f[3]*alphaDrag[4])+1.9364916731037085*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[9] += 2.5980762113533156*alphaDrag[3]*f[11]+3.968626966596886*alphaDrag[10]*f[10]+2.5980762113533156*alphaDrag[2]*f[9]+1.3228756555322954*alphaDrag[8]*f[8]+2.9580398915498085*alphaDrag[1]*f[7]+3.968626966596886*alphaDrag[6]*f[6]+2.9580398915498085*alphaDrag[0]*f[5]+1.3228756555322954*alphaDrag[4]*f[4]+3.968626966596886*(alphaDrag[3]*f[3]+alphaDrag[2]*f[2])+1.3228756555322954*(alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[10] += (0.5163977794943223*alphaDrag[6]+0.8660254037844386*alphaDrag[2])*f[10]+(0.5163977794943223*f[6]+0.8660254037844386*f[2])*alphaDrag[10]+(0.5163977794943223*alphaDrag[4]+0.8660254037844386*alphaDrag[0])*f[8]+(0.5163977794943223*f[4]+0.8660254037844386*f[0])*alphaDrag[8]+0.7606388292556648*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]); 
  out[11] += (2.3237900077244498*alphaDrag[6]+2.5980762113533156*alphaDrag[2])*f[11]+3.485685011586674*(alphaDrag[6]*f[10]+f[6]*alphaDrag[10])+2.5980762113533156*alphaDrag[3]*f[9]+1.1618950038622249*(alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+(2.6457513110645907*alphaDrag[4]+2.9580398915498085*alphaDrag[0])*f[7]+3.5496478698597698*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6])+2.9580398915498085*alphaDrag[1]*f[5]+1.1832159566199232*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+3.968626966596886*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3])+1.3228756555322954*(alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 

  return fabs(1.75*alphaDrag[0]-1.956559480312316*alphaDrag[4]); 

} 
