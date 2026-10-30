#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_1x1v_tensor_p3(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[1]; 

  double alphaDrag[16]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (1.4142135623730951*nuUSum[0]-1.4142135623730951*nuSum[0]*w[1])*rdvx2; 
  alphaDrag[1] = (1.4142135623730951*nuUSum[1]-1.4142135623730951*nuSum[1]*w[1])*rdvx2; 
  alphaDrag[2] = -(0.408248290463863*nuSum[0]*dxv[1]*rdvx2); 
  alphaDrag[3] = -(0.408248290463863*dxv[1]*nuSum[1]*rdvx2); 
  alphaDrag[4] = (1.4142135623730951*nuUSum[2]-1.4142135623730951*w[1]*nuSum[2])*rdvx2; 
  alphaDrag[6] = -(0.408248290463863*dxv[1]*nuSum[2]*rdvx2); 
  alphaDrag[8] = (1.4142135623730951*nuUSum[3]-1.4142135623730951*w[1]*nuSum[3])*rdvx2; 
  alphaDrag[11] = -(0.408248290463863*dxv[1]*nuSum[3]*rdvx2); 

  out[2] += 0.8660254037844386*(alphaDrag[11]*f[11]+alphaDrag[8]*f[8]+alphaDrag[6]*f[6]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.7606388292556648*(alphaDrag[6]*f[11]+f[6]*alphaDrag[11]+alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+0.7745966692414833*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 1.7320508075688772*alphaDrag[11]*f[13]+1.9364916731037085*(alphaDrag[8]*f[11]+f[8]*alphaDrag[11])+1.7320508075688772*(alphaDrag[6]*f[10]+alphaDrag[3]*f[7])+1.9364916731037085*(alphaDrag[4]*f[6]+f[4]*alphaDrag[6])+1.7320508075688772*alphaDrag[2]*f[5]+1.9364916731037085*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[6] += 0.5163977794943223*alphaDrag[11]*f[11]+0.7606388292556648*(alphaDrag[3]*f[11]+f[3]*alphaDrag[11])+0.5163977794943223*alphaDrag[8]*f[8]+0.7606388292556648*(alphaDrag[1]*f[8]+f[1]*alphaDrag[8])+0.5532833351724881*alphaDrag[6]*f[6]+0.8660254037844386*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6])+0.5532833351724881*alphaDrag[4]*f[4]+0.8660254037844386*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4])+0.7745966692414833*(alphaDrag[3]*f[3]+alphaDrag[1]*f[1]); 
  out[7] += 1.5212776585113297*alphaDrag[6]*f[13]+1.7008401285415224*alphaDrag[4]*f[11]+(1.5212776585113297*f[10]+1.7008401285415224*f[4])*alphaDrag[11]+1.5491933384829666*alphaDrag[3]*f[10]+1.7008401285415224*(alphaDrag[6]*f[8]+f[6]*alphaDrag[8])+1.5491933384829666*alphaDrag[6]*f[7]+1.7320508075688772*(alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[3]*(f[5]+f[4])+f[3]*alphaDrag[4])+1.9364916731037085*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[9] += 2.5980762113533156*(alphaDrag[11]*f[15]+alphaDrag[6]*f[14])+2.9580398915498085*alphaDrag[8]*f[13]+2.5980762113533156*alphaDrag[3]*f[12]+3.968626966596886*alphaDrag[11]*f[11]+2.9580398915498085*alphaDrag[4]*f[10]+2.5980762113533156*alphaDrag[2]*f[9]+1.3228756555322954*alphaDrag[8]*f[8]+2.9580398915498085*alphaDrag[1]*f[7]+3.968626966596886*alphaDrag[6]*f[6]+2.9580398915498085*alphaDrag[0]*f[5]+1.3228756555322954*alphaDrag[4]*f[4]+3.968626966596886*(alphaDrag[3]*f[3]+alphaDrag[2]*f[2])+1.3228756555322954*(alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[10] += (1.0327955589886446*alphaDrag[11]+1.5212776585113297*alphaDrag[3])*f[13]+(1.1547005383792517*alphaDrag[8]+1.7008401285415224*alphaDrag[1])*f[11]+(1.1547005383792517*f[8]+1.5212776585113297*f[7]+1.7008401285415224*f[1])*alphaDrag[11]+(1.1065666703449761*alphaDrag[6]+1.7320508075688772*alphaDrag[2])*f[10]+1.7008401285415224*(alphaDrag[3]*f[8]+f[3]*alphaDrag[8])+1.5491933384829666*alphaDrag[3]*f[7]+(1.2371791482634837*alphaDrag[4]+1.9364916731037085*alphaDrag[0])*f[6]+(1.7320508075688772*f[5]+1.2371791482634837*f[4])*alphaDrag[6]+1.9364916731037085*(f[0]*alphaDrag[6]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4])+1.7320508075688772*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[11] += (0.5163977794943223*alphaDrag[6]+0.8660254037844386*alphaDrag[2])*f[11]+(0.5163977794943223*f[6]+0.8660254037844386*f[2])*alphaDrag[11]+(0.5163977794943223*alphaDrag[4]+0.8660254037844386*alphaDrag[0])*f[8]+(0.5163977794943223*f[4]+0.8660254037844386*f[0])*alphaDrag[8]+0.7606388292556648*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]); 
  out[12] += 2.2819164877669946*alphaDrag[6]*f[15]+(2.2819164877669946*alphaDrag[11]+2.3237900077244498*alphaDrag[3])*f[14]+2.5980762113533156*alphaDrag[4]*f[13]+(2.3237900077244498*alphaDrag[6]+2.5980762113533156*alphaDrag[2])*f[12]+3.485685011586674*(alphaDrag[6]*f[11]+f[6]*alphaDrag[11])+(2.5980762113533156*alphaDrag[8]+2.6457513110645907*alphaDrag[1])*f[10]+2.5980762113533156*alphaDrag[3]*f[9]+1.1618950038622249*(alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+(2.6457513110645907*alphaDrag[4]+2.9580398915498085*alphaDrag[0])*f[7]+3.5496478698597698*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6])+2.9580398915498085*alphaDrag[1]*f[5]+1.1832159566199232*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+3.968626966596886*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3])+1.3228756555322954*(alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[13] += (1.0327955589886446*alphaDrag[6]+1.7320508075688772*alphaDrag[2])*f[13]+(1.1547005383792517*alphaDrag[4]+1.9364916731037085*alphaDrag[0])*f[11]+(1.0327955589886446*f[10]+1.7320508075688772*f[5]+1.1547005383792517*f[4]+1.9364916731037085*f[0])*alphaDrag[11]+1.5212776585113297*alphaDrag[3]*f[10]+(1.1547005383792517*alphaDrag[6]+1.9364916731037085*alphaDrag[2])*f[8]+(1.1547005383792517*f[6]+1.9364916731037085*f[2])*alphaDrag[8]+1.5212776585113297*alphaDrag[6]*f[7]+1.7008401285415224*(alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4]); 
  out[14] += (1.5491933384829666*alphaDrag[11]+2.2819164877669946*alphaDrag[3])*f[15]+(1.659850005517464*alphaDrag[6]+2.5980762113533156*alphaDrag[2])*f[14]+(1.7638342073763937*alphaDrag[8]+2.5980762113533156*alphaDrag[1])*f[13]+(2.2819164877669946*alphaDrag[11]+2.3237900077244498*alphaDrag[3])*f[12]+2.3664319132398464*alphaDrag[11]*f[11]+3.485685011586674*(alphaDrag[3]*f[11]+f[3]*alphaDrag[11])+(1.889822365046136*alphaDrag[4]+2.9580398915498085*alphaDrag[0])*f[10]+2.5980762113533156*alphaDrag[6]*f[9]+(0.7888106377466154*alphaDrag[8]+1.1618950038622249*alphaDrag[1])*f[8]+(2.5980762113533156*f[7]+1.1618950038622249*f[1])*alphaDrag[8]+2.6457513110645907*alphaDrag[1]*f[7]+2.53546276418555*alphaDrag[6]*f[6]+3.968626966596886*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6])+alphaDrag[4]*(2.9580398915498085*f[5]+0.8451542547285166*f[4])+1.3228756555322954*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4])+3.5496478698597698*alphaDrag[3]*f[3]+1.1832159566199232*alphaDrag[1]*f[1]; 
  out[15] += (1.5491933384829666*alphaDrag[6]+2.5980762113533156*alphaDrag[2])*f[15]+(1.5491933384829666*alphaDrag[11]+2.2819164877669946*alphaDrag[3])*f[14]+(1.7638342073763937*alphaDrag[4]+2.9580398915498085*alphaDrag[0])*f[13]+2.2819164877669946*alphaDrag[6]*f[12]+(2.3664319132398464*alphaDrag[6]+3.968626966596886*alphaDrag[2])*f[11]+(2.5980762113533156*f[9]+2.3664319132398464*f[6]+3.968626966596886*f[2])*alphaDrag[11]+(1.7638342073763937*alphaDrag[8]+2.5980762113533156*alphaDrag[1])*f[10]+(0.7888106377466154*alphaDrag[4]+1.3228756555322954*alphaDrag[0])*f[8]+(2.9580398915498085*f[5]+0.7888106377466154*f[4]+1.3228756555322954*f[0])*alphaDrag[8]+2.5980762113533156*alphaDrag[4]*f[7]+3.485685011586674*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6])+1.1618950038622249*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4]); 

  return fabs(1.75*alphaDrag[0]-1.956559480312316*alphaDrag[4]); 

} 
