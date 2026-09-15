#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_1x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[1]; 

  double alphaDrag[6]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (1.4142135623730951*nuUSum[0]-1.4142135623730951*nuSum[0]*w[1])*rdvx2; 
  alphaDrag[1] = (1.4142135623730951*nuUSum[1]-1.4142135623730951*nuSum[1]*w[1])*rdvx2; 
  alphaDrag[2] = -(0.408248290463863*nuSum[0]*dxv[1]*rdvx2); 
  alphaDrag[3] = -(0.408248290463863*dxv[1]*nuSum[1]*rdvx2); 

  out[2] += 0.8660254037844386*(alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[4] += 1.7320508075688772*(alphaDrag[3]*f[5]+alphaDrag[2]*f[4])+1.9364916731037085*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[5] += 1.7320508075688772*(alphaDrag[2]*f[5]+alphaDrag[3]*f[4])+1.9364916731037085*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 

  return fabs(1.25*alphaDrag[0]); 

} 
