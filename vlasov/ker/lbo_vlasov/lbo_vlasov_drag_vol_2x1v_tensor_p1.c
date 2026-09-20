#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_drag_vol_2x1v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[2]; 

  double alphaDrag[12]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (1.4142135623730951*nuUSum[0]-1.4142135623730951*nuSum[0]*w[2])*rdvx2; 
  alphaDrag[1] = (1.4142135623730951*nuUSum[1]-1.4142135623730951*nuSum[1]*w[2])*rdvx2; 
  alphaDrag[2] = (1.4142135623730951*nuUSum[2]-1.4142135623730951*nuSum[2]*w[2])*rdvx2; 
  alphaDrag[3] = -(0.408248290463863*nuSum[0]*dxv[2]*rdvx2); 
  alphaDrag[4] = (1.4142135623730951*nuUSum[3]-1.4142135623730951*w[2]*nuSum[3])*rdvx2; 
  alphaDrag[5] = -(0.408248290463863*nuSum[1]*dxv[2]*rdvx2); 
  alphaDrag[6] = -(0.408248290463863*dxv[2]*nuSum[2]*rdvx2); 
  alphaDrag[8] = -(0.408248290463863*dxv[2]*nuSum[3]*rdvx2); 

  out[3] += 0.6123724356957944*(alphaDrag[8]*f[8]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[5] += 0.6123724356957944*(alphaDrag[6]*f[8]+f[6]*alphaDrag[8]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.6123724356957944*(alphaDrag[5]*f[8]+f[5]*alphaDrag[8]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[7] += 1.224744871391589*(alphaDrag[8]*f[11]+alphaDrag[6]*f[10]+alphaDrag[5]*f[9])+1.369306393762915*(alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+1.224744871391589*alphaDrag[3]*f[7]+1.369306393762915*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[8] += 0.6123724356957944*(alphaDrag[3]*f[8]+f[3]*alphaDrag[8]+alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[9] += 1.224744871391589*(alphaDrag[6]*f[11]+alphaDrag[8]*f[10]+alphaDrag[3]*f[9])+1.369306393762915*(alphaDrag[2]*f[8]+f[2]*alphaDrag[8])+1.224744871391589*alphaDrag[5]*f[7]+1.369306393762915*(alphaDrag[4]*f[6]+f[4]*alphaDrag[6]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[10] += 1.224744871391589*(alphaDrag[5]*f[11]+alphaDrag[3]*f[10]+alphaDrag[8]*f[9])+1.369306393762915*(alphaDrag[1]*f[8]+f[1]*alphaDrag[8])+1.224744871391589*alphaDrag[6]*f[7]+1.369306393762915*(alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[4]*f[5]+f[4]*alphaDrag[5]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[11] += 1.224744871391589*(alphaDrag[3]*f[11]+alphaDrag[5]*f[10]+alphaDrag[6]*f[9])+1.369306393762915*alphaDrag[0]*f[8]+1.224744871391589*f[7]*alphaDrag[8]+1.369306393762915*(f[0]*alphaDrag[8]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4]); 

  return fabs(0.8838834764831842*alphaDrag[0]); 

} 
