#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_vol_3x3v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[6]: Cell-center coordinates. 
  // dxv[6]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuVtSqSum = &nuPrimMomsSum[24];

  const double rdvxSq4 = 4.0/(dxv[3]*dxv[3]); 
  const double rdvySq4 = 4.0/(dxv[4]*dxv[4]); 
  const double rdvzSq4 = 4.0/(dxv[5]*dxv[5]); 

  return fabs(3.181980515339463*nuVtSqSum[0]*rdvxSq4)+fabs(3.181980515339463*nuVtSqSum[0]*rdvySq4)+fabs(3.181980515339463*nuVtSqSum[0]*rdvzSq4); 

} 
