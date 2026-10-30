#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_diff_vol_2x2v_tensor_p1(const double *w, const double *dxv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds squared times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuVtSqSum = &nuPrimMomsSum[8];

  const double rdvxSq4 = 4.0/(dxv[2]*dxv[2]); 
  const double rdvySq4 = 4.0/(dxv[3]*dxv[3]); 

  return fabs(4.5*nuVtSqSum[0]*rdvxSq4)+fabs(4.5*nuVtSqSum[0]*rdvySq4); 

} 
