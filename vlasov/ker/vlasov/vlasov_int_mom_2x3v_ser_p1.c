#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_int_five_moments_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*0.03125; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += volFact*(5.656854249492382*f[0]*wx1+1.632993161855453*f[3]*dv1); 
  out[2] += volFact*(5.656854249492382*f[0]*wx2+1.632993161855453*f[4]*dv2); 
  out[3] += volFact*(5.656854249492382*f[0]*wx3+1.632993161855453*f[5]*dv3); 
  out[4] += volFact*(5.656854249492382*f[0]*wx3_sq+3.265986323710906*f[5]*dv3*wx3+5.656854249492382*f[0]*wx2_sq+3.265986323710906*f[4]*dv2*wx2+5.656854249492382*f[0]*wx1_sq+3.265986323710906*f[3]*dv1*wx1+0.421637021355784*f[18]*dv3_sq+0.4714045207910317*f[0]*dv3_sq+0.421637021355784*f[17]*dv2_sq+0.4714045207910317*f[0]*dv2_sq+0.421637021355784*f[16]*dv1_sq+0.4714045207910317*f[0]*dv1_sq); 
} 
