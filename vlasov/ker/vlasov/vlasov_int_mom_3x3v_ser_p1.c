#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_int_five_moments_3x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*dxv[5]*0.015625; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 8.0*f[0]*volFact; 
  out[1] += volFact*(8.0*f[0]*wx1+2.309401076758503*f[4]*dv1); 
  out[2] += volFact*(8.0*f[0]*wx2+2.309401076758503*f[5]*dv2); 
  out[3] += volFact*(8.0*f[0]*wx3+2.309401076758503*f[6]*dv3); 
  out[4] += volFact*(8.0*f[0]*wx3_sq+4.618802153517007*f[6]*dv3*wx3+8.0*f[0]*wx2_sq+4.618802153517007*f[5]*dv2*wx2+8.0*f[0]*wx1_sq+4.618802153517007*f[4]*dv1*wx1+0.5962847939999438*f[24]*dv3_sq+0.6666666666666666*f[0]*dv3_sq+0.5962847939999438*f[23]*dv2_sq+0.6666666666666666*f[0]*dv2_sq+0.5962847939999438*f[22]*dv1_sq+0.6666666666666666*f[0]*dv1_sq); 
} 
