#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_int_five_moments_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*0.0625; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += volFact*(4.0*f[0]*wx1+1.154700538379252*f[3]*dv1); 
  out[2] += volFact*(4.0*f[0]*wx2+1.154700538379252*f[4]*dv2); 
  out[3] += volFact*(4.0*f[0]*wx2_sq+2.309401076758503*f[4]*dv2*wx2+4.0*f[0]*wx1_sq+2.309401076758503*f[3]*dv1*wx1+0.2981423969999719*f[14]*dv2_sq+0.3333333333333333*f[0]*dv2_sq+0.2981423969999719*f[13]*dv1_sq+0.3333333333333333*f[0]*dv1_sq); 
} 
