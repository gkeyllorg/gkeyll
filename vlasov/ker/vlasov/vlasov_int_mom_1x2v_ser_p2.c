#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_int_five_moments_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[2]*dv1); 
  out[2] += volFact*(2.828427124746191*f[0]*wx2+0.8164965809277261*f[3]*dv2); 
  out[3] += volFact*(2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[3]*dv2*wx2+2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[2]*dv1*wx1+0.210818510677892*f[9]*dv2_sq+0.2357022603955158*f[0]*dv2_sq+0.210818510677892*f[8]*dv1_sq+0.2357022603955158*f[0]*dv1_sq); 
} 
